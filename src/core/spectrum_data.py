#!/usr/bin/env python
#
# spectrum_data.py
#
# Containers for spectrum data
#
# Author A R Back 
#
# 31/01/2014 <ab571@sussex.ac.uk> : First revision
# 20/02/2014 <ab571@sussex.ac.uk> : Updated __init__, more info extracted 
#                                   from filename. Updated get_histogram,
#                                   can be obtained from root file.
# 29/04/2014 <ab571@sussex.ac.uk> : Now initialised with half life - for
#                                   limit setting
# 09/05/2014 <ab571@sussex.ac.uk> : Made scale_histogram into two separate
#                                   methods, scale_by_t_half and scale_by_
#                                   mass. Also added new methods, get_
#                                   number_nuclei and set_events_by_t_half/
#                                   _mass
# 21/05/2014 <ab571@sussex.ac.uk> : Extraction from filename more intuitive
#                                   also moved extraction to separate
#                                   method, but still by default for 
#                                   correct filenames, to allow for 
#                                   integration with production macros.
# 05/06/2014 <ab571@sussex.ac.uk> : Changed handling of default values. 
#                                   Added scale_histogram method, to 
#                                   simplify existing scale methods
# 06/06/2014 <ab571@sussex.ac.uk> : Refactored for storage of multiple 
#                                   labelled histograms in memory
###########################################################################
from ROOT import TH1D
from ROOT import TFile

import rat

import isotope
import constants
import defaults
import file_manips
import list_manips
import spectrum_utils
import roi_utils
import generators

import re
import inspect
import os
import unittest

class SpectrumData(object):
    """ Base class, provides containers for spectrum data. """
    def __init__(self, path, t_half=None):
        """ Initialise the container. 

        MajoRat will automatically set relevant parameters, when the 
        filename is of the form:
        "RAT[version]_[^N_Evts^]_[generator]_[(type)]_[isotope]_[(level)]_[*mode*]_[(E_lo)]_[(E_hi)].root"
        otherwise the set_parameters method should be called.
        ^only required if setting parameters from filename^
        *required for decay0 generator*
        (optional)

        ***Files MUST be RAT generated ROOT files (ntuples are fine)***
        
        To set parameters from the filename it must specify at least 
         * RAT[version]
         * [N_Evts]
         * [generator]
         * [isotope]
        Defaults will be applied for the other options.
        :param path: path to RAT generated ROOT file
        :type path: str
        :param t_half: half life - required only for double beta spectra
        :type t_half: float
        """
        self._path = path
        print self._path
        self._dir, self._filename = file_manips.split_path(path)
        self._name, self._ext = file_manips.split_ext(self._filename)
        self._options = None
        self._prefix = ""
        self._rat_release = None
        self._n_events = None
        self._generator = None
        self._spectral_index = None
        self._label = None
        options = self._name.split("_")
        if (list_manips.item_containing("RAT", options) != None):
            if (len(options[list_manips.item_containing("RAT", options):]) > 3):
                # Assume filename is in MajoRat format
                # call set_parameters_from_filename by default
                self.set_parameters_from_filename(options)
        else:
            self._options = options
        if (self._label != None): # has been set
            self._label += "-Truth"
        self._histograms = {}
        self._unscaled_histograms = {}
        self._events = None
        self._t_half = t_half
    def set_parameters_from_filename(self, options):
        for item in options[:list_manips.item_containing("RAT", options)]:
            self._prefix += item + "_"
        main_options = options[list_manips.item_containing("RAT", options):]
        assert len(main_options) >= 4, \
            ("SpectrumData.set_parameters_from_filename: ERROR can't set parameters\n"
             "To set parameters from filename, it must contain at least:"
             "  RAT[version]\n"
             "  [N_Evts]\n"
             "  [generator]\n"
             "  [isotope]")
        option = iter(main_options)
        self._rat_release = option.next().replace("-",".")
        n_events = option.next()
        if (n_events.find("prod") >= 0):
            self._n_events = None # Set as default for production
        elif n_events.find("k") >= 0:
            evts = int(n_events[:n_events.find("k")]) * 1000
            self._n_events = evts
        elif n_events.find("M") >= 0:
            evts = int(n_events[:n_events.find("M")]) * 1e6
            self._n_events = evts
        else:
            evts = int(n_events[:re.search("[0-9][A-Za-z]*$", n_events).start()+1])
            self._n_events = evts
        generator = option.next()
        if (generator == "solar"):
            isotope_name = option.next()
            try:
                type_ = option.next()
            except StopIteration:
                # Set type as "nue", which is default
                type_ = "nue"
            # Assume energy range as default
            self._generator = generators.Solar(isotope_name, type_)
            self._label = self._generator.get_isotope()
        elif (generator == "decaychain"):
            isotope_name = option.next()
            # Assume energy range as default
            self._generator = generators.DecayChain(isotope_name)
            self._label = self._generator.get_isotope()
        elif (generator == "decay0"):
            self.set_decay0_options(option)
        else:
            print "SpectrumData.set_parameters_from_filename: ERROR generator "\
                "not recognised."
            sys.exit(1)
    def set_decay0_options(self, option):
        """ Method to set decay0 specific options which are a bit more 
        tricky.
        
        :param option: iterator over list options, created from filename
        :type option: iter
        """
        temp_option = option.next()
        if (re.search("[a-zA-Z]+[0-9]+", temp_option)): # isotope
            isotope_name = temp_option
            if (isotope_name == "Te130"):
                type_ = "2beta"
            else:
                type_ = "backg"
        else:
            type_ = temp_option
            isotope_name = option.next()
        if (type_ == "backg"):
            self._generator = generators.Backg(isotope_name)
            self._label = self._generator.get_isotope()
        elif (type_ == "2beta"):
            if (isotope_name == "Te130"): # Assume we want SNO+ specific instance
                isotope_ = isotope.SNOPlusTe(isotope_name)
            else: # Create standard Isotope instance
                isotope_ = isotope.Isotope(isotope_name)
            temp_option_1 = option.next()
            try:
                temp_option_2 = option.next()
            except StopIteration:
                # temp_option_1 is last parameter --> *mode*
                mode = float(temp_option_1)
                self._generator = generators.TwoBeta(isotope_, mode)
            level = float(temp_option_1)
            mode = float(temp_option_2)
            e_lo = float(option.next().replace("-","."))
            e_hi = float(option.next().replace("-","."))
            self._generator = generators.TwoBeta(isotope_, mode, level,
                                                 type_, e_lo, e_hi)
            self._spectral_index = spectrum_utils.get_spectral_index\
            (self._generator.get_mode())
            self._label = spectrum_utils.get_label(self._spectral_index)
    def make_histogram(self, hist_label="default", 
                       append=False, 
                       always_remake=False,
                       bin_width="default",
                       hist_path="default"):
        """ Method to generate make histogram 
        
        :param hist_label: self._histograms key (and ROOT label)
        :type hist_label: str
        :param append: if True adds to saved histogram
        :type append: bool
        :param always_remake: if True always fills new blank histogram
        :type always_remake: bool
        :param bin_width: width (MeV) of histogram bins
        :type bin_width: float
        :param hist_path: specify alternative path to saved histograms
        :type hist_path: str
        """
        if (hist_label == "default"):
            hist_label = self._label
        try:
            assert (self._histograms.get(hist_label) == None), \
                "histogram with label " + hist_label + " already exists!"
        except AssertionError as detail:
            print "SpectrumData.make_histogram: ERROR", detail
        try:
            assert (not(append and always_remake)),\
                "append and always_remake options cannot both be True"
        except AssertionError as detail:
            print "SpectrumData.make_histogram: ERROR", detail
            raise
        if (bin_width == "default"):
            bin_width = defaults.spectrum.get("bin_width")
        if (self._prefix.find("hist") >= 0):
            try:
                input_file = TFile(self._path, "READ")
                histogram = input_file.Get(hist_label)
                assert (isinstance(histogram, TH1D)), \
                    "no object " + hist_label + " of type ROOT.TH1D, found "\
                    "in file"
                histogram.SetDirectory(0) # Stop ROOT memory managing
            except AssertionError as detail:
                print "SpectrumData.make_histogram:", detail
                print " --> cannot create new histogram from this file"
                raise
            input_file.Close()
        elif (always_remake):
            histogram = self.get_blank_histogram(hist_label, bin_width)
            histogram.SetDirectory(0) # Stop ROOT memory managing
        else:
            if (hist_path == "default"):
                hist_path = self.get_default_hist_path()
            else:
                directory, file_ = file_manips.split_path(hist_path)
                if (len(file_) <= 0): # Alternative directory only
                    hist_path = directory + \
                        file_manips.strip_path(self.get_default_hist_path())
            try:
                input_file = TFile(hist_path, "READ")
                histogram = input_file.Get(hist_label)
                assert (isinstance(histogram, TH1D)), \
                    "no object " + hist_label + " of type ROOT.TH1D, found "\
                    "in file"
                histogram.SetDirectory(0) # Stop ROOT memory managing
            except AssertionError as detail:
                print "SpectrumData.make_histogram:", detail
                print " --> making new histogram"
                histogram = self.get_blank_histogram(hist_label, bin_width)
                always_remake = True # In this instance fill from scratch
        self._histograms[hist_label] = histogram
        assert isinstance(self._histograms[hist_label], TH1D), \
            "SpectrumData.make_histogram: histogram is not a TH1D object"
        if (append != always_remake): # Exclusive or, can't have both True 
            self.fill_histogram(hist_label)
        self._histograms[hist_label].SetDirectory(0) # Stop ROOT memory managing
        if (self._unscaled_histograms.get(hist_label) == None):
            self._unscaled_histograms[hist_label] = \
                self._histograms[hist_label].Clone()
            self._unscaled_histograms[hist_label].SetDirectory(0)
    def get_blank_histogram(self, hist_label="default",
                            bin_width="default",
                            no_label=False):
        """ Method to return a blank spectrum histogram to be filled.
        
        :param hist_label: label 
        :param bin_width: supply a bin width (MeV)
        :type bin_width: float
        :returns: blank histogram
        :rtype: ROOT.TH1D
        """
        if (bin_width == "default"):
            bin_width = defaults.spectrum.get("bin_width")
        lower = self._generator.get_e_lo()
        upper = self._generator.get_e_hi()
        n_bins = int((upper-lower) / bin_width)
        if no_label:
            histogram = TH1D("default_hist", "default_hist", n_bins, lower, upper)
        else:
            histogram = TH1D(hist_label, hist_label, n_bins, lower, upper)
        return histogram
    def get_default_hist_path(self):
        """
        :returns: file "hist_ ..." filename, where (if it exists)
                  histograms would be saved.
        :rtype: str
        """
        hist_path = os.environ.get("MAJORAT_DATA") + "/hist_" + self._name\
            + self._ext
        return str(hist_path)
    def fill_histogram(self, hist_label="default"):
        """ Method to read DS, extract total KE for each event and fill 
        histogram.
        
        :param hist_label: histogram key in self._histograms
        :type hist_label: str
        """
        if (hist_label == "default"):
            hist_label = self._label
        histogram = self._histograms.get(hist_label)
        assert isinstance(histogram, TH1D), \
            "SpectrumData.fill_histogram: histogram is not a TH1D object"
        for ds, run in rat.dsreader(self._path):
            for event in range(0, ds.GetEVCount()):
                mc = ds.GetMC();
                truth_energy = 0
                for particle in range (0, mc.GetMCParticleCount()):
                    # Check they are electrons
                    if (mc.GetMCParticle(particle).GetPDGCode() == 11):
                        truth_energy += mc.GetMCParticle(particle).GetKE()
                histogram.Fill(truth_energy)
    def get_histogram(self, hist_label="default",
                      always_recreate=False,
                      bin_width="default"):
        """ 
        :param hist_label: histogram key in self._histograms
        :type hist_label: str
        :returns: histogram with label supplied (if None, calls make_histogram)
        :rtype: ROOT.TH1D
        """
        if (hist_label == "default"):
            hist_label = self._label
        if (bin_width == "default"):
            bin_width = defaults.spectrum.get("bin_width")
        if (self._histograms.get(hist_label) == None):
            append=False
            self.make_histogram(hist_label, append, always_recreate, bin_width)
        return self._histograms.get(hist_label)
    def get_unscaled_histogram(self, hist_label="default"):
        """ 
        :param hist_label: histogram key in self._unscaled_histograms
        :type hist_label: str
        :returns: clone of unscaled histogram with label supplied
        :rtype: ROOT.TH1D
        """
        if (hist_label == "default"):
            hist_label = self._label
        try:
            assert (self._unscaled_histograms.get(hist_label) != None), \
                "no unscaled histogram with label " + hist_label + " found!"
            # clone so we never alter the original unscaled version
            histogram = self._unscaled_histograms.get(hist_label).Clone()
        except AssertionError as detail:
            print "SpectrumData.get_unscaled_histogram: ERROR", detail
        return histogram
    def get_number_nuclei(self):
        """ Method to return the number of nuclei. The calculation of
        number of nuclei is dependent on the object type of self._isotope.
        This method checks the object type of self._isotope, then
        calculates the number of nuclei by the appropriate method. Returns
        number of nuclei.
        """
        if isinstance(self._generator.get_isotope(), isotope.SNOPlusTe):
            number_nuclei = self._generator.get_isotope().get_number_nuclei\
                (defaults.snoplus.get("fv_radius"))
        else:
            number_nuclei = self._generator.get_isotope().get_number_nuclei\
                (constants.isotope_mass.get\
                     (self._generator.get_isotope().get_name()))
        return number_nuclei
    def set_events_by_t_half(self, t_half="default", livetime="default"):
        """ Method that allows you to set self._events by half life. Assumes
        default number of nuclei.
        
        :param livetime: livetime in years
        :type livetime: float
        :param t_half: half life in years (default = self._t_half)
        :type t_half: float
        """
        if (livetime == "default"):
            livetime = defaults.livetime_general
        if (t_half == "default"):
            t_half = self._t_half
        assert (self._t_half != None),\
            "SpectrumData.set_events_by_t_half: ERROR t_half = None"
        self._events = self._generator.get_isotope().get_decays_from_t_half\
            (t_half, livetime)
        print "SpectrumData.set_events_by_t_half: number of events set at",\
            self._events
    def set_events_by_mass(self, effective_mass, livetime="default"):
        """ Method that allows you to set self._events by effective mass.
        Assumes default number of nuclei.

        :param effective_mass: effective double beta mass in eV (no default)
        :type effective_mass: float
        :param livetime: livetime in years
        :type livetime: float
        """
        if (livetime == "default"):
            livetime = defaults.livetime_general
        self._events = self._generator.get_isotope().get_decays_from_mass\
            (effective_mass, livetime)
        print "SpectrumData.set_events_by_mass: number of events set at",\
            self._events
    def scale_histogram(self, scaling_factor,
                        hist_label="default",
                        e_lo="default", 
                        e_hi="default"):
        """ Method for scaling any histogram by a scaling factor.
        
        :param scaling_factor: proposed integral of histogram after scaling
        :type scaling_factor: float
        :param hist_label: histogram key in self._histograms
        :type hist_label: str
        :param e_lo: lower energy integration limit
        :type e_lo: float
        :param e_hi: upper energy integration limit
        :type e_hi: float
        """
        full_spectrum = False
        if (hist_label == "default"):
            hist_label = self._label
        histogram = self.get_histogram(hist_label)
        if (e_lo == "default") or (e_hi == "default"):
            full_spectrum = True
        n_bins = histogram.GetNbinsX()
        # work out mean bin width
        bin_width_sum = 0.0
        for bin_ in range(0, n_bins):
            bin_width_sum += histogram.GetBinWidth(bin_)
        bin_width = bin_width_sum / n_bins
        if not(full_spectrum):
            # set integration limits
            from_bin = int(e_lo/bin_width)
            to_bin = int(e_hi/bin_width)
        try:
            if (full_spectrum):
                histogram.Scale(scaling_factor/histogram.Integral())
            else:
                histogram.Scale(scaling_factor/histogram.Integral\
                                    (from_bin, to_bin))
        except ZeroDivisionError as detail:
            print "SpectrumData.scale_histogram: WARNING", detail
            print " --> reverting to unscaled histogram then re-scaling"
            histogram = self.get_unscaled_histograms()
            if (full_spectrum):
                histogram.Scale(scaling_factor/histogram.Integral())
            else:
                histogram.Scale(scaling_factor/histogram.Integral\
                                    (from_bin, to_bin))
        print "SpectrumData.scale_histogram: histogram now contains:"
        if (full_spectrum):
            print " -->", histogram.Integral(),"events"
        else:
            print " -->", histogram.Integral(from_bin, to_bin),"events in range"
    def scale_by_t_half(self, t_half="default",
                        hist_label="default",
                        livetime="default"):
        """ Scaling is not done automatically in the constructor. This
        method scales the histogram according to the number of events,
        calculated from the half life.
        
        :param t_half: half life in years (default = self._t_half)
        :type t_half: float
        :param hist_label: label of histogram to scale
        :type hist_label: str
        :param livetime: livetime in years
        :type livetime: float
        """
        if (t_half == "default"):
            t_half = self._t_half
        if (hist_label == "default"):
            hist_label = self._label
        if (livetime == "default"):
            livetime = defaults.livetime_general
        self.set_events_by_t_half(t_half, livetime)
        self.scale_histogram(self._events, hist_label)
    def scale_by_mass(self, effective_mass,
                      hist_label="default",
                      livetime="default"):
        """ Scaling is not done automatically in the constructor. This
        method scales the histogram according to the number of events,
        calculated from the effective mass.
        
        :param effective_mass: effective double beta mass in eV (no default)
        :type effective_mass: float
        :param hist_label: label of histogram to scale
        :type hist_label: str
        :param livetime: livetime in years (default = 1.0)
        :type livetime: float
        """
        if (hist_label == "default"):
            hist_label = self._label
        if (livetime == "default"):
            livetime = defaults.livetime_general
        self.set_events_by_mass(effective_mass, livetime)
        self.scale_histogram(self._events, hist_label)
    def scale_by_events(self, events,
                        hist_label="default",
                        e_lo="default",
                        e_hi="default"):
        """ Scaling is not done automatically in the constructor. This
        method scales the histogram according to the number of events
        supplied.

        :param events: number of events in full spectrum
        :type events: float
        :param hist_label: label of histogram to scale
        :type hist_label: str
        :param e_lo: lower energy integration limit
        :type e_lo: float
        :param e_hi: upper energy integration limit
        :type e_hi: float
        """
        if (hist_label == "default"):
            hist_label = self._label
        self.scale_histogram(events, hist_label, e_lo, e_hi)
    def scale_by_roi_events(self, events, roi,
                            hist_label="default"):
        """ Scaling is not done automatically in the constructor. This
        method scales the histogram according to the number of events
        supplied, in the ROI supplied

        :param events: number of events in full spectrum
        :type events: float
        :param roi: region of interest to consider
        :type roi: roi_utils.RoIUtils
        :param hist_label: label of histogram to scale
        :type hist_label: str
        """
        assert (isinstance(roi, roi_utils.RoIUtils)), \
            "SpectrumData.scale_by_roi_events: ERROR, no valid "\
            "roi_utils.RoIUtils object"
        if (hist_label == "default"):
            hist_label = self._label
        e_lo = roi.get_lower_limit()
        e_hi = roi.get_upper_limit()
        self.scale_histogram(events, hist_label, e_lo, e_hi)

###########################################################################
if __name__ == "__main__":
    from ROOT import THStack
    from ROOT import TLegend
    from ROOT import TCanvas

    import argparse

    parser = argparse.ArgumentParser\
        (description=("Generate energy spectrum for the single (RAT-generated)"
                      " Root file specified"))
    parser.add_argument("root_file", help="path to RAT-generated Root file")
    args = parser.parse_args()
    print args

    t_half = 5.0e25 # years # KLZ limit from Gando et al.

    spectrum = SpectrumData(args.root_file, t_half)
    spectrum.scale_by_t_half(0.0)
    print spectrum._events
    spectrum.scale_by_mass(0.0)
    print spectrum._events
    histogram = spectrum.get_histogram()
    histogram.Draw()
    raw_input("RETURN to exit")
