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
###########################################################################
from ROOT import TH1D
from ROOT import TFile

import rat

import isotope
import constants
import file_manips
import list_manips
import spectrum_utils
import generators

import re
import inspect

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
        self._dir, self._filename = file_manips.cut_path(path)
        self._name = file_manips.strip_ext(self._filename)
        self._prefix = None
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
        self._histogram = None
        self._unscaled_histogram = None
        self._events = None
        self._t_half = t_half
    def set_parameters_from_filename(self, options):
        self._prefix = ""
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
        if n_events.find("k") >= 0:
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
            self._generator = generators.TwoBeta(isotope_, mode, type_,
                                                level, e_lo, e_hi)
            self._spectral_index = spectrum_utils.get_spectral_index\
            (self._generator.get_mode())
            self._label = spectrum_utils.get_label(self._spectral_index)
    def get_histogram(self):
        """ Accesses the DS, reads particle MC Truth energies and sums. 
        Histograms summed KE to produce an energy spectrum. Returns 
        histogram.
        """
        if (self._histogram == None):
            if (self._prefix.find("hist") >= 0):
                input_file = TFile(self._path, "READ")
                input_file.ls()
                self._histogram = input_file.Get(self._label+"-Truth")
                print "SpectrumData.get_histogram: found saved histogram " \
                    + self._label+"-Truth"
                self._histogram.SetDirectory(0)
                input_file.Close()
            else:
                self._histogram = TH1D(self._label+"-Truth", self._label+"-Truth",
                                       30, 0.0, 3.0)
                for ds, run in rat.dsreader(self._path):
                    for iEV in range(0, ds.GetEVCount()):
                        mc = ds.GetMC();
                        KE = 0
                        for particle in range (0, mc.GetMCParticleCount()):
                            # Check they are electrons
                            if (mc.GetMCParticle(particle).GetPDGCode() == 11):
                                KE += mc.GetMCParticle(particle).GetKE()
                        self._histogram.Fill(KE)
            if (self._unscaled_histogram == None):
                self._unscaled_histogram = self._histogram.Clone()
        else:
            print ("SpectrumData.get_histogram: SpectrumData._histogram object "
                   "already exists\n --> re-using!")
        return self._histogram
    def get_number_nuclei(self):
        """ Method to return the number of nuclei. The calculation of
        number of nuclei is dependent on the object type of self._isotope.
        This method checks the object type of self._isotope, then
        calculates the number of nuclei by the appropriate method. Returns
        number of nuclei.
        """
        if isinstance(self._generator.get_isotope(), isotope.SNOPlusTe):
            number_nuclei = self._generator.get_isotope().get_number_nuclei\
                (constants.snoplus.get("fv_radius"))
        else:
            number_nuclei = self._generator.get_isotope().get_number_nuclei\
                (constants.isotope_mass.get\
                     (self._generator.get_isotope().get_name()))
        return number_nuclei
    def set_events_by_t_half(self, t_half=None):
        """ Method that allows you to set self._events.  By default
        self._t_half is used unless an alternative half life is supplied.
        """
        if (t_half == None):
            t_half = self._t_half
        assert (self._t_half != None),\
            "SpectrumData.set_events_by_t_half: ERROR t_half = None"
        number_nuclei = self.get_number_nuclei()
        self._events = self._generator.get_isotope().get_decays_from_t_half\
            (t_half, number_nuclei)
        print "SpectrumData.set_events_by_t_half: number of events set at",\
            self._events
    def set_events_by_mass(self, effective_mass):
        """ Method that allows you to set self._events. By default 
        self._t_half is used unless an alternative half life is supplied.
        """
        number_nuclei = self.get_number_nuclei()
        self._events = self._generator.get_isotope().get_decays_from_mass\
            (effective_mass, number_nuclei)
        print "SpectrumData.set_events_by_mass: number of events set at",\
            self._events
    def scale_by_t_half(self, t_half=None):
        """ Scaling is not done automatically in the constructor. This
        method scales the histogram according to the number of events,
        which is calculated from the half life. By default self._t_half
        is used unless an alternative half life is supplied. 
        """
        if (t_half == None):
            t_half = self._t_half
        if (self._histogram == None):
            self._histogram = self.get_histogram()
        self.set_events_by_t_half(t_half)
        try:
            self._histogram.Scale(self._events/self._histogram.Integral())
        except ZeroDivisionError as detail:
            print "SpectrumData.scale_by_t_half: WARNING", detail
            self._histogram = self._unscaled_histogram.Clone()
            self._histogram.Scale(self._events/self._histogram.Integral())
            print " --> reverting to unscaled histogram then re-scaling"
        print "SpectrumData.scale_by_t_half: histogram integral is",\
            self._histogram.Integral()
    def scale_by_mass(self, effective_mass):
        """ Scaling is not done automatically in the constructor. This
        method scales the histogram according to the number of events,
        which is calculated from the half life. By default self._t_half
        is used unless an alternative half life is supplied. 
        """
        if (self._histogram == None):
            self._histogram = self.get_histogram()
        self.set_events_by_mass(effective_mass)
        try:
            self._histogram.Scale(self._events/self._histogram.Integral())
        except ZeroDivisionError as detail:
            print "SpectrumData.scale_by_mass: WARNING", detail
            self._histogram = self._unscaled_histogram.Clone()
            self._histogram.Scale(self._events/self._histogram.Integral())
            print " --> reverting to unscaled histogram then re-scaling"
        print "SpectrumData.scale_by_mass: histogram integral is",\
            self._histogram.Integral()

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
    print args.root_file

    t_half = 5.0e25 # years # KLZ limit from Gando et al.

    spectrum = SpectrumData(args.root_file, t_half)
    spectrum.scale_by_t_half(0.0)
    print spectrum._events
    spectrum.scale_by_mass(0.0)
    print spectrum._events
    histogram = spectrum.get_histogram()
    histogram.Draw()
    raw_input("RETURN to exit")
