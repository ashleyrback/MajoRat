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
###############################################################################
from ROOT import TH1D

import rat

import isotope
import constants
import file_manips
import list_manips
import spectrum_utils

import re
import inspect

class SpectrumData(object):
    """ Base class, provides containers for spectrum data. """
    def __init__(self, path, t_half):
        """ Initialise the container. Parameters are obtained from the
        filename in the path supplied. Files are assumed to be in the 
        form:
        "RAT[version]_[N_Evts]_[generator]_[type]_[isotope]_[level]_[mode]_[E_lo]_[E_hi].root"
        They do not have to specify all options, but at the very least
        must be root files and specify RAT[version], [N_Evts], 
        [generator] and [mode], defaults will be applied for the other
        options.
        The class should be initialised with a half-life as well,
        instead of getting it from the spectral index, half life 
        dictionary.
        """
        self._path = path
        self._dir, self._filename = file_manips.cut_path(path)
        self._name = file_manips.strip_ext(self._filename)
        options = self._name.split("_")

        self._prefix = ""
        for item in options[:list_manips.item_containing("RAT", options)]:
            self._prefix += item

        main_options = options[list_manips.item_containing("RAT", options):]
        assert len(main_options) >= 4, ("SpectrumData.__init__: supplied"
                                        "filename must contain at least:\n"
                                        "  RAT[version]\n"
                                        "  [N_Evts]\n"
                                        "  [generator]\n"
                                        "  [mode]")

        i_opt = iter(main_options)
        option = i_opt.next()
        self._rat_release = option.replace("-",".")

        option = i_opt.next()
        if option.find("k") >= 0:
            evts = int(option[:option.find("k")]) * 1000
            self._N_evts = evts
        elif option.find("M") >= 0:
            evts = int(option[:option.find("M")]) * 1e6
            self._N_evts = evts
        else:
            evts = int(option[:re.search("[0-9][A-Za-z]*$", option).start()+1])
            self._N_evts = evts

        option = i_opt.next()
        self._generator = option

        self._options_set = False
        if self._generator == "decay0":
            option = i_opt.next()
            if re.search("[a-zA-Z]+", option):
                self._type = option
                option = i_opt.next()
                if option == "Te130": # Assume we want SNO+ specific instance
                    self._isotope = isotope.SNOPlusTe(option)
                else: # Create standard Isotope instance
                    self._isotope = isotope.Isotope(option)
                option = i_opt.next()
                self._level = int(option)
                option = i_opt.next()
                self._mode = int(option)
                option = i_opt.next()
                if re.search("[0-9]+\.*[0-9]*", option):
                    option = option.replace("-", ".")
                    self._E_lo = float(option)
                    option = i_opt.next()
                if re.search("[0-9]+\.*[0-9]*", option):
                    option = option.replace("-", ".")
                    self._E_hi = float(option)
                    self._options_set = True
            elif re.search("[0-9]+", option):
                self._type = "2beta" # default
                self._isotope = isotope.SNOPlusTe("Te130") #default
                self._level = 0 #default
                self._mode = int(option)
                option = i_opt.next()
                if re.search("[0-9]+\.*[0-9]*", option):
                    option = option.replace("-", ".")
                    self._E_lo = float(option)
                    option = i_opt.next()
                if re.search("[0-9]+\.*[0-9]*", option):
                    option = option.replace("-", ".")
                    self._E_hi = float(option)
                    self._options_set = True
            assert (self._options_set == True), ("SpectrumData.__init__: unable"
                                                 " to create SpectrumData "
                                                 "object,\n --> filename has "
                                                 "incorrect fortmat.")
        
        self._spectral_index = spectrum_utils.get_spectral_index(self._mode)
        self._label = spectrum_utils.get_label(self._spectral_index)
        self._histogram = None
        self._events = None
        self._t_half = t_half
    def get_histogram(self):
        """ Accesses the DS, reads particle MC Truth energies and sums. 
        Histograms summed KE to produce an energy spectrum. Returns 
        histogram.
        """
        if self._histogram == None:
            if (self._prefix.find("hist_") >= 0):
                input_file = TFile(self._path)
                self._histogram = input_file.Get(self._label+"-Truth")
                print "hist = ", self._histogram
                self._histogram.SetDirectory(0)
                input_file.close()
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
        else:
            print ("SpectrumData.get_histogram: SpectrumData._histogram object"
                   "already exists\n --> re-using!")
        return self._histogram
    def scale_histogram(self):
        """ Scaling is not done automatically in the constructor. This
        method scales the histogram according to the number of events,
        which is calculated from the half life.
        """
        if self._events == None:
            if isinstance(self._isotope, isotope.SNOPlusTe):
                self._events = self._isotope.get_number_decays\
                    (self._t_half, self._isotope.get_number_nuclei\
                         (constants.snoplus.get("fv_radius")))
            else:
                self._events = self._isotope.get_number_decays\
                    (self._t_half, self._isotope.get_number_nuclei\
                         (constants.isotope_mass.get(self._name)))
        else:
            if self._histogram == None:
                self._histogram == get_histogram
            else:
                self._histogram.Scale(self._events/self._histogram.Integral())
            

################################################################################
if __name__ == "__main__":
    from ROOT import THStack
    from ROOT import TLegend
    from ROOT import TCanvas

    import argparse

    parser = argparse.ArgumentParser(description=("Generate energy spectrum for"
                                                  " the single (RAT-generated) "
                                                  "Root file specified"))
    parser.add_argument("root_file", help="path to RAT-generated Root file")
    args = parser.parse_args()
    print args
    print args.root_file

    t_half = 5.0e25 # years # KLZ limit from Gando et al.

    spectrum = SpectrumData(args.root_file, t_half)
    spectrum.scale_histogram()
    histogram = spectrum.get_histogram()
    histogram.Draw()
    raw_input("RETURN to exit")


