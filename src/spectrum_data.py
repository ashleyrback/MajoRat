#!/usr/bin/env python
#
# spectrum_data.py
#
# Containers for spectrum data
#
# Author A R Back - 31/01/2014 <ab571@sussex.ac.uk> : First revision
#        A R Back - 20/02/2014 <ab571@sussex.ac.uk> : Updated __init__, more
#                                                     info extracted from
#                                                     filename. Updated 
#                                                     get_histogram, can be 
#                                                     obtained from root file.   
###############################################################################
from ROOT import TH1D

import rat

import file_manips
import list_manips
from detector_parameters import mass_Te_in_FV
from detector_parameters import number_Te
from detector_parameters import n_decays

import re
import inspect

mode_to_index = {1:0, 4:5, 5:1, 6:2, 7:3, 8:7}

spectrum_label = {0 : "0#nu#beta#beta",
                  1 : "0#nu#beta#beta#chi^{0} (n=1)",
                  2 : "0#nu#beta#beta#chi^{0} 'bulk' (n=2)",
                  3 : "0#nu#beta#beta#chi^{0}(#chi^{0}) (n=3)",
                  5 : "2#nu#beta#beta (n=5)",
                  7 : "0#nu#beta#beta#chi^{0}#chi^{0} (n=7)"}

number_Te = number_Te(mass_Te_in_FV(0.003, 4.0))

half_lives = {0 : 5.0e25,
              1 : 3.0e24,
              2 : 1.0e24,
              3 : 5.0e23,
              5 : 1.0e21,
              7 : 1.0e22}

class SpectrumData(object):
    """ Base class, provides containers for spectrum data. """
    def __init__(self, path):
        """ Initialise the container. Parameters are obtained from the
        filename in the path supplied. Files are assumed to be in the 
        form:
        "RAT[version]_[N_Evts]_[generator]_[type]_[isotope]_[level]_[mode]_[E_lo]_[E_hi].root"
        They do not have to specify all options, but at the very least
        must be root files and specify RAT[version], [N_Evts], 
        [generator] and [mode], defaults will be applied for the other
        options.
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
                self._isotope = option
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
                self._isotope = "Te130" #default
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
        
        self._spectral_index = mode_to_index.get(self._mode)
        self._label = spectrum_label.get(self._spectral_index)
        self._histogram = None
        self._t_half = half_lives.get(self._spectral_index)
        self._events = n_decays(self._t_half, number_Te)

    def get_histogram(self):
        """ Accesses the DS, reads particle MC Truth energies and sums. 
        Histograms summed KE to produce an energy spectrum. Returns histogram.
        """
        if self._histogram == None:
            if (self._prefix.find("hist_") >= 0):
                input_file = TFile(self._path)
                self._histogram = input_file.Get(self._label+"-Truth")
                print "hist = ", self._histogram
                self._histogram.SetDirectory(0)
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
                self._histogram.Scale(self._events/self._histogram.Integral())
        else:
            print ("SpectrumData.get_histogram: SpectrumData._histogram object"
                   "already exists\n --> re-using!")
        return self._histogram

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

    spectrum = SpectrumData(args.root_file)
    histogram = spectrum.get_histogram()
    histogram.Draw()
    raw_input("RETURN to exit")


