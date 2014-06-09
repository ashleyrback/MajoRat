#!/usr/bin/env python
#
# production.py
#
# Produce histogram of spectrum from SNO+ production ntuples
#
# Author A R Back 
#
# 15/05/2014 <ab571@sussex.ac.uk> : First revision
# 06/06/2014 <ab571@sussex.ac.uk> : Refactored for storage of multiple 
#                                   labelled histograms in memory
#
###########################################################################
from ROOT import TH1D
from ROOT import TFile
from ROOT import TTree
from ROOT import TChain
from ROOT import TRandom3
from ROOT import gDirectory

import rat

from write_spectrum import WriteSpectrum
import constants
import defaults
import isotope
import generators
import spectrum_utils
import list_manips

import math
import re
import os
import traceback

class Production(WriteSpectrum):
    """ Derived class, special case of Reconstructed, for first analysis of
    SNO+ production ntuples. Generates a histogram with Gaussian smearing,
    which is then saved to a new Root file.
    """
    def __init__(self, path):
        """ Overloaded version of init, specific to SNO+ production data

        :param path: path to SNO+ production ntuple
        :type path: str
        """
        super(Production, self).__init__(path)
        if (self._label != None): # has been set
            if (self._label.find("-Truth") >= 0):
                self._label = self._label[:self._label.find("-Truth")]
        self._majorat_name = None
    def set_parameters(self):
        """ Method to initialise spectrum parameters for non MajoRat style
        filenames.

        """
        # Assume rat release is 4.5 for now
        self._rat_release = "RAT4.5"
        self._majorat_name = self._rat_release.replace(".", "-")
        # For production data we don't really care about number of events
        self._n_events = "prod"
        self._majorat_name += "_" + self._n_events
        assert (self._options != None), "Production.set_parameters: ERROR "\
            "self._options = None\n --> Parameters already initialised in "\
            "SpectrumData.set_parameters_from_filename()"
        option = iter(self._options)
        first = option.next()
        first = first[first.find("TeLoaded")+len("TeLoaded"):]
        if (re.search("[a-zA-Z]+[0-9]+", first)): # isotope
            isotope_name = first
        else:
            isotope_name = None
        macro_name = ""
        next_option = first
        while (re.search(r"^r[0-9]+$", next_option) == None):
            if (next_option == "Solar"):
                macro_name += "_" + "solar"
            else:
                macro_name += "_" + next_option
            next_option = option.next()
        macro_name = macro_name[1:]
        rat_root = os.environ.get("RATROOT")
        macro_dir = rat_root + "/mac/production/teloaded"
        macro_path = macro_dir + "/" + macro_name + ".mac"
        for line in open(macro_path):
            if (line.find("/generator/add") >= 0):
                generator_options = line.split(" ")
        generator_add = iter(generator_options)
        next_option = generator_add.next()
        next_option = generator_add.next()
        if (next_option != "combo"):
            generator = next_option
        else:
            next_option = generator_add.next()
            assert (next_option.find(":") >= 0), \
                "Production.set_parameters: ERROR unable to determine generator"
            generator = next_option.split(":")[0]
        if (isotope_name != None):
            if (generator == "solar"):
                # Assume default type nue
                self._generator = generators.Solar(isotope_name)
                self._label = self._generator.get_isotope()
                self._majorat_name += "_" + self._generator.get_generator()
                self._majorat_name += "_" + self._generator.get_isotope()
            elif (generator == "decaychain"):
                self._generator = generators.DecayChain(isotope_name)
                self._label = self._generator.get_isotope()
                self._majorat_name += "_" + self._generator.get_generator()
                self._majorat_name += "_" + self._generator.get_isotope()
            elif (generator == "decay0"):
                for line in open(macro_path):
                    if (line.find("/generator/vtx/set") >= 0):
                        generator_vtx_set = line.split(" ")
                type_ = generator_vtx_set[1]
                #assert (generator_vtx_set[2][:-1] == isotope_name),\
                #    "Production.set_parameters: ERROR isotope name mismatch "\
                #    "generator_vtx_set[2] = " + str(generator_vtx_set[2]) + \
                #    " isotope = " + str(isotope_name)
                if (type_ == "backg"):
                    self._generator = generators.Backg(isotope_name)
                    self._label = self._generator.get_isotope()
                    self._majorat_name += "_" + self._generator.get_generator()
                    self._majorat_name += "_" + self._generator.get_type()
                    self._majorat_name += "_" + self._generator.get_isotope()
                elif (type_ == "2beta"):
                    if (isotope_name == "Te130"): # Assume we want SNO+ specific instance
                        isotope_ = isotope.SNOPlusTe(isotope_name)
                    else: # Create standard Isotope instance
                        isotope_ = isotope.Isotope(isotope_name)
                    if (len(generator_vtx_set) < 5): # All defaults
                        mode = 1
                        self._generator = generators.TwoBeta(isotope_, mode)
                        self._majorat_name += "_" + \
                            self._generator.get_generator()
                        self._majorat_name += "_" + self._generator.get_type()
                        self._majorat_name += "_" + \
                            self._generator.get_isotope().get_name()
                        self._majorat_name += "_" + \
                            str(self._generator.get_mode())
                    elif (len(generator_vtx_set) < 7): # Default energies
                        level = int(generator_vtx_set[3])
                        mode = int(generator_vtx_set[4])
                        self._generator = generators.TwoBeta(isotope_,
                                                             mode,
                                                             level)
                        self._majorat_name += "_" + \
                            self._generator.get_generator()
                        self._majorat_name += "_" + self._generator.get_type()
                        self._majorat_name += "_" + \
                            self._generator.get_isotope().get_name()
                        self._majorat_name += "_" + \
                            str(self._generator.get_level())
                        self._majorat_name += "_" + \
                            str(self._generator.get_mode())
                        self._majorat_name += "_" + \
                            str(self._generator.get_e_lo()).replace(".", "-")
                        self._majorat_name += "_" + \
                            str(self._generator.get_e_hi()).replace(".", "-")
                    elif (len(generator_vtx_set) >= 7):
                        level = int(generator_vtx_set[3])
                        mode = int(generator_vtx_set[4])
                        type_ = "2beta"
                        e_lo = float(generator_vtx_set[5])
                        e_hi = float(generator_vtx_set[6])
                        self._generator = generators.TwoBeta(isotope_, mode,
                                                             level, type_,
                                                             e_lo, e_hi)
                        self._majorat_name += "_" + \
                            self._generator.get_generator()
                        self._majorat_name += "_" + self._generator.get_type()
                        self._majorat_name += "_" + \
                            self._generator.get_isotope().get_name()
                        self._majorat_name += "_" + \
                            str(self._generator.get_level())
                        self._majorat_name += "_" + \
                            str(self._generator.get_mode())
                        self._majorat_name += "_" + \
                            str(self._generator.get_e_lo()).replace(".", "-")
                        self._majorat_name += "_" + \
                            str(self._generator.get_e_hi()).replace(".", "-")
                    self._spectral_index = spectrum_utils.get_spectral_index\
                        (self._generator.get_mode())
                    self._label = spectrum_utils.get_label(self._spectral_index)
        else:
            # add special cases
            pass
    def get_default_hist_path(self):
        """
        :returns: file "hist_ ..." filename, where (if it exists)
                  histograms would be saved.
        :rtype: str
        """
        hist_path = os.environ.get("MAJORAT_DATA") + "/hist_"\
            + self._majorat_name + self._ext
        return str(hist_path)
    def fill_histogram(self, hist_label="default"):
        """ Overloads SpectrumData.fill_histogram. Method to read ntuple,
        DS, extract total KE for each event and fill histogram.
        
        :param hist_label: histogram key in self._histograms
        :type hist_label: str
        """
        if (hist_label == "default"):
            hist_label = self._label
        apply_fv_cut = True # by default
        reconstruct_position = False # by default
        use_nhit_energy = False # by default
        include_zero_energy_events = False # by default
        ntuple_options = hist_label.split("-")
        if (list_manips.item_containing("no_fv_cut", ntuple_options) != None):
            apply_fv_cut = False
        if (list_manips.item_containing("reco_pos", ntuple_options) != None):
            reconstruct_position = True
        if (list_manips.item_containing("nhit_energy", ntuple_options) != None):
            use_nhit_energy = True
        if (list_manips.item_containing("zero_energy", ntuple_options) != None):
            include_zero_energy_events = True
        histogram = self._histograms.get(hist_label)
        assert isinstance(histogram, TH1D), \
            "SpectrumData.fill_histogram: histogram is not a TH1D object"
        chain = TChain("output") # ntuple is always called output
        chain.Add(self._path)
        for event in chain:
            if apply_fv_cut:
                if reconstruct_position:
                    pos_x = event.posx
                    pos_y = event.posy
                    pos_z = event.posz
                else:
                    pos_x = event.mcPosx
                    pos_y = event.mcPosy
                    pos_z = event.mcPosz
                radius = math.fabs(math.sqrt(pos_x**2+pos_y**2+pos_z**2))
                if use_nhit_energy:
                    energy = event.nhits / defaults.snoplus.get\
                        ("energy_resolution")
                else:
                    energy = event.energy
                fv_radius = defaults.snoplus.get("fv_radius")
                fv_radius *= 1000 # convert to mm
                if include_zero_energy_events:
                    if (radius < fv_radius):
                        histogram.Fill(energy)
                else:
                    if (radius < fv_radius) and (energy > 0.0):
                        histogram.Fill(energy)
            else:
                histogram.Fill(event.energy)

if __name__ == "__main__":
    from ROOT import TCanvas

    import argparse

    parser = argparse.ArgumentParser\
        (description=("Generate energy spectrum for the single SNO+ production "
                      "ntuple Root file specified"))
    parser.add_argument("ntuple_file1", help="path to first SNO+ production "
                        "ntuple Root file")
    parser.add_argument("ntuple_file2", help="path to second SNO+ production "
                        "ntuple Root file")
    args = parser.parse_args()
    print args
    
    spectrum1 = Production(args.ntuple_file1)
    spectrum1.set_parameters()
    hist_label = "default"
    append = False # Make new for first one
    always_remake = True # Make from scratch for first one
    spectrum1.make_histogram(hist_label, append, always_remake)
    spectrum1.write_histograms()
    spectrum2 = Production(args.ntuple_file2)
    spectrum2.set_parameters()
    hist_label = "default"
    append = True # Now we want to append to existing
    spectrum2.make_histogram(hist_label, append)
    histogram = spectrum2.get_histogram(hist_label)
    histogram.Draw()
    print "production.py: histogram contains", histogram.GetEntries(), "entries"
    spectrum2.write_histograms()
    raw_input("RETURN to exit")
