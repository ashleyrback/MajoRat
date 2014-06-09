#!/usr/bin/env python
#
# production.py
#
# Produce histogram of spectrum from SNO+ production ntuples
#
# Author A R Back 
#
# 15/05/2014 <ab571@sussex.ac.uk> : First revision
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
import isotope
import generators
import spectrum_utils

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
    def fill_histogram(self):
        """ Overloads SpectrumData.fill_histogram. Method to read ntuple,
        DS, extract total KE for each event and fill histogram.

        """
        chain = TChain("output") # ntuple is always called output
        chain.Add(self._path)
        for event in chain:
            if (event.energy != 0.0):
                self._histogram.Fill(event.energy)

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
    append = False # Make new for first one
    always_remake = True # Make from scratch for first one
    spectrum1.make_histogram(append, always_remake)
    spectrum1.write_histogram()
    spectrum2 = Production(args.ntuple_file2)
    spectrum2.set_parameters()
    append = True # Now we want to append to existing
    spectrum2.make_histogram(append)
    histogram = spectrum2.get_histogram()
    histogram.Draw()
    print "production.py: histogram contains", histogram.GetEntries(), "entries"
    spectrum2.write_histogram()
    raw_input("RETURN to exit")
