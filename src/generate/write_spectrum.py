#!/usr/bin/env python
#
# write_spectrum.py
#
# Produce histogram of spectrum from Rat generated root file
#
# Author A R Back 
#
# 31/01/2014 <ab571@sussex.ac.uk> : First revision
# 29/04/2014 <ab571@sussex.ac.uk> : Now initialised with half life, as in
#                                   base class
# 06/06/2014 <ab571@sussex.ac.uk> : Refactored for storage of multiple 
#                                   labelled histograms in memory
###########################################################################
from ROOT import TFile
from ROOT import TObject

import rat

from spectrum_data import SpectrumData
import file_manips

class WriteSpectrum(SpectrumData):
    """ Derived class, special case of SpectrumData, for first analysis of 
    RAT generated Root file. Alows for easy generation of histograms, which
    are then saved to a new Root file.
    """
    def __init__(self, path, t_half=None):
        """ Initialises the class, extracts information from filename """
        super(WriteSpectrum, self).__init__(path, t_half)
    def write_histogram(self, hist_label="default",
                        hist_path="default"):
        """ Writes the histogram that has been created to a separate ROOT 
        file prefixed with "hist_"

        :param hist_label: self._histograms key (and ROOT label)
        :type hist_label: str
        :param hist_path: directory/path for ROOT file (default is 
                          $MAJORAT_DATA)
        :type hist_path: str
        """
        if (hist_label == "default"):
            hist_label = self._label
        histogram = self._histograms.get(hist_label)
        assert (histogram != None),\
            "WriteSpectrum.write_histogram: histogram " + hist_label +" not found"
        if (hist_path == "default"):
            hist_path = self.get_default_hist_path()
        else:
            directory, file_ = file_manips.split_path(hist_path)
            if (len(file_) <= 0): # Alternative directory only
                hist_path = directory + \
                    file_manips.strip_path(self.get_default_hist_path())
        output_file = TFile(hist_path, "UPDATE")
        print "WriteSpectrum.write_histogram: writing to ..."
        print " -->", hist_path
        output_file.cd()
        histogram.Write("", TObject.kOverwrite)
        output_file.Write()
        output_file.Close()
    def write_histograms(self, hist_path="default"):
        """ Writes all histograms in memory to a separate ROOT file 
        prefixed with "hist_".
        
        :param hist_path: directory/path for ROOT file (default is 
                          $MAJORAT_DATA)
        :type hist_path: str
        """
        assert (len(self._histograms.items()) > 0), \
            "WriteSpectrum.write_histograms: no histograms created yet"
        for label, histogram in self._histograms.iteritems():
            print "WriteSpectrum.write_histograms: writing ... ", label
            self.write_histogram(label, hist_path)
    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=("Generate energy spectrum for"
                                                  " the single (RAT-generated) "
                                                  "Root file specified"))
    parser.add_argument("root_file", help="path to RAT-generated Root file")
    args = parser.parse_args()
    print args

    t_half = 5.0e25 # years # KLZ limit from Gando et al.

    spectrum = WriteSpectrum(args.root_file)
    spectrum_isotope = spectrum.get_isotope()
    spectrum_isotope.set_scintillator_masses()
    spectrum_isotope.set_number_nuclei(spectrum_isotope.get_mass())
    zero_nu = spectrum_isotope.get_zero_nu()
    spectrum_isotope.set_effective_mass(zero_nu.half_life_to_mass(t_half))
    spectrum_isotope.set_counts_by_mass()
    spectrum.scale_histogram(spectrum_isotope.get_counts())
    histogram = spectrum.get_histogram()
    histogram.Draw()
    raw_input("RETURN to exit")
    spectrum.write_histograms()
