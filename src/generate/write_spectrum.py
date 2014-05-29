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

    def write_histogram(self, hist_path="default"):
        """ Writes the histogram that has been created to a separate Root 
        file prefixed with "hist_"
        """
        assert (self._histogram != None), ("GetSpectrum.write_histogram: "
                                           "histogram needs to be defined "
                                           "before it can be \nwritten to file")
        if (hist_path == "default"):
            hist_path = self.get_default_hist_path()
        else:
            directory, file_ = file_manips.split_path(hist_path)
            if (len(file_) <= 0): # Alternative directory only
                hist_path = directory + \
                    file_manips.strip_path(self.get_default_hist_path())
        output_file = TFile(hist_path, "UPDATE")
        output_file.cd()
        self._histogram.Write("", TObject.kOverwrite)
        output_file.Write()
        output_file.ls()
        output_file.Close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=("Generate energy spectrum for"
                                                  " the single (RAT-generated) "
                                                  "Root file specified"))
    parser.add_argument("root_file", help="path to RAT-generated Root file")
    args = parser.parse_args()
    print args
    print args.root_file

    spectrum = WriteSpectrum(args.root_file)
    histogram = spectrum.get_histogram()
    histogram.Draw()
    raw_input("RETURN to exit")
    spectrum.write_histogram()
