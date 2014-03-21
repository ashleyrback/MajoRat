#!/usr/bin/env python
#
# write_spectrum.py
#
# Produce histogram of spectrum from Rat generated root file
#
# Author A R Back - 31/01/2014 <ab571@sussex.ac.uk> : First revision
###############################################################################
from ROOT import TFile

import rat

from spectrum_data import SpectrumData

class WriteSpectrum(SpectrumData):
    """ Derived class, special case of SpectrumData, for first analysis of 
    RAT generated Root file. Alows for easy generation of histograms, which
    are then saved to a new Root file.
    """
    def __init__(self, path):
        """ Initialises the class, extracts information from filename """
        super(WriteSpectrum, self).__init__(path)

    def write_histogram(self, prefix):
        """ Writes the histogram that has been created to a separate Root file
        """
        assert (self._histogram != None), ("GetSpectrum.write_histogram: "
                                           "histogram needs to be defined "
                                           "before it can be \nwritten to file")
        filename = prefix + self._filename
        path = self._dir + filename
        output_file = TFile(path, "RECREATE")
        output_file.cd()
        self._histogram.Write()
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
