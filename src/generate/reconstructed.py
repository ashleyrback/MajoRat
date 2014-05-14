#!/usr/bin/env python
#
# reconstructed.py
#
# Produce histogram of spectrum from Rat generated root file
#
# Author A R Back 
#
# 31/01/2014 <ab571@sussex.ac.uk> : First revision
# 29/04/2014 <ab571@sussex.ac.uk> : Now initialised with half life, as base
#                                   class
###########################################################################
from ROOT import TH1D
from ROOT import TFile
from ROOT import TRandom3

import rat

from write_spectrum import WriteSpectrum
import constants

import math

class Reconstructed(WriteSpectrum):
    """ Derived class, special case of WriteSpectrum, for first analysis of 
    RAT generated Root file. Generates a histogram with Gaussian smearing, 
    which is then saved to a new Root file.
    """
    def __init__(self, path, t_half):
        """ Initialises the class, extracts information from filename """
        super(Reconstructed, self).__init__(path, t_half)

    def write_histogram(self, prefix="hist_"):
        """ Writes the histogram that has been created to a separate Root file
        """
        super(Reconstructed, self).write_histogram(prefix)

    def get_histogram(self):
        """ Overloads get_histogram() method so that it generates a gaussian
        smeared histogram.
        """
        print self._prefix
        energy_resolution = constants.snoplus.get("energy_resolution")
        if (self._histogram == None):
            if (self._prefix.find("hist") >= 0):
                input_file = TFile(self._path, "READ")
                input_file.ls()
                self._histogram = input_file.Get(self._label)
                print "Reconstructed.get_histogram: found saved histogram " \
                    + self._label
                self._histogram.SetDirectory(0)
                input_file.Close()
            else:
                seed = 12345
                random = TRandom3(seed)
                self._histogram = TH1D(self._label, self._label, 30, 0.0, 3.0)
                self._histogram.SetDirectory(0)
                for ds, run in rat.dsreader(self._path):
                    for iEV in range(0, ds.GetEVCount()):
                        mc = ds.GetMC();
                        KE = 0
                        for particle in range (0, mc.GetMCParticleCount()):
                            # Check they are electrons
                            if (mc.GetMCParticle(particle).GetPDGCode() == 11):
                                KE += mc.GetMCParticle(particle).GetKE()
                        # Gausian smearing of events
                        E_true = KE
                        mu = E_true
                        sigma = math.sqrt(E_true*energy_resolution)
                        sigma /= energy_resolution
                        E_reco = random.Gaus(mu, sigma)
                        self._histogram.Fill(E_reco)
            print "Integral of histogram = ", self._histogram.Integral()
            print "Unscaled histogram = ", self._histogram_unscaled
            if (self._histogram_unscaled == None):
                self._histogram_unscaled = self._histogram.Clone()
            print "Unscaled histogram = ", self._histogram_unscaled
        else:
            print ("SpectrumData.get_histogram: SpectrumData._histogram object "
                   "already exists\n --> re-using!")
        print self._histogram
        return self._histogram


if __name__ == "__main__":
    from ROOT import TCanvas
    from ROOT import THStack
    from ROOT import TLegend

    import constants

    import argparse
    
    parser = argparse.ArgumentParser(description=("Generate energy spectrum for"
                                                  " the single (RAT-generated) "
                                                  "Root file specified"))
    parser.add_argument("root_file", help="path to RAT-generated Root file")
    args = parser.parse_args()
    print args
    print args.root_file

    t_half = constants.half_life.get("Xe136").get(0)

    spectrum = WriteSpectrum(args.root_file, t_half)
    truth_hist = spectrum.get_histogram()
    spectrum.write_histogram()

    reco_spectrum = Reconstructed(args.root_file, t_half)
    reco_hist = reco_spectrum.get_histogram()
    reco_spectrum.write_histogram()

    title = "Truth/reconstructed comparison " + reco_spectrum._label
    stack = THStack(title, title)
    stack.Add(truth_hist)
    truth_hist.SetLineColor(2)
    stack.Add(reco_hist)
    reco_hist.SetLineColor(4)

    legend = TLegend(0.7, 0.85, 0.98, 0.98)
    legend.AddEntry(truth_hist, truth_hist.GetTitle(), "l")
    legend.AddEntry(reco_hist, reco_hist.GetTitle(), "l")

    c1 = TCanvas()
    c1.cd()
    stack.Draw("nostack")
    stack.GetXaxis().SetTitle("E (MeV)")
    stack.GetYaxis().SetTitle("Events")
    legend.Draw()
    c1.Draw()
    c1.Print("truth_reco_comparison_n"+
             str(reco_spectrum._spectral_index)+".png")
    raw_input("RETURN to exit")
