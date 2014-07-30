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
# 06/06/2014 <ab571@sussex.ac.uk> : Refactored for storage of multiple 
#                                   labelled histograms in memory
###########################################################################
from ROOT import TH1D
from ROOT import TFile
from ROOT import TRandom3

import rat

from write_spectrum import WriteSpectrum
import constants
import defaults

import math

class Reconstructed(WriteSpectrum):
    """ Derived class, special case of WriteSpectrum, for first analysis of 
    RAT generated Root file. Generates a histogram with Gaussian smearing, 
    which is then saved to a new Root file.
    """
    def __init__(self, path, t_half=None):
        """ Initialises the class, extracts information from filename """
        super(Reconstructed, self).__init__(path, t_half)
        if (self._label != None): # has been set
            if (self._label.find("-Truth") >= 0):
                self._label = self._label[:self._label.find("-Truth")]
    def fill_histogram(self, hist_label="default"):
        """ Overloads SpectrumData.fill_histogram, to read DS, extract
        total KE for each event, then apply a gaussian smearing before
        filling the histogram.
        
        """
        if (hist_label == "default"):
            hist_label = self._label
        histogram = self._histograms.get(hist_label)
        assert isinstance(histogram, TH1D), \
            "SpectrumData.fill_histogram: histogram is not a TH1D object"
        energy_resolution = defaults.analysis.get("production_label").\
            get(self._rat_release).get("energy_resolution")
        seed = 12345
        random = TRandom3(seed)
        for ds, run in rat.dsreader(self._path):
            for event in range(0, ds.GetEVCount()):
                mc = ds.GetMC();
                truth_energy = 0
                for particle in range (0, mc.GetMCParticleCount()):
                    # Check they are electrons
                    if (mc.GetMCParticle(particle).GetPDGCode() == 11):
                        truth_energy += mc.GetMCParticle(particle).GetKE()
                # Gausian smearing of events
                mu = truth_energy
                sigma = math.sqrt(truth_energy*energy_resolution)
                sigma /= energy_resolution
                reconstructed_energy = random.Gaus(mu, sigma)
                histogram.Fill(reconstructed_energy)

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
    parser.add_argument("-b", "--bin_width", help="specify histogram bin width",
                        type=float, default=0.02)
    parser.add_argument("-e", "--e_lo", help="change default lower energy limit",
                        type=float, default=None)
    parser.add_argument("-E", "--e_hi", help="change default higher energy limit",
                        type=float, default=None)
    parser.add_argument("-s", "--save_image", help="save truth/reconstructed "
                        "comparison as png", action="store_true")
    args = parser.parse_args()
    print args

    t_half = constants.zero_nu_lifetimes.get("Xe136").get(0)

    spectrum = WriteSpectrum(args.root_file, t_half)
    if (args.e_lo != None):
        spectrum.set_e_lo(args.e_lo)
    if (args.e_hi != None):
        spectrum.set_e_hi(args.e_hi)
    hist_label = "default"
    always_remake = True
    if (args.bin_width != 0.02):
        truth_hist = spectrum.get_histogram(hist_label,
                                            always_remake,
                                            args.bin_width)
    else:
        truth_hist = spectrum.get_histogram(hist_label, always_remake)
    spectrum.write_histograms()

    reco_spectrum = Reconstructed(args.root_file, t_half)
    if (args.e_lo != None):
        reco_spectrum.set_e_lo(args.e_lo)
    if (args.e_hi != None):
        reco_spectrum.set_e_hi(args.e_hi)
    hist_label = "default"
    always_remake = True
    if (args.bin_width != 0.02):
        reco_hist = reco_spectrum.get_histogram(hist_label,
                                                always_remake,
                                                args.bin_width)
    else:
        reco_hist = reco_spectrum.get_histogram(hist_label, 
                                                always_remake)
    reco_spectrum.write_histograms()

    title = "Truth/reconstructed comparison " + reco_spectrum._label
    stack = THStack(title, title)
    stack.Add(truth_hist)
    truth_hist.SetLineColor(2)
    stack.Add(reco_hist)
    reco_hist.SetLineColor(4)

    legend = TLegend(0.7, 0.85, 0.98, 0.98)
    legend.AddEntry(truth_hist, truth_hist.GetTitle(), "l")
    legend.AddEntry(reco_hist, reco_hist.GetTitle(), "l")

    c1 = TCanvas("canvas", "canvas")
    c1.cd()
    stack.Draw("nostack")
    stack.GetXaxis().SetTitle("E (MeV)")
    stack.GetYaxis().SetTitle("Events")
    legend.Draw()
    c1.Draw()
    if args.save_image:
        c1.Print("truth_reco_comparison_n"+
                 str(reco_spectrum._spectral_index)+".png")
    raw_input("RETURN to exit")

