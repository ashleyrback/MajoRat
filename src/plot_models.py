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
from ROOT import THStack
from ROOT import TLegend
from ROOT import TCanvas

from reconstructed import Reconstructed
import constants

import os

if __name__ == "__main__":
    spectra = [("hist_RAT4-5_1k_decay0_2beta_Te130_0_4_0_3-5.root", (1,1)), 
               ("hist_RAT4-5_1k_decay0_2beta_Te130_0_1_0_3-5.root", (1,2)),
               ("hist_RAT4-5_1k_decay0_2beta_Te130_0_5_0_3-5.root", (2,1)),
               ("hist_RAT4-5_1k_decay0_2beta_Te130_0_6_0_3-5.root", (4,1)),
               ("hist_RAT4-5_1k_decay0_2beta_Te130_0_7_0_3-5.root", (6,1)),
               ("hist_RAT4-5_1k_decay0_2beta_Te130_0_8_0_3-5.root", (8,1))]

    stack = THStack("Energy Spectra", "Energy Spectra")
    legend = TLegend(0.6, 0.75, 0.98, 0.98)

    t_half = constants.half_life.get("Xe136").get(0)

    for spectrum_file, style in spectra:
        spectrum = Reconstructed\
            (os.environ.get("MAJORAT_DATA")+"/"+spectrum_file, t_half)
        print spectrum
        hist = spectrum.get_histogram()
        spectrum.scale_by_t_half\
            (constants.half_life.get("Xe136").get(spectrum._spectral_index))
        print hist
        color, line = style
        hist.SetLineColor(color)
        hist.SetLineStyle(line)
        stack.Add(hist)
        legend.AddEntry(hist, hist.GetTitle(), "l")

    c1 = TCanvas("canvas", "canvas")
    c1.cd()
    stack.Draw("nostack")
    stack.GetXaxis().SetTitle("E (MeV)")
    stack.GetYaxis().SetTitle("Events/(0.1 MeV)")
    legend.Draw()
    c1.SetLogy()
    c1.Draw()
    c1.Print("EnergySpectrum.png")    

    raw_input("RETURN to exit")
