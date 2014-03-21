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

from reconstructed_spectrum import Reconstructed

spectra = [("hist_RAT4-5_1k_decay0_2beta_Te130_0_4_ELT3-5.root", (1,1)), 
           ("hist_RAT4-5_1k_decay0_2beta_Te130_0_1_ELT3-5.root", (1,2)),
           ("hist_RAT4-5_1k_decay0_2beta_Te130_0_5_ELT3-5.root", (2,1)),
           ("hist_RAT4-5_1k_decay0_2beta_Te130_0_6_ELT3-5.root", (4,1)),
           ("hist_RAT4-5_1k_decay0_2beta_Te130_0_7_ELT3-5.root", (6,1)),
           ("hist_RAT4-5_1k_decay0_2beta_Te130_0_8_ELT3-5.root", (8,1))]

data_dir = "/home/ashley/snoplus/data/"

stack = THStack("Energy Spectra", "Energy Spectra")
legend = TLegend(0.6, 0.75, 0.98, 0.98)

for spectrum_file, style in spectra:
    spectrum = Reconstructed(data_dir+spectrum_file)
    print spectrum
    hist = spectrum.get_histogram()
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
c1.Draw()
c1.Print("EnergySpectrum.png")    

raw_input("RETURN to exit")
