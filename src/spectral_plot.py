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

from production import Production
import constants
import defaults
import roi_utils

import os

if __name__ == "__main__":
    stack = THStack("SNO+ spectral plot", "SNO+ spectral plot")
    legend = TLegend(0.6, 0.75, 0.98, 0.98)

    roi = roi_utils.RoIUtils("reconstructed_energy")
    livetime = defaults.snoplus.get("livetime")
    print "spectral_plot.__main__: livetime =", livetime

    zero_nu = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_2beta_Te130_0_1_0-0_3-5.ntuple.root")
    effective_mass = 0.27 # eV
    zero_nu.scale_by_mass(effective_mass, livetime)
    zero_nu_hist = zero_nu.get_histogram()
    zero_nu_hist.Draw()
    zero_nu_hist.SetLineColor(1)
    stack.Add(zero_nu_hist)
    legend.AddEntry(zero_nu_hist, zero_nu_hist.GetTitle(), "l")
    raw_input("RETURN to exit")

    two_nu = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_2beta_Te130_0_4_0-0_3-5.ntuple.root")
    total_events = 3.8e6
    two_nu.scale_by_events(total_events*livetime)
    two_nu_hist = two_nu.get_histogram()
    two_nu_hist.Draw()
    two_nu_hist.SetLineColor(2)
    stack.Add(two_nu_hist)
    legend.AddEntry(two_nu_hist, two_nu_hist.GetTitle(), "l")
    raw_input("RETURN to exit")

    solar = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_solar_B8.ntuple.root")
    events_in_roi = 4.04
    solar.scale_by_roi_events(events_in_roi*livetime, roi)
    solar_hist = solar.get_histogram()
    solar_hist.Draw()
    solar_hist.SetLineColor(3)
    stack.Add(solar_hist)
    legend.AddEntry(solar_hist, solar_hist.GetTitle(), "l")
    raw_input("RETURN to exit")

    c1 = TCanvas("canvas", "canvas")
    c1.cd()
    stack.Draw("nostack")
    stack.GetXaxis().SetTitle("E (MeV)")
    stack.GetYaxis().SetTitle("Events/(0.02 MeV)")
    legend.Draw()
    c1.SetLogy()
    c1.Draw()
    c1.Print("spectral_plot.png")    

    raw_input("RETURN to exit")
