#!/usr/bin/env python
#
# spectrum_plot.py
#
# Script to generate the SNO+ spectral plot
#
# Author A R Back 
#
# 05/06/2014 <ab571@sussex.ac.uk> : First revision
###############################################################################
from ROOT import TH1D
from ROOT import THStack
from ROOT import TLegend
from ROOT import TCanvas
from ROOT import TPad
from ROOT import gStyle

from production import Production
import set_limit
import constants
import defaults
import roi_utils

import os

def my_decreasing_range(start, end, step):
    while start >= end:
        yield start
        start -= step

if __name__ == "__main__":
    stack = THStack("SNO+ spectral plot", "SNO+ spectral plot")
    legend = TLegend(0.0, 0.0, 1.0, 1.0)
    legend.SetFillColor(0)
    gStyle.SetOptStat("")

    roi = roi_utils.RoIUtils("reconstructed_energy")
    livetime = defaults.ll_analysis.get("livetime")
    print "spectral_plot.__main__: livetime =", livetime

    two_nu = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_2beta_Te130_0_4_0-0_3-5.ntuple.root")
    total_events = 3.8e6
    hist_label = two_nu._label + "-reco_pos"
    two_nu.scale_by_events(total_events*livetime, hist_label)
    two_nu_hist = two_nu.get_histogram(hist_label)
    two_nu_hist.Draw()
    two_nu_hist.SetLineWidth(2)
    two_nu_hist.SetLineColor(2)
    stack.Add(two_nu_hist)
    legend.AddEntry(two_nu_hist, two_nu_hist.GetTitle(), "l")

    solar = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_solar_B8.ntuple.root")
    hist_label = solar._label + "-reco_pos"
    events_in_roi = 4.04
    solar.scale_by_roi_events(events_in_roi*livetime, roi)
    solar_hist = solar.get_histogram()
    solar_hist.Draw()
    solar_hist.SetLineWidth(2)
    solar_hist.SetLineColor(3)
    stack.Add(solar_hist)
    legend.AddEntry(solar_hist, solar_hist.GetTitle(), "l")

    Bi214 = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_backg_Bi214.ntuple.root")
    hist_label = Bi214._label + "-reco-pos"
    events_in_roi = 1.33
    Bi214.scale_by_events(events_in_roi*livetime, roi)
    Bi214_hist = Bi214.get_histogram(hist_label)
    Bi214_hist.Draw()
    Bi214_hist.SetLineWidth(2)
    Bi214_hist.SetFillColor(4)
    
 
    sum_background_hist = two_nu_hist + solar_hist
    sum_background_hist.SetTitle("Sum, background")
    sum_background_hist.SetLineColor(1)
    sum_background_hist.SetLineStyle(2)
    sum_background_hist.SetLineWidth(1)
    stack.Add(sum_background_hist)
    legend.AddEntry(sum_background_hist, sum_background_hist.GetTitle(), "l")
    
    no_limit_set = True
    effective_mass_limit = 0.0
    for effective_mass in my_decreasing_range(0.065, 0.045, 0.001):
        zero_nu = Production\
            (os.environ.get("MAJORAT_DATA")+\
                 "/hist_RAT4-5_prod_decay0_2beta_Te130_0_1_0-0_3-5.ntuple.root")
        hist_label = "default"
        zero_nu.scale_by_mass(effective_mass, hist_label, livetime)
        zero_nu_hist = zero_nu.get_histogram()
        zero_nu_hist.Draw()
        zero_nu_hist.SetLineWidth(2)
        zero_nu_hist.SetLineColor(1)

        sum_hist = sum_background_hist + zero_nu_hist
        sum_hist.SetTitle("Sum, " + zero_nu._label + "-" + \
                          str(int(effective_mass*1000)) + " meV")
        sum_hist.SetLineStyle(1)
        sum_hist.SetLineWidth(1)

        reduced_chi_squared = set_limit.log_likelihood_mass(sum_hist,
                                                            sum_background_hist,
                                                            zero_nu)
        print reduced_chi_squared
        if no_limit_set:
            effective_mass_limit = effective_mass
            if (reduced_chi_squared < 2.71): # 90% C.L. for delta-chi-squared
                no_limit_set = False
                stack.Add(zero_nu_hist)
                legend.AddEntry(zero_nu_hist, zero_nu_hist.GetTitle(), "l")
                stack.Add(sum_hist)
                legend.AddEntry(sum_hist, sum_hist.GetTitle(), "l")
                c1 = TCanvas("canvas", "canvas", 1200, 900)
                c1.cd()
                p1 = TPad("p1", "Region of Interest", 0.0, 0.5, 0.7, 1.0, 0)
                p1.Draw()
                p2 = TPad("p2", "Key", 0.7, 0.5, 1.0, 1.0, 0)
                p2.SetBottomMargin(0.15)
                p2.SetTopMargin(0.15)
                p2.SetLeftMargin(0.05)
                p2.SetRightMargin(0.05)
                p2.Draw()
                p3 = TPad("p3", "Full Spectrum", 0.0, 0.0, 1.0, 0.5, 0)
                p3.Draw()
                p3.cd()
                p3.SetLogy()
                stack.Draw("nostack")
                stack.GetXaxis().SetTitle("E (MeV)")
                stack.GetYaxis().SetTitle("Events/(0.02 MeV)")
                stack.SetMinimum(0.1)
                p1.cd()
                p1.SetLogy()
                stack_zoomed = stack.Clone()
                x_axis = stack_zoomed.GetXaxis()
                bin_width = x_axis.GetBinWidth(0) # assume uniform
                e_lo = 2.0 # MeV
                e_hi = 3.5 # MeV
                x_axis.SetRange(int(e_lo/bin_width), int(e_hi/bin_width))
                stack_zoomed.SetMaximum(5e3)
                stack_zoomed.Draw("nostack")
                p2.cd()
                legend.Draw()
                c1.Draw()
                c1.Print("spectral_plot.png")
        print "spectral_plot.py: effective mass =", effective_mass_limit*1e3,\
            "at 90% C.L"
    raw_input("RETURN to exit")
