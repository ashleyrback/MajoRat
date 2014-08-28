#!/usr/bin/env python
#
# MajoRat.py
#
# MajoRat tool to compare with SNO+ neutrinoless double beta decay 
# sensitivity. Prepares sig & background histograms, calculates 90% CLs on 
# signal. Produces spectral plot.
#
# Author A R Back 
#
# 08/07/2014 <ab571@sussex.ac.uk> : First revision
###########################################################################
import ROOT
from ROOT import THStack
from ROOT import TLegend

from production import Production
from set_limit import SetLimit
import roi_utils
import defaults
import plots

import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser\
        (description=("MajoRat: tool for calculating confidence limits for "
                      "different neutrinoless double beta decay signals in "
                      "SNO+ and generating spectral plots to visualise them."))
    parser.add_argument("-m", "--mode", help="spectral_plot mode or "
                        "limit_setting mode. Plots are also produced in "
                        "limit_setting mode.")
    args = parser.parse_args()
                        
    # Create stack
    stack = THStack("SNO+_spectral_plot", "SNO+ spectral plot")
    # Create legend
    legend = TLegend(0.0, 0.0, 1.0, 1.0)
    legend.SetFillColor(0)
    
    roi = roi_utils.RoIUtils("reconstructed_energy")
    if (args.mode == "limit_setting"):
        livetime = defaults.ll_analysis.get("livetime")
    else:
        livetime = defaults.spectral_plot.get("livetime")
    print "SNOPlusSensitivity.__main__: livetime =", livetime
    
    # Backgrounds ##############################################################
    two_nu = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_2beta_Te130_0_4_0-0_3-5.ntuple.root")
    hist_label = two_nu._label + "-reco_pos"
    Te130_2nu = two_nu.get_isotope()
    scintillator_cocktail = "default"
    apply_fv_cut = True
    Te130_2nu.set_scintillator_masses(scintillator_cocktail, apply_fv_cut)
    Te130_2nu.set_number_nuclei(Te130_2nu.get_mass())
    Te130_2nu.set_half_life(Te130_2nu.get_two_nu().get_half_life())
    Te130_2nu.set_counts(livetime)
    two_nu.scale_histogram(Te130_2nu.get_counts(), hist_label)
    two_nu_hist = two_nu.get_histogram(hist_label)
    two_nu_hist.SetLineWidth(2)
    two_nu_hist.SetLineColor(2)
    stack.Add(two_nu_hist)
    legend.AddEntry(two_nu_hist, two_nu_hist.GetTitle(), "l")

    solar = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_solar_B8.ntuple.root")
    hist_label = solar._label + "-reco_pos"
    total_events = 4.04
    solar.scale_histogram(total_events*livetime)
    solar_hist = solar.get_histogram(hist_label)
    solar_hist.SetLineWidth(2)
    solar_hist.SetLineColor(3)
    stack.Add(solar_hist)
    legend.AddEntry(solar_hist, solar_hist.GetTitle(), "l")

    Bi214 = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_backg_Bi214.ntuple.root")
    hist_label = Bi214._label + "-reco_pos"
    events_in_roi = 1.33
    Bi214.scale_by_roi_events(events_in_roi*livetime, roi)
    Bi214_hist = Bi214.get_histogram()
    Bi214_hist.SetLineWidth(2)
    Bi214_hist.SetLineColor(4)
    stack.Add(Bi214_hist)
    legend.AddEntry(Bi214_hist, Bi214_hist.GetTitle(), "l")

    Tl208 = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_backg_Tl208.ntuple.root")
    hist_label = Tl208._label + "-reco_pos"
    events_in_roi = 0.002
    Tl208.scale_by_roi_events(events_in_roi*livetime, roi)
    Tl208_hist = Tl208.get_histogram()
    Tl208_hist.SetLineWidth(2)
    Tl208_hist.SetLineColor(ROOT.kOrange+7)
    stack.Add(Tl208_hist)
    legend.AddEntry(Tl208_hist, Tl208_hist.GetTitle(), "l")

    Bi212po212 = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_backg_Bi212po212.ntuple.root")
    hist_label = Bi212po212._label + "-reco_pos"
    events_in_roi = 0.454
    Bi212po212.scale_by_roi_events(events_in_roi*livetime, roi)
    Bi212po212_hist = Bi212po212.get_histogram()
    Bi212po212_hist.SetLineWidth(2)
    Bi212po212_hist.SetLineColor(ROOT.kCyan-7)
    stack.Add(Bi212po212_hist)
    legend.AddEntry(Bi212po212_hist, Bi212po212_hist.GetTitle(), "l")
    
    sum_background_hist = two_nu_hist \
        + solar_hist\
        + Bi214_hist\
        + Tl208_hist\
        + Bi212po212_hist
    sum_background_hist.SetDirectory(0)
    sum_background_hist.SetTitle("Sum, background")
    sum_background_hist.SetLineColor(1)
    sum_background_hist.SetLineStyle(2)
    sum_background_hist.SetLineWidth(1)
    stack.Add(sum_background_hist)
    legend.AddEntry(sum_background_hist, sum_background_hist.GetTitle(), "l")
    ############################################################################

    # Signal ###################################################################
    zero_nu = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_2beta_Te130_0_1_0-0_3-5.ntuple.root")
    if (args.mode == "limit_setting"):
        limit_setter = SetLimit(sum_background_hist, zero_nu)
        limit_setter.set_mass_limit(2.71,0.055, 0.075, 0.001)
        zero_nu_hist = limit_setter.get_signal_hist()
        sum_hist = limit_setter.get_sum_hist()
    else:
        effective_mass = defaults.spectral_plot.get("effective_mass")
        hist_label = "default"
        zero_nu.scale_by_mass(effective_mass, hist_label, livetime)
        zero_nu_hist = zero_nu.get_histogram()
        zero_nu_hist.SetLineWidth(2)
        zero_nu_hist.SetLineColor(1)
        sum_hist = sum_background_hist + zero_nu_hist
        sum_hist.SetLineStyle(1)
        sum_hist.SetLineWidth(1)
        sum_hist.SetTitle("Sum, " + zero_nu._label + "-" + \
                          str(int(effective_mass*1000)) + " meV")
    stack.Add(zero_nu_hist)
    legend.AddEntry(zero_nu_hist, zero_nu_hist.GetTitle(), "l")
    stack.Add(sum_hist)
    legend.AddEntry(sum_hist, sum_hist.GetTitle(), "l")
    ############################################################################

    # Spectral plot ############################################################
    plots.make_spectral_plot("SNO+_spectral_plot", stack, legend)
