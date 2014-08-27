#!/usr/bin/env python
#
# KLZComp.py
#
# MajoRat tool for comparison with KamLAND-Zen (KLZ). Prepares sig & 
# combined KLZ background histograms, calculates 90% CLs on signal. 
# Produces spectral plot.
#
# Author A R Back 
#
# 28/07/2014 <ab571@sussex.ac.uk> : First revision
###########################################################################
import ROOT
from ROOT import TFile
from ROOT import THStack
from ROOT import TLegend

from reconstructed import Reconstructed
from set_limit import SetLimit
import roi_utils
import defaults
import plots
import isotope

import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser\
        (description=("KLZComp: tool for calculating sensitivity to "
                      "different neutrinoless double beta decay signals "
                      "compared to the combined KamLAND-Zen background."))
    parser.add_argument("-m", "--mode", help="spectral_plot mode or "
                        "limit_setting mode. Plots are also produced in "
                        "limit_setting mode.")
    args = parser.parse_args()
                        
    # Create stack
    stack = THStack("klz_spectral_plot", "KamLAND-Zen spectral plot (reproduced)")
    # Create legend
    legend = TLegend(0.0, 0.0, 1.0, 1.0)
    legend.SetFillColor(0)
    
    livetime = isotope.convert_to_years("112.3 d")
    print "KLZComp.__main__: livetime =", livetime

    Xe136_mass = 125.0 # +/- 7kg 
    
    # Backgrounds ##############################################################
    background_file = TFile\
        (os.environ.get("MAJORAT_DATA")+"/KLZ_combined_backgrounds.root")
    total_background_hist = background_file.Get("Total B.G.")
    assert isinstance(total_background_hist, ROOT.TH1D), \
        "KLZComp.__main__: error - histogram with key 'Total B.G.' not found"
    total_background_hist.SetDirectory(0)
    total_background_hist.SetTitle("Total B.G.")
    total_background_hist.SetLineColor(ROOT.kBlack)
    total_background_hist.SetLineStyle(2)
    total_background_hist.SetLineWidth(2)
    total_background_hist.SetStats(0)
    stack.Add(total_background_hist)
    legend.AddEntry(total_background_hist, 
                    total_background_hist.GetTitle(), "l")

    two_nu = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_100k_decay0_2beta_Xe136_0_4_0-5_4-75.root")
    Xe136_n5 = two_nu.get_isotope()
    Xe136_n5.set_number_nuclei(Xe136_mass)
    Xe136_n5.set_half_life("2.30e21 y")
    Xe136_n5.set_counts(livetime)
    two_nu.scale_histogram(Xe136_n5.get_counts())
    two_nu_hist = two_nu.get_histogram()
    two_nu_hist.SetDirectory(0)
    two_nu_hist.SetLineWidth(2)
    two_nu_hist.SetLineColor(ROOT.kMagenta)
    two_nu_hist.SetStats(0)
    stack.Add(two_nu_hist)
    legend.AddEntry(two_nu_hist, two_nu_hist.GetTitle(), "l")

    assert(two_nu_hist.GetBinWidth(1) == total_background_hist.GetBinWidth(1)),\
        "different bin widths"

    sum_background_hist = total_background_hist \
        + two_nu_hist
    sum_background_hist.SetDirectory(0)
    sum_background_hist.SetTitle("Sum, background (2#nu#beta#beta + B.G.)")
    sum_background_hist.SetLineColor(ROOT.kRed)
    sum_background_hist.SetLineStyle(1)
    sum_background_hist.SetLineWidth(2)
    stack.Add(sum_background_hist)
    legend.AddEntry(sum_background_hist, sum_background_hist.GetTitle(), "l")
    ############################################################################

    # Signal ###################################################################
    zero_nu = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Xe136_0_1_0-5_4-75.root")
    Xe136_0nu = zero_nu.get_isotope()
    Xe136_0nu.set_number_nuclei(Xe136_mass)
    n1 = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Xe136_0_5_0-5_4-75.root")
    Xe136_n1 = n1.get_isotope()
    Xe136_n1.set_number_nuclei(Xe136_mass)
    n2 = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Xe136_0_6_0-5_4-75.root")
    Xe136_n2 = n2.get_isotope()
    Xe136_n2.set_number_nuclei(Xe136_mass)
    n3 = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Xe136_0_7_0-5_4-75.root")
    Xe136_n3 = n3.get_isotope()
    Xe136_n3.set_number_nuclei(Xe136_mass)
    n7 = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Xe136_0_8_0-5_4-75.root")
    Xe136_n7 = n7.get_isotope()
    Xe136_n7.set_number_nuclei(Xe136_mass)
    if (args.mode == "limit_setting"):
        zero_nu_limit_setter = SetLimit(sum_background_hist, zero_nu, livetime)
        zero_nu_limit_setter.set_t_half_limit(2.0e24, 1.2e25, 0.1e24, 2.71)
        zero_nu_hist = zero_nu_limit_setter.get_signal_hist()
        n1_limit_setter = SetLimit(sum_background_hist, n1, livetime)
        n1_limit_setter.set_t_half_limit(0.1e24, 1.0e25, 0.1e24, 2.71)
        n1_hist = n1_limit_setter.get_signal_hist()
        n2_limit_setter = SetLimit(sum_background_hist, n2, livetime)
        n2_limit_setter.set_t_half_limit(0.5e24, 1.5e24, 0.1e23, 2.71)
        n2_hist = n2_limit_setter.get_signal_hist()
        n3_limit_setter = SetLimit(sum_background_hist, n3, livetime)
        n3_limit_setter.set_t_half_limit(1.0e23, 1.0e24, 0.1e23, 2.71)
        n3_hist = n3_limit_setter.get_signal_hist()
        n7_limit_setter = SetLimit(sum_background_hist, n7, livetime)
        n7_limit_setter.set_t_half_limit(1.0e22, 2.0e23, 0.1e21, 2.71)
        n7_hist = n7_limit_setter.get_signal_hist()
        sum_hist = sum_background_hist + \
            zero_nu_hist +\
            n1_hist +\
            n2_hist+\
            n3_hist+\
            n7_hist
        sum_hist.SetTitle("Sum, 90 % CLs")
    else:
        Xe136_0nu.set_half_life("6.2e24 y")
        Xe136_0nu.set_counts(livetime)
        zero_nu.scale_histogram(Xe136_0nu.get_counts())
        zero_nu_hist = zero_nu.get_histogram()
        zero_nu_hist.SetDirectory(0)
        Xe136_n1.set_half_life("2.6e24 y")
        Xe136_n1.set_counts(livetime)
        n1.scale_histogram(Xe136_n1.get_counts())
        n1_hist = n1.get_histogram()
        n1_hist.SetDirectory(0)
        Xe136_n2.set_half_life("1.0e24 y")
        Xe136_n2.set_counts(livetime)
        n2.scale_histogram(Xe136_n2.get_counts())
        n2_hist = n2.get_histogram()
        n2_hist.SetDirectory(0)
        Xe136_n3.set_half_life("4.5e23 y")
        Xe136_n3.set_counts(livetime)
        n3.scale_histogram(Xe136_n3.get_counts())
        n3_hist = n3.get_histogram()
        n3_hist.SetDirectory(0)
        Xe136_n7.set_half_life("1.1e22 y")
        Xe136_n7.set_counts(livetime)
        n7.scale_histogram(Xe136_n7.get_counts())
        n7_hist = n7.get_histogram()
        n7_hist.SetDirectory(0)
        sum_hist = sum_background_hist + \
            zero_nu_hist +\
            n1_hist +\
            n2_hist+\
            n3_hist+\
            n7_hist
        sum_hist.SetTitle("Sum")
    zero_nu_hist.SetLineWidth(2)
    zero_nu_hist.SetLineColor(ROOT.kAzure+9)
    zero_nu_hist.SetLineStyle(3)
    n1_hist.SetLineWidth(2)
    n1_hist.SetLineColor(ROOT.kAzure-0)
    n1_hist.SetLineStyle(9)
    n2_hist.SetLineWidth(2)
    n2_hist.SetLineColor(ROOT.kAzure-2)
    n2_hist.SetLineStyle(5)
    n3_hist.SetLineWidth(2)
    n3_hist.SetLineColor(ROOT.kBlue-3)
    n3_hist.SetLineStyle(10)
    n7_hist.SetLineWidth(2)
    n7_hist.SetLineColor(ROOT.kBlue-4)
    n7_hist.SetLineStyle(7)
    sum_hist.SetLineStyle(1)
    sum_hist.SetLineWidth(1)
    sum_hist.SetLineColor(ROOT.kBlack)
    stack.Add(zero_nu_hist)
    stack.Add(n1_hist)
    stack.Add(n2_hist)
    stack.Add(n3_hist)
    stack.Add(n7_hist)
    legend.AddEntry(zero_nu_hist, zero_nu_hist.GetTitle(), "l")
    legend.AddEntry(n1_hist, n1_hist.GetTitle(), "l")
    legend.AddEntry(n2_hist, n2_hist.GetTitle(), "l")
    legend.AddEntry(n3_hist, n3_hist.GetTitle(), "l")
    legend.AddEntry(n7_hist, n7_hist.GetTitle(), "l")
    stack.Add(sum_hist)
    legend.AddEntry(sum_hist, sum_hist.GetTitle(), "l")
    ############################################################################

    # Spectral plot ############################################################
    if (args.mode == "limit_setting"):
        plots.make_KLZ_spectral_plot("KLZ_90CL_spectral_plot", stack, legend)
    else:
        plots.make_KLZ_spectral_plot("KLZ_spectral_plot", stack, legend)
