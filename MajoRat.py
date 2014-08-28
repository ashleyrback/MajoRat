#!/usr/bin/env python
#
# MajoRat.py
#
# Main MajoRat program - for analysing SNO+ Majoron modes. Prepares Majoron
# modes & relevant SNO+ background histograms, calculates 90% CLs on each 
# Majoron mode. Produces spectral plot.
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
from production import Production
from set_limit import SetLimit
import roi_utils
import defaults
import plots
import isotope

import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser\
        (description=("MajoRat: tool for calculating sensitivity to different "
                      "Majoron neutrinoless double beta decay signals in SNO+"))
    parser.add_argument("-m", "--mode", help="spectral_plot mode or "
                        "limit_setting mode. Plots are also produced in "
                        "limit_setting mode.")
    args = parser.parse_args()
                        
    # Create stack
    stack = THStack("SNO+_majorons_spectral_plot", "SNO+ Majorons spectral plot")
    # Create legend
    legend = TLegend(0.0, 0.0, 1.0, 1.0)
    legend.SetFillColor(0)
    
    #roi = roi_utils.RoIUtils("reconstructed_energy")
    #if (args.mode == "limit_setting"):
    #    livetime = defaults.ll_analysis.get("livetime")
    #else:
    #    livetime = defaults.spectral_plot.get("livetime")
    livetime = isotope.convert_to_years("112.3 d")
    print "MajoRat.__main__: livetime =", livetime

    # All spectra ##############################################################
    apply_fv_cut = True
    scintillator_cocktail = "default"
    ############################################################################

    # Backgrounds ##############################################################
    two_nu = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_decay0_2beta_Te130_0_4_0-0_3-5.ntuple.root")
    hist_label = two_nu._label + "-reco_pos"
    Te130_2nu = two_nu.get_isotope()
    Te130_2nu.set_scintillator_masses(scintillator_cocktail, apply_fv_cut)
    Te130_2nu.set_number_nuclei(Te130_2nu.get_mass())
    Te130_2nu.set_half_life("2.30e21 y")
    Te130_2nu.set_counts(livetime)
    two_nu.scale_histogram(Te130_2nu.get_counts(), hist_label)
    two_nu_hist = two_nu.get_histogram(hist_label)
    two_nu_hist.SetLineWidth(2)
    two_nu_hist.SetLineColor(ROOT.kRed)
    stack.Add(two_nu_hist)
    legend.AddEntry(two_nu_hist, two_nu_hist.GetTitle(), "l")

    solar = Production\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_prod_solar_B8.ntuple.root")
    hist_label = solar._label + "-reco_pos"
    total_events = 847 # 0-5 MeV range
    total_events *= (3.5/5.0) # Convert to 0-3.5 MeV range
    solar.scale_histogram(total_events*livetime)
    solar_hist = solar.get_histogram(hist_label)
    solar_hist.SetLineWidth(2)
    solar_hist.SetLineColor(ROOT.kGreen)
    stack.Add(solar_hist)
    legend.AddEntry(solar_hist, solar_hist.GetTitle(), "l")

    sum_background_hist = two_nu_hist + solar_hist
    sum_background_hist.SetDirectory(0)
    sum_background_hist.SetTitle("Sum, background (2#nu#beta#beta + B.G.)")
    sum_background_hist.SetLineColor(ROOT.kBlack)
    sum_background_hist.SetLineStyle(1)
    sum_background_hist.SetLineWidth(2)
    stack.Add(sum_background_hist)
    legend.AddEntry(sum_background_hist, sum_background_hist.GetTitle(), "l")
    ############################################################################

    # Signal ###################################################################
    n1 = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Te130_0_5_0-0_3-5.root")
    Te130_n1 = n1.get_isotope()
    Te130_n1.set_scintillator_masses(scintillator_cocktail, apply_fv_cut)
    Te130_n1.set_number_nuclei(Te130_n1.get_mass())
    n2 = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Te130_0_6_0-0_3-5.root")
    Te130_n2 = n2.get_isotope()
    Te130_n2.set_scintillator_masses(scintillator_cocktail, apply_fv_cut)
    Te130_n2.set_number_nuclei(Te130_n2.get_mass())
    n3 = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Te130_0_7_0-0_3-5.root")
    Te130_n3 = n3.get_isotope()
    Te130_n3.set_scintillator_masses(scintillator_cocktail, apply_fv_cut)
    Te130_n3.set_number_nuclei(Te130_n3.get_mass())
    n7 = Reconstructed\
        (os.environ.get("MAJORAT_DATA")+\
             "/hist_RAT4-5_10k_decay0_2beta_Te130_0_8_0-0_3-5.root")
    Te130_n7 = n7.get_isotope()
    Te130_n7.set_scintillator_masses(scintillator_cocktail, apply_fv_cut)
    Te130_n7.set_number_nuclei(Te130_n7.get_mass())
    if (args.mode == "limit_setting"):
        n1_limit_setter = SetLimit(sum_background_hist, n1, livetime)
        n1_limit_setter.set_t_half_limit(0.1e24, 1.0e27, 0.1e24, 2.71)
        n1_hist = n1_limit_setter.get_signal_hist()
        n2_limit_setter = SetLimit(sum_background_hist, n2, livetime)
        n2_limit_setter.set_t_half_limit(0.1e23, 1.0e26, 0.1e23, 2.71)
        n2_hist = n2_limit_setter.get_signal_hist()
        n3_limit_setter = SetLimit(sum_background_hist, n3, livetime)
        n3_limit_setter.set_t_half_limit(0.1e23, 1.0e26, 0.1e23, 2.71)
        n3_hist = n3_limit_setter.get_signal_hist()
        n7_limit_setter = SetLimit(sum_background_hist, n7, livetime)
        n7_limit_setter.set_t_half_limit(1.0e23, 2.0e25, 0.2e22, 2.71)
        n7_hist = n7_limit_setter.get_signal_hist()
        sum_hist = sum_background_hist + \
            n1_hist +\
            n2_hist+\
            n3_hist+\
            n7_hist
        sum_hist.SetTitle("Sum, 90 % CLs")
    else:
        Te130_n1.set_half_life("2.6e24 y")
        Te130_n1.set_counts(livetime)
        n1.scale_histogram(Te130_n1.get_counts())
        n1_hist = n1.get_histogram()
        n1_hist.SetDirectory(0)
        Te130_n2.set_half_life("1.0e24 y")
        Te130_n2.set_counts(livetime)
        n2.scale_histogram(Te130_n2.get_counts())
        n2_hist = n2.get_histogram()
        n2_hist.SetDirectory(0)
        Te130_n3.set_half_life("4.5e23 y")
        Te130_n3.set_counts(livetime)
        n3.scale_histogram(Te130_n3.get_counts())
        n3_hist = n3.get_histogram()
        n3_hist.SetDirectory(0)
        Te130_n7.set_half_life("1.1e22 y")
        Te130_n7.set_counts(livetime)
        n7.scale_histogram(Te130_n7.get_counts())
        n7_hist = n7.get_histogram()
        n7_hist.SetDirectory(0)
        sum_hist = sum_background_hist + \
            n1_hist +\
            n2_hist+\
            n3_hist+\
            n7_hist
        sum_hist.SetTitle("Sum")
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
    stack.Add(n1_hist)
    stack.Add(n2_hist)
    stack.Add(n3_hist)
    stack.Add(n7_hist)
    legend.AddEntry(n1_hist, n1_hist.GetTitle(), "l")
    legend.AddEntry(n2_hist, n2_hist.GetTitle(), "l")
    legend.AddEntry(n3_hist, n3_hist.GetTitle(), "l")
    legend.AddEntry(n7_hist, n7_hist.GetTitle(), "l")
    stack.Add(sum_hist)
    legend.AddEntry(sum_hist, sum_hist.GetTitle(), "l")
    ############################################################################

    # Spectral plot ############################################################
    if (args.mode == "limit_setting"):
        plots.make_majoron_spectral_plot("SNO+_90CL_Majoron_spectral_plot", stack, legend)
    else:
        plots.make_majoron_spectral_plot("SNO+_Majoron_spectral_plot", stack, legend)
