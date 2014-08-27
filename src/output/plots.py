#!/usr/bin/env python
#
# plots.py
#
# Module containing macros for all the standard plots output by MajoRat
#
# Author A R Back 
#
# 09/07/2014 <ab571@sussex.ac.uk> : First revision
########################################################################### 
import ROOT
from ROOT import TCanvas
from ROOT import TPad
from ROOT import gStyle

def make_spectral_plot(image_name,
                       stack,
                       legend):
    """
    :param image_name: name of output image (ext. ".png" will be added)
    :type image_name: str
    :param stack: stack of histograms to plot
    :type stack: THStack
    :param legend: legend to accompany histograms in stack
    :type legend: TLegend
    """
    ROOT.gROOT.SetStyle("Plain")
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadColor(0)
    gStyle.SetCanvasColor(0)
    gStyle.SetOptTitle(0)
    gStyle.SetLabelSize(0.06, "xyz")
    gStyle.SetTitleSize(0.06, "xyz")
    gStyle.SetOptStat("")
    c1 = TCanvas("canvas", "canvas", 1200, 900)
    c1.cd()
    p1 = TPad("p1", "Region of Interest", 0.0, 0.5, 0.7, 1.0, 0)
    p1.SetBottomMargin(0.15)
    p1.SetLeftMargin(0.15)
    p1.Draw()
    p2 = TPad("p2", "Key", 0.7, 0.5, 1.0, 1.0, 0)
    p2.SetBottomMargin(0.15)
    p2.SetTopMargin(0.15)
    p2.SetLeftMargin(0.05)
    p2.SetRightMargin(0.05)
    p2.Draw()
    p3 = TPad("p3", "Full Spectrum", 0.0, 0.0, 1.0, 0.5, 0)
    p3.SetBottomMargin(0.15)
    p3.SetLeftMargin(0.15)
    p3.Draw()
    p3.cd()
    gStyle.SetOptStat(0)
    p3.SetLogy()
    stack.Draw("nostack")
    stack.GetXaxis().SetTitle("E (MeV)")
    stack.GetYaxis().SetTitle("Events/(0.02 MeV)")
    stack.SetMinimum(0.1)
    p1.cd()
    gStyle.SetOptStat(0)
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
    image_path = image_name + ".png"
    c1.Print(image_path)
    raw_input("RETURN to exit")
def make_KLZ_spectral_plot(image_name,
                           stack,
                           legend):
    """
    :param image_name: name of output image (ext. ".png" will be added)
    :type image_name: str
    :param stack: stack of histograms to plot
    :type stack: THStack
    :param legend: legend to accompany histograms in stack
    :type legend: TLegend
    """
    ROOT.gROOT.SetStyle("Plain")
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadColor(0)
    gStyle.SetCanvasColor(0)
    gStyle.SetOptTitle(0)
    gStyle.SetLabelSize(0.06, "xyz")
    gStyle.SetTitleSize(0.06, "xyz")
    c1 = TCanvas("canvas", "canvas", 1200, 600)
    c1.cd()
    p1 = TPad("p1", "Full Spectrum", 0.0, 0.0, 1.0, 1.0, 0)
    p1.SetBottomMargin(0.15)
    p1.SetLeftMargin(0.15)
    p1.Draw()
    p2 = TPad("p2", "Key", 0.5, 0.6, 0.9, 0.9, 0)
    p2.SetBottomMargin(0.2)
    p2.SetTopMargin(0.2)
    p2.Draw()
    p2.cd()
    legend.SetTextSize(0.08)
    legend.SetNColumns(2)
    legend.SetEntrySeparation(0.01)
    legend.Draw()
    p1.cd()
    gStyle.SetOptStat(0)
    p1.SetLogy()
    stack.Draw("nostack")
    stack.GetXaxis().SetTitle("E (MeV)")
    stack.GetYaxis().SetTitle("Events/(0.05 MeV)")
    stack.SetMinimum(0.1)
    c1.Draw()
    image_path = image_name + ".png"
    c1.Print(image_path)
    raw_input("RETURN to exit")
def make_majoron_spectral_plot(image_name,
                               stack,
                               legend):
    """
    :param image_name: name of output image (ext. ".png" will be added)
    :type image_name: str
    :param stack: stack of histograms to plot
    :type stack: THStack
    :param legend: legend to accompany histograms in stack
    :type legend: TLegend
    """
    ROOT.gROOT.SetStyle("Plain")
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadColor(0)
    gStyle.SetCanvasColor(0)
    gStyle.SetOptTitle(0)
    gStyle.SetLabelSize(0.06, "xyz")
    gStyle.SetTitleSize(0.06, "xyz")
    c1 = TCanvas("canvas", "canvas", 1500, 600)
    c1.cd()
    p1 = TPad("p1", "Full Spectrum", 0.0, 0.0, 0.7, 1.0, 0)
    p1.SetBottomMargin(0.15)
    p1.SetLeftMargin(0.15)
    p1.Draw()
    p2 = TPad("p2", "Key", 0.7, 0.0, 1.0, 1.0, 0)
    p2.SetBottomMargin(0.2)
    p2.SetTopMargin(0.2)
    p2.Draw()
    p2.cd()
    legend.SetTextSize(0.05)
    legend.SetEntrySeparation(0.01)
    legend.Draw()
    p1.cd()
    gStyle.SetOptStat(0)
    p1.SetLogy()
    stack.Draw("nostack")
    stack.GetXaxis().SetTitle("E (MeV)")
    stack.GetYaxis().SetTitle("Events/(0.02 MeV)")
    stack.SetMinimum(0.1)
    c1.Draw()
    image_path = image_name + ".png"
    c1.Print(image_path)
    raw_input("RETURN to exit")
