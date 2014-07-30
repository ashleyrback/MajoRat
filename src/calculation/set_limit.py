#!/usr/bin/env python
#
# set_limit.py
#
# Limit setting script using summed delta log likelihood
#
# Author A R Back 
#
# 29/04/2014 <ab571@sussex.ac.uk> : First revision
###########################################################################
from ROOT import TCanvas
from ROOT import TH1D
from ROOT import TColor
from ROOT import gROOT
from ROOT import gPad
from ROOT import TGaxis
from ROOT import TAxis

from reconstructed import Reconstructed
import constants
import defaults
import statistics
import physics

import array
import argparse

def my_range(start, end, step):
    while start <= end:
        yield start
        start += step
def my_decreasing_range(start, end, step):
    while start >= end:
        yield start
        start -= step

class SetLimit(object):
    """ Base class for handling limit setting code """
    def __init__(self, sum_background_hist,
                 signal,
                 livetime="default"):
        """ Initialise SetLimit class. Note signal must be a full SpectrumData
        object so that it can be scaled properly.

        :param sum_background_hist: summed background for chi-squared fitting
        :type sum_background_hist: TH1D
        :param signal: signal spectrum to set limit for (add to background)
        :type signal: SpectrumData (Production)
        :param livetime: livetime to use for setting limit
        :type livetime: float
        """
        self._sum_background_hist = sum_background_hist
        self._signal = signal
        if (livetime == "default"):
            livetime = defaults.ll_analysis.get("livetime")
        self._livetime = livetime
        self._limit = None
        self._signal_hist = None
        self._sum_hist = None
    def set_mass_limit(self, required_delta_chi_squared="default",
                       outer_mass_lo="default",
                       outer_mass_hi="default",
                       outer_mass_step="default"):
        """ Outer loop. Iterrates over values of effective mass, producing a
        plot of delta-chi-squared vs. effective mass for each. The y-intercept
        (i.e. value of bin at 0 meV effective mass) is noted. The limit is the
        the last effective mass value before a delta-chi-squared that falls the
        required-delta-chi-squared threshold.

        :param required_delta_chi_squared: threshold for required CL
        :type required_delta_chi_squared: float
        :param outer_mass_lo: outer loop lower effective mass limit
        :type outer_mass_lo: float
        :param outer_mass_hi: outer loop upper effective mass limit
        :type outer_mass_hi: float
        :param outer_mass_step: step size to use in iteration
        :type outer_mass_step:
        """
        if (required_delta_chi_squared == "default"):
            required_delta_chi_squared = defaults.ll_analysis.get\
                ("delta_chi_squared")
        if (outer_mass_lo == "default"):
            outer_mass_lo = defaults.ll_analysis.get("outer_mass_lo")
        if (outer_mass_hi == "default"):
            outer_mass_hi = defaults.ll_analysis.get("outer_mass_hi")
        if (outer_mass_step == "default"):
            outer_mass_step = defaults.ll_analysis.get("outer_mass_step")
        no_limit_set = True
        for effective_mass in my_decreasing_range(outer_mass_hi,
                                                  outer_mass_lo,
                                                  outer_mass_step):
            hist_label = self._signal._label + "-reco_pos"
            self._signal.scale_by_mass(effective_mass, 
                                       hist_label, 
                                       self._livetime)
            signal_hist = self._signal.get_histogram(hist_label)
            signal_hist.SetDirectory(0)
            assert (type(signal_hist) is TH1D),\
                "SetLimit.set_limit: signal_hist is not of type TH1D"
            signal_hist.SetTitle(self._signal._label + "-" + \
                                     str(int(effective_mass*1000)) + "meV")
            sum_hist = self._sum_background_hist + signal_hist
            sum_hist.SetDirectory(0)
            assert (type(sum_hist) is TH1D),\
                "SetLimit.set_limit: sum_hist is not of type TH1D"
            sum_hist.SetTitle("Sum, " + self._signal._label + "-" + \
                                  str(int(effective_mass*1000)) + " meV")
            delta_chi_squared = self.get_delta_chi_squared_mass(sum_hist)
            if no_limit_set:
                if (delta_chi_squared >= required_delta_chi_squared):
                    self._limit = effective_mass * 1000.0
                    self.set_signal_hist(signal_hist)
                    self.set_sum_hist(sum_hist)
                else:
                    print "SetLimit.set_limit: set limit at " + str(self._limit) + "meV"
                    no_limit_set = False
    def set_t_half_limit(self,
                         t_half_lo="default",
                         t_half_hi="default",
                         t_half_step="default",
                         required_delta_chi_squared="default"):
        """
        :param t_half_lo: t_half lower limit
        :type t_half_lo: float
        :param t_half_hi: t_half upper limit
        :type t_half_hi: float
        :param t_half_step: t_half step
        :type t_half_step: float
        :param required_delta_chi_squared: e.g. (2.71 = 90 % CL,
                                                 25.0 = 5 sigma)
        :type required_delta_chi_squared: float
        :returns: t_half at limit
        :rtype: float
        """
        if (t_half_lo == "default"):
            t_half_lo = defaults.ll_analysis.get("inner_t_half_lo")
        assert (t_half_lo != 0.0),\
            "SetLimit.get_t_half_limit: error inner_t_half_lo == 0.0"
        if (t_half_hi == "default"):
            t_half_hi = defaults.ll_analysis.get("inner_t_half_hi")
        if (t_half_step == "default"):
            t_half_step = defaults.ll_analysis.get("inner_t_half_step")
        if (required_delta_chi_squared == "default"):
            required_delta_chi_squared = defaults.ll_analysis.get\
                ("delta_chi_squared")
        sum_hist = self._sum_background_hist # just background for the moment
        sum_hist.SetDirectory(0)
        assert (type(sum_hist) is TH1D),\
            "SetLimit.set_limit: sum_hist is not of type TH1D"
        t_half_limit = self.get_t_half_limit(sum_hist, 
                                             t_half_lo,
                                             t_half_hi,
                                             t_half_step,
                                             required_delta_chi_squared)
        self._limit = t_half_limit
        print "SetLimit.set_limit: set limit at " + str(self._limit) + " yr"
        raw_input("RETURN to continue...")
        signal_isotope = self._signal.get_isotope()
        signal_isotope.set_half_life(self._limit)
        signal_isotope.set_counts(self._livetime)
        self._signal.scale_histogram(signal_isotope.get_counts())
        signal_hist = self._signal.get_histogram()
        signal_hist.SetDirectory(0)
        assert (type(signal_hist) is TH1D),\
            "SetLimit.set_limit: signal_hist is not of type TH1D"
        # Format half life string
        if (int(self._limit/1.0e23) > 0):
            half_life_string = str(int(((self._limit/1.0e24)*100.0)+0.5)/100.0) 
            half_life_string += "*10^{24} y"
        elif (int(self._limit/1.0e20) > 0):
            half_life_string = str(int(((self._limit/1.0e21)*100.0)+0.5)/100.0) 
            half_life_string += "*10^{21} y"
        signal_hist.SetTitle(self._signal._label + " - " + half_life_string)
        self.set_signal_hist(signal_hist)
        sum_hist = sum_hist + signal_hist
        sum_hist.SetTitle("Sum, " + self._signal._label + "-" + \
                              str(int(self._limit/1.0e24)) + " 10^{24}yr")
        self.set_sum_hist(sum_hist)
    def set_signal_hist(self, signal_hist):
        """ Set signal_hist as class member, format.

        :param signal_hist: histogram to set as class member
        :type signal_hist: TH1D
        """
        self._signal_hist = signal_hist.Clone()
        self._signal_hist.SetDirectory(0)
        signal_hist_id = hex(id(signal_hist))
        self_signal_hist_id = hex(id(self._signal_hist))
        try:
            assert(signal_hist_id != self_signal_hist_id),\
                "clone of signal hist has same memory address as original"
        except AssertionError as detail:
            print "SetLimit.set_signal_hist: error", detail
            raise
        self._signal_hist.SetLineWidth(2)
        self._signal_hist.SetLineColor(1)
    def get_signal_hist(self):
        """
        :returns: signal histogram at limit
        :rtype: TH1D
        """
        assert (self._signal_hist != None),\
            "SetLimit.get_signal_hist: error signal_hist not set"
        return self._signal_hist
    def set_sum_hist(self, sum_hist):
        """ Set sum_hist as class member, format.

        :param sum_hist: histogram to set as class member
        :type sum_hist: TH1D
        """
        self._sum_hist = sum_hist.Clone()
        self._sum_hist.SetDirectory(0)
        sum_hist_id = hex(id(sum_hist))
        self_sum_hist_id = hex(id(self._sum_hist))
        try:
            assert(sum_hist_id != self_sum_hist_id),\
                "clone of signal hist has same memory address as original"
        except AssertionError as detail:
            print "SetLimit.set_sum_hist: error", detail
            raise
        self._sum_hist.SetLineColor(1)
        self._sum_hist.SetLineWidth(1)
        self._sum_hist.SetLineStyle(1)
    def get_sum_hist(self):
        """
        :returns: sum histogram at limit
        :rtype: TH1D
        """
        assert (self._sum_hist != None),\
            "SetLimit.get_sum_hist: error sum_hist not set"
        return self._sum_hist
    def get_delta_chi_squared_mass(self, data_hist,
                                   inner_mass_hi="default",
                                   inner_mass_step="default"):
        """ 
        :param data_hist: summed background plus scaled signal
        :type data_hist: TH1D
        :param inner_mass_hi: inner effective mass upper limit
        :type inner_mass_hi: float
        :param inner_mass_step: inner effective mass step
        :type inner_mass_step: float
        :returns: delta chi squared y-intercept (i.e. at 0 meV)
        :rtype: float
        """
        inner_mass_lo = defaults.ll_analysis.get("inner_mass_lo")
        assert (inner_mass_lo == 0.0),\
            "SetLimit.get_delta_chi_squared_mass: error inner_mass_lo != 0.0"
        if (inner_mass_hi == "default"):
            inner_mass_hi = defaults.ll_analysis.get("inner_mass_hi")
        if (inner_mass_step == "default"):
            inner_mass_step = defaults.ll_analysis.get("inner_mass_step")
        # Book histograms
        n_bins = int((inner_mass_hi-inner_mass_lo) / inner_mass_step)
        label = "Log-likelihood (m_{#beta#beta}) " + self._signal._label
        ll_hist = TH1D(label, label, n_bins, inner_mass_lo, inner_mass_hi)
        label = "#Delta#chi^{2} (m_{#beta#beta}) " + self._signal._label
        delta_chi_squared_hist = TH1D(label, label, n_bins, 
                                      inner_mass_lo, inner_mass_hi) 
        # Loop over effective_mass values
        for effective_mass in my_range(inner_mass_lo, 
                                       inner_mass_hi, 
                                       inner_mass_step):
            hist_label = self._signal._label + "-reco_pos"
            livetime = defaults.ll_analysis.get("livetime") 
            self._signal.scale_by_mass(effective_mass, hist_label, livetime)
            signal_hist = self._signal.get_histogram(hist_label)
            mc_hist = self._sum_background_hist + signal_hist
            assert (type(mc_hist) == TH1D),\
                "SetLimit.get_delta_chi_squared_mass: "\
                "error mc_hist not type TH1D"
            assert (mc_hist.Integral() > 0.0),\
                "SetLimit.get_delta_chi_squared_mass: "\
                "error mc_hist contains no events"
            assert (type(data_hist) == TH1D),\
                "SetLimit.get_delta_chi_squared_mass: "\
                "error data_hist not type TH1D"
            assert (data_hist.Integral() > 0.0),\
                "SetLimit.get_delta_chi_squared_mass: "\
                "error data_hist contains no events"
            # From histograms get summed delta log-likelihood
            summed_ll = statistics.sum_ll(data_hist, mc_hist)
            # Hack to remove binning errors
            effective_mass_shifted = effective_mass * (1+(inner_mass_step/0.1))
            ll_hist.Fill(effective_mass_shifted, summed_ll)
        # Delta log likelihood
        offset = ll_hist.GetMinimum()
        for bin_ in range(1, ll_hist.GetNbinsX()+1):
            log_likelihood = ll_hist.GetBinContent(bin_)
            delta_ll = log_likelihood - offset
            delta_chi_squared_hist.SetBinContent(bin_, 2.0*delta_ll)
        c1 = TCanvas("#Delta#chi^{2}", "#Delta#chi^{2} (m_{#beta#beta)",
                     750, 600)
        c1.cd()
        delta_chi_squared_hist.SetMarkerColor(9)
        delta_chi_squared_hist.SetMarkerStyle(7)
        delta_chi_squared_hist.SetMarkerSize(2.0)
        delta_chi_squared_hist.Draw("p")
        delta_chi_squared_hist.SetXTitle("m_{#beta#beta}")
        delta_chi_squared_hist.SetYTitle("#Delta#chi^{2}")
        image_name = data_hist.GetTitle() + "_delta_chi_squared_mass_hist.png"
        image_name = image_name.replace(",","")
        image_name = image_name.replace(" ","_")
        c1.Print(image_name)
        first_nonzero_bin = delta_chi_squared_hist.FindFirstBinAbove(0.0)
        return delta_chi_squared_hist.GetBinContent(first_nonzero_bin)
    def get_t_half_limit(self, data_hist,
                         inner_t_half_lo="default",
                         inner_t_half_hi="default",
                         inner_t_half_step="default",
                         required_delta_chi_squared="default"):
        """ 
        :param data_hist: summed background plus scaled signal
        :type data_hist: TH1D
        :param inner_t_half_hi: inner t_half upper limit
        :type inner_t_half_hi: float
        :param inner_t_half_step: inner t_half step
        :type inner_t_half_step: float
        :returns: t_half at limit
        :rtype: float
        """
        if (inner_t_half_lo == "default"):
            inner_t_half_lo = defaults.ll_analysis.get("inner_t_half_lo")
        assert (inner_t_half_lo != 0.0),\
            "SetLimit.get_t_half_limit: error inner_t_half_lo == 0.0"
        if (inner_t_half_hi == "default"):
            inner_t_half_hi = defaults.ll_analysis.get("inner_t_half_hi")
        if (inner_t_half_step == "default"):
            inner_t_half_step = defaults.ll_analysis.get("inner_t_half_step")
        if (required_delta_chi_squared == "default"):
            required_delta_chi_squared = defaults.ll_analysis.get\
                ("delta_chi_squared")
        # Book histograms
        n_bins = int((inner_t_half_hi-inner_t_half_lo) / inner_t_half_step)
        label = "LL (T_{1/2}) " + self._signal._label
        ll_hist = TH1D(label, label, (n_bins+1), inner_t_half_lo,
                       inner_t_half_hi)
        label = "#Delta#chi^{2} (T_{1/2}) " + self._signal._label
        delta_chi_squared_hist = TH1D(label, label, (n_bins+1), 
                                      inner_t_half_lo, inner_t_half_hi)
        assert (ll_hist.GetNbinsX() == delta_chi_squared_hist.GetNbinsX()),\
            "SetLimit.get_t_half_limit: error - number of bins mismatch"\
            " --> ll_hist: " + str(ll_hist.GetNbinsX()) + " "\
            "delta_chi_squared_hist: " + str(delta_chi_squared_hist.GetNbinsX())
        ll_axis = ll_hist.GetXaxis()
        delta_chi_axis = delta_chi_squared_hist.GetXaxis()
        bins = ll_axis.GetNbins()
        new_bins = []
        for bin_ in range(1, bins+1):
            t_half = ll_axis.GetBinLowEdge(bin_)
            new_bins.append(1.0/t_half)
        new_bins.sort()
        new_bins_array = array.array("d", new_bins)
        ll_axis.Set(len(new_bins_array)-1, new_bins_array)
        delta_chi_axis.Set(len(new_bins_array)-1, new_bins_array)
        assert (ll_hist.GetNbinsX() == delta_chi_squared_hist.GetNbinsX()),\
            "SetLimit.get_t_half_limit: error - number of bins mismatch"\
            " --> ll_hist: " + str(ll_hist.GetNbinsX()) + " "\
            "delta_chi_squared_hist: " + str(delta_chi_squared_hist.GetNbinsX())
        #assert (ll_hist.GetBinWidth(1) == inner_t_half_step),\
        #    "SetLimit.get_t_half_limit: error - wrong bin width\n"\
        #    " --> bin width = " + str(ll_hist.GetBinWidth(1)/1.0e23)\
        #    + "e23 not " + str(inner_t_half_step/1.0e23) + "e23"
        #assert (delta_chi_squared_hist.GetBinWidth(1) == inner_t_half_step),\
        #    "SetLimit.get_t_half_limit: error - wrong bin width"\
        #    " --> bin width = " + str(delta_chi_squared_hist.GetBinWidth(1)/1.0e23)\
        #    + "e23 not " + str(inner_t_half_step/1.0e23) + "e23"
        #assert (positive_delta_chi_hist.GetBinWidth(1) == inner_t_half_step),\
        #    "SetLimit.get_t_half_limit: error - wrong bin width"\
        #    " --> bin width = " + str(positive_delta_chi_hist.GetBinWidth(1)/1.0e23)\
        #    + "e23 not " + str(inner_t_half_step/1.0e23) + "e23"
        # Loop over t_half values
        for t_half in my_range(inner_t_half_lo, 
                               inner_t_half_hi, 
                               inner_t_half_step):
            #hist_label = self._signal._label + "-reco_pos"
            signal_isotope = self._signal.get_isotope()
            signal_isotope.set_half_life(t_half)
            signal_isotope.set_counts(self._livetime)
            self._signal.scale_histogram(signal_isotope.get_counts())
            signal_hist = self._signal.get_histogram()
            mc_hist = self._sum_background_hist + signal_hist
            assert (type(mc_hist) == TH1D),\
                "SetLimit.get_t_half_limit: "\
                "error mc_hist not type TH1D"
            assert (mc_hist.Integral() > 0.0),\
                "SetLimit.get_t_half_limit: "\
                "error mc_hist contains no events"
            assert (type(data_hist) == TH1D),\
                "SetLimit.get_t_half_limit: "\
                "error data_hist not type TH1D"
            assert (data_hist.Integral() > 0.0),\
                "SetLimit.get_t_half_limit: "\
                "error data_hist contains no events"
            # From histograms get summed delta log-likelihood
            summed_ll = statistics.sum_ll(data_hist, mc_hist)
            # Hack to remove binning errors
            t_half_shifted = t_half + (inner_t_half_step*0.1)
            ll_hist.Fill(1.0/t_half_shifted, summed_ll)
        # Delta log likelihood
        y_offset = ll_hist.GetMinimum()
        for bin_ in range(1, int(ll_hist.GetEntries())+1):
            log_likelihood = ll_hist.GetBinContent(bin_)
            delta_chi_squared = (log_likelihood-y_offset)
            delta_chi_squared_hist.SetBinContent(bin_, delta_chi_squared)
        assert (ll_hist.GetEntries() == delta_chi_squared_hist.GetEntries()),\
            "SetLimit.get_t_half_limit: error - number of entries mismatch"\
            " --> ll_hist: " + str(ll_hist.GetEntries()) + " "\
            "delta_chi_squared_hist: " + str(delta_chi_squared_hist.GetEntries())
        c1 = TCanvas("#Delta#chi^{2}", "#Delta#chi^{2} (m_{#beta#beta})",
                     800, 600)
        c1.cd()
        delta_chi_squared_hist.SetMarkerColor(9)
        delta_chi_squared_hist.SetMarkerStyle(7)
        delta_chi_squared_hist.SetMarkerSize(2.0)
        delta_chi_squared_hist.Draw("p")
        delta_chi_squared_hist.SetXTitle("T_{1/2}")
        delta_chi_squared_hist.SetYTitle("#Delta#chi^{2}")
        image_name = signal_hist.GetName() + "_delta_chi_squared_vs_t_half_hist.png"
        image_name = image_name.replace(" ","_")
        image_name = image_name.replace("(","_")
        image_name = image_name.replace(") ","_")
        c1.Print(image_name)
        first_bin_above = delta_chi_squared_hist.FindFirstBinAbove(2.71)
        try:
            assert (first_bin_above != -1),\
                "cannot set limit, no bin content above 2.71"
        except AssertionError as detail:
            print "SetLimit.get_t_half_limit: error -", detail
            raise
        raw_input("RETURN to continue...")
        return 1.0/delta_chi_squared_hist.GetBinLowEdge(first_bin_above)
    
    
                
def log_likelihood_half_life(data_hist,
                             sum_background_hist,
                             signal, # spectrum not histogram, so we can scale
                             t_half_lo="default",
                             t_half_hi="default",
                             t_half_step="default"):
    if (t_half_lo == "default"):
        t_half_lo = defaults.ll_analysis.get("t_half_lo")
    if (t_half_hi == "default"):
        t_half_hi = defaults.ll_analysis.get("t_half_hi")
    if (t_half_step == "default"):
        t_half_step = defaults.ll_analysis.get("t_half_step")
    # Book histograms
    n_bins = int((t_half_hi-t_half_lo) / t_half_step)
    label = "LL T_{1/2} limit setting " + signal._label
    ll_hist_t_half = TH1D(label, label, 20, 0.0, 1.0e25)
    label = "#DeltaLL T_{1/2} limit setting " + signal._label
    delta_ll_hist_t_half = TH1D(label, label, 20, 0.0, 1.0e25)
    c4 = TCanvas("t_half_scaling", "T_{1/2} scaling")
    data_hist_copy = data_hist.Clone()
    x_axis = data_hist_copy.GetXaxis()
    bin_width = x_axis.GetBinWidth(0) # assume uniform
    e_lo = 2.0 # MeV
    e_hi = 3.5 # MeV
    x_axis.SetRange(int(e_lo/bin_width), int(e_hi/bin_width))
    c4.cd()
    c4.SetLogy()
    data_hist_copy.Draw()

    # Loop over half life values
    for t_half in my_range(t_half_lo, t_half_hi, t_half_step):
        hist_label = "default"
        # Assume spectral plot default livetime 
        livetime = defaults.snoplus.get("livetime") 
        signal.scale_by_t_half(t_half, hist_label, livetime)
        signal_hist = signal.get_histogram(hist_label)
        mc_hist = background_hist + signal_hist
        mc_hist.SetLineColor(4)
        c1.cd()
        mc_hist.Draw("same")

        # From histograms get summed delta log-likelihood
        summed_ll_signal = statistics.sum_ll(data_hist, mc_hist)
        ll_hist_t_half.Fill(t_half, summed_ll_signal)
    image_name = signal_hist.GetTitle() + "_data_vs_mc.png"
    c4.Print(image_name)

    # Half life
    c5 = TCanvas("ll_t_half", "Log-likelihood T_{1/2}", 750, 600)
    TGaxis.SetMaxDigits(3)
    #c5.Divide(2)
    c5.cd()
    ll_hist_t_half.Draw("l")
    ll_hist_t_half.SetXTitle("T_{1/2}")
    ll_hist_t_half.SetYTitle("LL")
    #c5.cd(2)
    #ll_hist_t_half.Draw("l")
    #ll_hist_t_half.SetXTitle("T_{1/2}")
    #ll_hist_t_half.SetYTitle("LL")
    #gPad.SetLogy()
    image_name = signal_hist.GetTitle() + "_ll_hist_t_half.png"
    c5.Print(image_name)

    # Delta log likelihood
    offset = ll_hist_t_half.GetMinimum()
    for bin_ in range(1, ll_hist_t_half.GetNbinsX()+1):
        log_likelihood = ll_hist_t_half.GetBinContent(bin_)
        delta_ll = log_likelihood - offset
        delta_ll_hist_t_half.SetBinContent(bin_, delta_ll)
    c6 = TCanvas("#Deltall_t_half", "#DeltaLL T_{1/2}", 750, 600)
    #c6.Divide(2)
    c6.cd()
    delta_ll_hist_t_half.Draw("l")
    delta_ll_hist_t_half.SetXTitle("T_{1/2}")
    delta_ll_hist_t_half.SetYTitle("#DeltaLL")
    #c3.cd(2)
    #delta_ll_hist_t_half.Draw("l")
    #delta_ll_hist_t_half.SetXTitle("T_{1/2}")
    #delta_ll_hist_t_half.SetYTitle("#DeltaLL")
    #gPad.SetLogy()
    image_name = signal_hist.GetTitle() + "_delta_ll_hist_t_half.png"
    c6.Print(image_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=("Set a 90% confidence limit "
                                                  "based on supplied signal "
                                                  "file and background file"))
    parser.add_argument("signal_file", help="path to RAT-generated Root file")
    parser.add_argument("background_file", help="path to RAT-generated Root "
                        "file")
    args = parser.parse_args()
    
    # Background spectrum should not change so define here, based on KLZ half life
    # Assuming just the two neutrino background
    te_double_beta = physics.DoubleBeta()
    t_half_te_2nu = te_double_beta.get_half_life()
    background = Reconstructed(args.background_file, t_half_te_2nu)
    background.scale_by_t_half()
    background_hist = background.get_histogram()
    
    xe_converter = physics.ZeroNuConverter("Xe136")
    te_converter = physics.ZeroNuConverter()

    t_half_xe = constants.half_life.get("Xe136").get(0) # KLZ limit for 0vBB
    xe_effective_mass = xe_converter.half_life_to_mass(t_half_xe)
    t_half_te = te_converter.mass_to_half_life(xe_effective_mass)
    print "set_limit.py: - converting to Te half life"
    print " --> Te half life = " + str(t_half_te)
    signal = Reconstructed(args.signal_file, t_half_te)
    signal.scale_by_mass(0)
    signal_hist = signal.get_histogram()

    # Define data histogram
    data_hist = background_hist + signal_hist
    data_hist.SetLineColor(2)
    data_hist.SetFillColor(2)

    # Redefine bins for mass histogram
    """

    axis = delta_ll_hist_mass.GetXaxis()
    bins = axis.GetNbins()
    new_bins = [0.0]
    for bin_ in range(2, bins+1):
        t_half = axis.GetBinLowEdge(bin_)
        new_bins.append(te_converter.half_life_to_mass(t_half))
    new_bins.sort()
    print new_bins
    new_bins_array = array.array("d", new_bins)
    print new_bins_array
    print len(new_bins_array)-1
    axis.Set(len(new_bins_array)-1, new_bins_array)
    """
