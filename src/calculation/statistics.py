#!/usr/bin/env python
#
# statistics.py
#
# Some useful statistical methods
#
# Author A R Back
#
# 14/02/2014 <ab571@sussex.ac.uk> : First revision
# 30/04/2014 <ab571@sussex.ac.uk> : Added sum_ll method
###########################################################################
import math

def log_likelihood(data, mc):
    """ Returns the value of the log likelihood function for a given
    number of MC events and a given number of data events
    """
    try:
        ll = -2*(mc - data + data*math.log(data/mc))
    except ZeroDivisionError as detail:
        ll = 0.0
        print ("statistics.log_likelihood: runtime error:\n",
               detail)
    except ValueError as detail:
        ll = 0.0
        print ("statistics.log_likelihood: runtime error:\n",
               detail)
    return ll
def sum_ll(data_hist, mc_hist):
    """ Returns the summed log-likelihood for a given data spectrum,
    compared to a mc spectrum. "Data" should be a spectrum where T_half is
    unknown (i.e. cycling through), "mc" should be the spectrum where you 
    know T_half (fixed)
    """
    assert (data_hist.GetNbinsX() == mc_hist.GetNbinsX()), \
        "Signal and background histograms differ in N_bins"
    ll_signal_sum = 0.0
    for n_bin in range(1, data_hist.GetNbinsX()):
        energy = data_hist.GetBinLowEdge(n_bin)
        data_events = data_hist.GetBinContent(n_bin)
        mc_events = mc_hist.GetBinContent(n_bin)
        ll_signal = math.fabs(log_likelihood(data_events, mc_events))
        ll_signal_sum += ll_signal
    return ll_signal_sum
def sum_delta_ll(sig_hist, bkg_hist):
    """ Method to return the summed delta log-likelihood (ll) """
    assert (sig_hist.GetNbinsX() == bkg_hist.GetNbinsX()), \
        "Signal and background histograms differ in N_bins"
    delta_ll_total = 0.0
    for n_bin in range(1, sig_hist.GetNbinsX()):
        assert (sig_hist.GetBinLowEdge(n_bin) == bkg_hist.GetBinLowEdge(n_bin))\
            , "Bin energy differs between sig and bkg histograms"
        energy = sig_hist.GetBinLowEdge(n_bin)
        sig_events = sig_hist.GetBinContent(n_bin)
        bkg_events = bkg_hist.GetBinContent(n_bin)
        mc = bkg_events
        data = sig_events+bkg_events
        ll_mc = math.fabs(log_likelihood(mc, data))
        data = bkg_events
        ll_data = math.fabs(log_likelihood(mc, data))
        delta_ll = ll_mc-ll_data
        delta_ll_total += delta_ll
    return delta_ll_total
