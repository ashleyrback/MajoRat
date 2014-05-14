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
def my_range(start, end, step):
    while start <= end:
        yield start
        start += step
def my_decreasing_range(start, end, step):
    while start >= end:
        yield start
        start -= step

if __name__ == "__main__":
    from ROOT import TCanvas
    from ROOT import TH1D
    from ROOT import TColor
    from ROOT import gROOT
    from ROOT import gPad
    from ROOT import TGaxis
    from ROOT import TAxis

    from reconstructed import Reconstructed
    import constants
    import statistics
    import physics
    
    import array
    import argparse
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
    print "set_limit.py:58 - converting to Te half life"
    print " --> Te half life = " + str(t_half_te)
    signal = Reconstructed(args.signal_file, t_half_te)
    signal.scale_by_mass(0)
    signal_hist = signal.get_histogram()

    # Define data histogram
    data_hist = background_hist + signal_hist
    data_hist.SetLineColor(2)
    data_hist.SetFillColor(2)

    label = "LL T_{1/2} limit setting " + signal._label
    ll_hist_t_half = TH1D(label, label, 20, 0.0, 1.0e25)
    label = "LL m_{#beta#beta} limit setting " + signal._label
    ll_hist_mass = TH1D(label, label, 20, 0.0, 1.0e25)

    label = "#DeltaLL T_{1/2} limit setting " + signal._label
    delta_ll_hist_t_half = TH1D(label, label, 20, 0.0, 1.0e25)
    label = "#DeltaLL m_{#beta#beta} limit setting " + signal._label
    delta_ll_hist_mass = TH1D(label, label, 20, 0.0, 1.0e25)

    # Redefine bins for mass histogram
    axis = ll_hist_mass.GetXaxis()
    bins = axis.GetNbins()
    new_bins = [0.0]
    for bin_ in range(2, bins+1):
        t_half = axis.GetBinLowEdge(bin_)
        print t_half
        new_bins.append(te_converter.half_life_to_mass(t_half))
    new_bins.sort()
    print new_bins
    new_bins_array = array.array("d", new_bins)
    print new_bins_array
    print len(new_bins_array)-1
    axis.Set(len(new_bins_array)-1, new_bins_array)

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

    # Loop over half life values
    for t_half in my_range(0.0, 1.0e25, 0.5e24):
        # Calculate effective mass
        effective_mass = te_converter.half_life_to_mass(t_half)

        # Get reconstructed histograms corresponding to t_half
        signal.scale_by_mass(effective_mass)
        signal_hist = signal.get_histogram()
        mc_hist = background_hist + signal_hist
        mc_hist.SetLineColor(4)
        c1 = TCanvas()
        c1.cd()
        data_hist.Draw()
        mc_hist.Draw("same")
        image_name = "data_vs_mc_" + str(t_half).replace(".","-") + ".png"
        c1.SetLogy()
        c1.Print(image_name)

        # From histograms get summed delta log-likelihood
        summed_ll_signal = statistics.sum_ll(data_hist, mc_hist)
        ll_hist_t_half.Fill(t_half, summed_ll_signal)
        ll_hist_mass.Fill(effective_mass, summed_ll_signal)

    # Half life
    c2 = TCanvas("c2", "c2", 1500, 600)
    TGaxis.SetMaxDigits(3)
    c2.Divide(2)
    c2.cd(1)
    ll_hist_t_half.Draw("l")
    ll_hist_t_half.SetXTitle("T_{1/2}")
    ll_hist_t_half.SetYTitle("LL")
    c2.cd(2)
    ll_hist_t_half.Draw("l")
    ll_hist_t_half.SetXTitle("T_{1/2}")
    ll_hist_t_half.SetYTitle("LL")
    gPad.SetLogy()
    image_name = "n" + str(signal._spectral_index) + "_ll_hist_t_half.png"
    c2.Print(image_name)

    # Delta log likelihood
    offset = ll_hist_t_half.GetMinimum()
    for bin_ in range(1, ll_hist_t_half.GetNbinsX()+1):
        log_likelihood = ll_hist_t_half.GetBinContent(bin_)
        delta_ll = log_likelihood - offset
        delta_ll_hist_t_half.SetBinContent(bin_, delta_ll)
    c3 = TCanvas("c3", "c3", 1500, 600)
    c3.Divide(2)
    c3.cd(1)
    delta_ll_hist_t_half.Draw("l")
    delta_ll_hist_t_half.SetXTitle("T_{1/2}")
    delta_ll_hist_t_half.SetYTitle("#DeltaLL")
    c3.cd(2)
    delta_ll_hist_t_half.Draw("l")
    delta_ll_hist_t_half.SetXTitle("T_{1/2}")
    delta_ll_hist_t_half.SetYTitle("#DeltaLL")
    gPad.SetLogy()
    image_name = "n" + str(signal._spectral_index) + "_delta_ll_hist_t_half.png"
    c3.Print(image_name)

    # Effective mass
    c4 = TCanvas("c4", "c4", 750, 600)
    c4.cd()
    ll_hist_mass.Draw("l")
    ll_hist_mass.SetXTitle("m_{#beta#beta}")
    ll_hist_mass.SetYTitle("LL")
    image_name = "n" + str(signal._spectral_index) + "_ll_hist_mass.png"
    c4.Print(image_name)

    # Delta log likelihood
    offset = ll_hist_mass.GetMinimum()
    for bin_ in range(1, ll_hist_mass.GetNbinsX()+1):
        log_likelihood = ll_hist_mass.GetBinContent(bin_)
        delta_ll = log_likelihood - offset
        delta_ll_hist_mass.SetBinContent(bin_, delta_ll)
    c5 = TCanvas("c5", "c5", 750, 600)
    c5.cd()
    delta_ll_hist_mass.Draw("l")
    delta_ll_hist_mass.SetXTitle("m_{#beta#beta}")
    delta_ll_hist_mass.SetYTitle("#DeltaLL")
    image_name = "n" + str(signal._spectral_index) + "_ll_hist_mass.png"
    c5.Print(image_name)
    raw_input("RETURN to exit")
