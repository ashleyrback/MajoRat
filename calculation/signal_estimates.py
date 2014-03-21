#!/usr/bin/env python
#
# signal_esitmates.py
#
# Back-of-the-envelope-style calculations of expected signals
#
# Author A R Back - 12/02/2014 <ab571@sussex.ac.uk> : First revision
###############################################################################
from ROOT import TH1D

from spectrum_data import SpectrumData
from reconstructed import Reconstructed


class SignalEstimates(SpectrumData):
    """ Derived class, adapts spectrum data class for looking at estimated 
    signals and backgrounds.
    """
    def __init__(self, filepath):
        """ Initialise the container """
        super(SignalEstimates, self).__init__(filepath)
        self._events_above_2MeV = self.get_fraction_events(2.0)

    def get_fraction_events(self, elo=0.0, ehi=3.0):
        """ Integrates the normalised energy spectrum between the limits E_lo
        and E_hi. Returns the fraction of events within these limits.
        """
        spectrum = self.get_histogram()
        n_bins = spectrum.GetNbinsX()
        bin_width = spectrum.GetBinWidth(0)
        integrate_from = int(elo/bin_width)
        integrate_to = int(ehi/bin_width)
        integral = spectrum.Integral()
        try:
            assert (integral == self._events), ("SignalEstimates.get_fraction_"
                                                "events:\n --> spectrum"
                                                "integral not normalised")
        except AssertionError as detail:
            spectrum.Scale(self._events/spectrum.Integral())
            print "Handling run-time error: ", detail
        fraction = spectrum.Integral(integrate_from, integrate_to)
        return fraction

    def set_events_above_cut(self, cut_energy):
        events_above_cut = self.get_fraction_events(cut_energy)
        self._events_above_cut = (events_above_cut, cut_energy)
        return None

    def get_blank_histogram(self, title="", xaxis="", yaxis=""):
        """ Returns a blank histogram in Energy """
        label = self._label
        hist = TH1D(label,label, 300, 0.0, 3.0)
        hist.SetTitle(title)
        hist.SetXTitle(xaxis)
        hist.SetYTitle(yaxis)
        return hist

print __name__
if __name__ == "__main__":
    from ROOT import TH1D
    from ROOT import TGraph
    from ROOT import TCanvas

    from statistics import log_likelihood

    import math
    
    def my_range(start, end, step):
        while start <= end:
            yield start
            start += step

    filenames = ["hist_RAT4-5_1k_decay0_2beta_Te130_0_1_ELT3-5.root",
                 "hist_RAT4-5_1k_decay0_2beta_Te130_0_4_ELT3-5.root",
                 "hist_RAT4-5_1k_decay0_2beta_Te130_0_5_ELT3-5.root",
                 "hist_RAT4-5_1k_decay0_2beta_Te130_0_6_ELT3-5.root",
                 "hist_RAT4-5_1k_decay0_2beta_Te130_0_7_ELT3-5.root",
                 "hist_RAT4-5_1k_decay0_2beta_Te130_0_8_ELT3-5.root"]
    data_directory = "/home/ashley/snoplus/data/"
    filepath = data_directory + filenames[0]
    zero_nu = SignalEstimates(filepath)
    filepath = data_directory + filenames[1]
    two_nu = SignalEstimates(filepath)

    print zero_nu._label,zero_nu._t_half,zero_nu._events,zero_nu._events_above_2MeV
    print two_nu._label,two_nu._t_half,two_nu._events,two_nu._events_above_2MeV

    f = open("table_2nv0nu.tex", "w")
    f.write(r"\documentclass[18pt]{article}""\n")
    f.write(r"\usepackage{pdflscape}""\n")
    f.write(r"\begin{document}""\n")
    f.write(r"\begin{landscape}""\n")
    f.write(r"\begin{tabular}{r|c|c|c}""\n")
    f.write(r"Model & $T_{1/2}$ & Events & Events $\geq 2$ MeV \\ ""\n")
    f.write(r"\hline""\n")
    cell = " {0} & ".format(zero_nu._label)
    f.write(cell)
    cell = " {0:.1e} &".format(zero_nu._t_half)
    f.write(cell)
    cell = " {0:.2e} &".format(zero_nu._events)
    f.write(cell)
    cell = " {0:.0f} &".format(zero_nu._events_above_2MeV)
    f.write(cell)
    f.write(r"\\""\n")
    cell = " {0} & ".format(two_nu._label)
    f.write(cell)
    cell = " {0:.1e} &".format(two_nu._t_half)
    f.write(cell)
    cell = " {0:.2e} &".format(two_nu._events)
    f.write(cell)
    cell = " {0:.0f} &".format(two_nu._events_above_2MeV)
    f.write(cell)
    f.write(r"\\""\n")
    f.write(r"\end{tabular}""\n")
    f.write(r"\end{landscape}""\n")
    f.write(r"\end{document}""\n")
    f.close()

    rows = []
    for filename in filenames[2:]:
        filepath = data_directory + filename
        sig = SignalEstimates(filepath)
        bkg = two_nu
        
        # sig/sqrt(bkg) estimates
        sig_bkg = sig._events/math.sqrt(bkg._events)
        sig_bkg_above_2MeV = sig._events_above_2MeV/math.sqrt(bkg._events_above_2MeV)
        
        # Log-likelihood analysis
        sig_hist = sig.get_histogram()
        bkg_hist = bkg.get_histogram()
        assert (sig_hist.GetNbinsX() == bkg_hist.GetNbinsX()), \
            "Signal and background histograms differ in N_bins"
        ll_total = 0.0
        delta_ll_total = 0.0
#        ll_plot = TGraph()
#        ll_hist = sig.get_blank_histogram("#DeltaLL analysis", 
#                                          "E (MeV)", "#DeltaLL")
        for n_bin in range(1, sig_hist.GetNbinsX()):
            assert (sig_hist.GetBinLowEdge(n_bin) == bkg_hist.GetBinLowEdge(n_bin)), \
                "Bin energy differs between sig and bkg histograms"
            energy = sig_hist.GetBinLowEdge(n_bin)
            sig_events = sig_hist.GetBinContent(n_bin)
            bkg_events = bkg_hist.GetBinContent(n_bin)
            mc = bkg_events
            data = sig_events + bkg_events
            ll_s = math.fabs(log_likelihood(mc, data))
            data = bkg_events
            ll_0 = math.fabs(log_likelihood(mc, data))
            ll_total += ll_s
            delta_ll = ll_s - ll_0
#            ll_plot.SetPoint(n_bin, energy, delta_ll)
#            ll_hist.Fill(energy, delta_ll)
            delta_ll_total += delta_ll
        root_delta_ll_total = math.sqrt(delta_ll_total)
#        c1 = TCanvas()
#        c1.cd()
#        ll_plot.Draw("a*")
#        ll_plot.GetXaxis().SetTitle("E")
#        ll_plot.GetYaxis().SetTitle("#Deltall")
#        image_name = "n" + str(sig._spectral_index) + "_delta_ll_plot.png"
#        c1.Print(image_name)

#        c2 = TCanvas()
#        c2.cd()
#        ll_hist.Draw("p")
#        ll_hist.SetMarkerStyle(7)
#        ll_hist.SetMarkerSize(2.0)
#        image_name = "n" + str(sig._spectral_index) + "_delta_ll_hist.png"
#        c2.Print(image_name)
        
        row = ["n=" + str(sig._spectral_index),sig._t_half,sig._events,sig._events_above_2MeV,
               sig_bkg,sig_bkg_above_2MeV,delta_ll_total,root_delta_ll_total]
        rows.append(row)

        # Histogram sig/sqrt(bkg) as a function of E_cut
#        hist = sig.get_blank_histogram("Sig/\sqrt{Bkg} vs. E_{cut}",
#                                       "E_{cut} (MeV)",
#                                       "Sig/\sqrt{Bkg}")
#        sig_bkg_plot = TGraph()
#        for cut_energy in my_range(0.0, 3.0, 0.1):
#            sig.set_events_above_cut(cut_energy)
#            bkg.set_events_above_cut(cut_energy)
#            sig_events, sig_energy = sig._events_above_cut
#            bkg_events, bkg_energy = bkg._events_above_cut
#            try:
#                sig_over_root_bkg = sig_events/math.sqrt(bkg_events)
#                print ("E_cut = " + str(cut_energy) + ", sig/sqrt(bkg) = "
#                       + str(sig_over_root_bkg))
#                hist.Fill(cut_energy, sig_over_root_bkg)
#                sig_bkg_plot.SetPoint(int(cut_energy/0.1), cut_energy, 
#                                      sig_over_root_bkg) 
#            except ZeroDivisionError as detail:
#                print "Handling run-time error: ", detail

    f = open("table.tex", "w")
    f.write(r"\documentclass[18pt]{article}""\n")
    f.write(r"\usepackage{pdflscape}""\n")
    f.write(r"\begin{document}""\n")
    f.write(r"\begin{landscape}""\n")
    f.write(r"\begin{tabular}{r|c|c|c|c|c|c|c}""\n")
    f.write(r"Model & $T_{1/2}$ & Events & Events $\geq 2$ MeV & $\frac{S}{\sqrt{B}}$ & $\frac{S}{\sqrt{B}} \geq 2$ MeV & $\Delta LL$ & $\sqrt{\Delta LL}$ \\ ""\n")
    f.write(r"\hline""\n")
    for row in rows:
        cell = " {0} & ".format(row[0])
        f.write(cell)
        cell = " {0:.1e} &".format(row[1])
        f.write(cell)
        cell = " {0:.2e} &".format(row[2])
        f.write(cell)
        cell = " {0:.0f} &".format(row[3])
        f.write(cell)
        cell = " {0:4.2f} &".format(row[4])
        f.write(cell)
        cell = " {0:3.1f} &".format(row[5])
        f.write(cell)
        cell = " {0:4.2f} &".format(row[6])
        f.write(cell)
        cell = " {0:4.2f}".format(row[7])
        f.write(cell)
        f.write(r"\\""\n")
    f.write(r"\end{tabular}""\n")
    f.write(r"\end{landscape}""\n")
    f.write(r"\end{document}""\n")
    f.close()
