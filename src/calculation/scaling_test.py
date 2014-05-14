#!/usr/bin/env python
#
# scaling_test.py
#
# Script which checks that the histograms are being scaled correctly
#
# Author A R Back 
#
# 30/04/2014 <ab571@sussex.ac.uk> : First revision
###########################################################################
if __name__ == "__main__":
    from ROOT import TCanvas
    from ROOT import TH1D

    from reconstructed import Reconstructed
    import constants
    import statistics
    
    import argparse
    parser = argparse.ArgumentParser(description=("Set a 90% confidence limit "
                                                  "based on supplied signal "
                                                  "file and background file"))
    parser.add_argument("signal_file", help="path to RAT-generated Root file")
    parser.add_argument("background_file", help="path to RAT-generated Root "
                        "file")
    args = parser.parse_args()

    # Background
    t_half = constants.half_life.get("Xe136").get(5) # KLZ limit for SI=5
    background = Reconstructed(args.background_file, t_half)
    background_hist = background.get_histogram()

    # Signal
    t_half = constants.half_life.get("Xe136").get(0) # KLZ limit for 0vBB
    signal = Reconstructed(args.signal_file, t_half)
    signal_hist = signal.get_histogram()

    c1 = TCanvas()
    c1.cd()
    signal_hist.Draw()
    background_hist.Draw("same")
    c1.Print("unscaled.png")
    raw_input("RETURN to exit")

    background.scale_histogram()
    background_hist = background.get_histogram()

    signal.scale_histogram()
    signal_hist = signal.get_histogram()

    c2 = TCanvas()
    c2.cd()
    signal_hist.Draw()
    background_hist.Draw("same")
    c2.Print("scaled.png")
    raw_input("RETURN to exit")
    
 
