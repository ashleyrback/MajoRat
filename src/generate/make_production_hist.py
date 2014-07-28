#!/usr/bin/env python
#
# make_production_hist.py
#
# Script to process production ntuples and generate spectra
#
# Author A R Back 
#
# 29/05/2014 <ab571@sussex.ac.uk> : First revision
#
###########################################################################
from production import Production

import argparse
import os
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate energy spectrum "
                                     "for the SNO+ production ntuple file " 
                                     "specified and write histogram to file")
    parser.add_argument("ntuple_path", help="path to SNO+ production "
                        "ntuple file(s)")
    parser.add_argument("-d", "--is_directory", help="use to specify that "
                        "the path supplied is a directory", 
                        action="store_true")
    parser.add_argument("-z", "--zero_energy", help="include zero energy "
                        "events", action="store_true")
    parser.add_argument("-p", "--reco_pos", help="use reconstructed position "
                        "in FV cut", action="store_true")
    parser.add_argument("-n", "--nhit_energy", help="use energy calculated "
                        "from nhits", action="store_true")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-a", "--append", help="append to existing histogram",
                       action="store_true")
    group.add_argument("-r", "--always_remake", help="always generate new "
                       "histogram from ntuple (even if saved histogram exists",
                       action="store_true")
    args = parser.parse_args()

    # make list of files
    file_list = []
    if args.is_directory:
        for root, dirs, files in os.walk(args.ntuple_path):
            for file in files:
                match = re.search(r".ntuple.root$", file)
                if match:
                    file_list.append(os.path.join(root, file))
    else:
        file_list.append(args.ntuple_path)
    first_file = True
    for file in file_list:
        spectrum = Production(file)
        spectrum.set_parameters()
        # make list of histogram labels
        hist_labels = [spectrum._label, 
                       spectrum._label+"-no_fv_cut"] # produced by default
        if args.reco_pos:
            hist_labels.append(spectrum._label+"-reco_pos")
        if args.nhit_energy:
            hist_labels.append(spectrum._label+"-nhit_energy")
            hist_labels.append(spectrum._label+"-reco_pos-nhit_energy")
        if args.zero_energy:
            zero_energy_labels = []
            for hist_label in hist_labels:
                hist_label += "-zero_energy"
                zero_energy_labels.append(hist_label)
            hist_labels += zero_energy_labels
        print hist_labels
        for hist_label in hist_labels:
            if args.is_directory: # default is to always_remake for first file
                                  # then append for subsequent files
                if first_file:
                    first_file = False
                    append = False
                    always_remake = True
                    spectrum.make_histogram(hist_label, append, always_remake)
                else:
                    append = True
                    spectrum.make_histogram(hist_label, append)
            else: # use command line options set
                spectrum.make_histogram(hist_label,
                                        args.append,
                                        args.always_remake)
        spectrum.write_histograms()
