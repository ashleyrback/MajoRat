#!/usr/bin/env python
#
# make_python_header.py
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate energy spectrum "
                                     "for the SNO+ production ntuple file " 
                                     "specified and write histogram to file")
    parser.add_argument("ntuple_file", help="path to SNO+ production "
                        "ntuple file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-a", "--append", help="append to existing histogram",
                       action="store_true")
    group.add_argument("-r", "--always_remake", help="always generate new "
                       "histogram from ntuple (even if saved histogram exists",
                       action="store_true")
    args = parser.parse_args()

    spectrum = Production(args.ntuple_file)
    spectrum.set_parameters()
    spectrum.make_histogram(args.append, args.always_remake)
    spectrum.write_histogram()
