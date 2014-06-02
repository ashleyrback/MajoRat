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
import argapse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate energy spectrum "
                                     "for the SNO+ production ntuple file " 
                                     "specified and write histogram to file")
    parser.add_argument("ntuple_path", help="path to SNO+ production "
                        "ntuple file(s)")
    parser.add_argument("-d", "--is_directory", help="use to specify that "
                        "the path supplied is a directory", 
                        action="store_true")
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
        spectrum = Production(args.ntuple_file)
        spectrum.set_parameters()
        if args.is_directory: # default is to always_remake for first file
                              # then append for subsequent files
            if first_file:
                first_file = False
                append = False
                always_remake = True
                spectrum.make_histogram(append, always_remake)
            else:
                append = True
                spectrum.make_histogram(append)
        else: # use command line options set
            spectrum.make_histogram(args.append, args.always_remake)
        spectrum.write_histogram()
