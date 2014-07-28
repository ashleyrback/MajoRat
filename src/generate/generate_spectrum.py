#!/usr/bin/env python
#
# generate_spectrum.py
#
# Script to run a rat macro with command line arguments, based on a template
# macro
#
# Author A R Back 
# 31/01/2014 <ab571@sussex.ac.uk> : First revision
# 21/03/2014 <ab571@sussex.ac.uk> : Added inheritance from SpectrumData 
#                                   class
###########################################################################
from ROOT import TFile

from spectrum_data import SpectrumData
import file_manips

import re
import subprocess
import socket
import sys

line_to_edit = re.compile(r"^#/.*") # line beginning with "#/"

class GenerateSpectrum(SpectrumData):
    """ Class to generate RAT montecarlo for use with other scripts that
    inherit from the SpectrumData class.
    """
    def __init__(self, path, template):
        """ Use __init__ method of SpectrumData class. GenerateSpectrum class is 
        initiated by supplying the full path to the root file that you wish to
        generate.
        """
        super(GenerateSpectrum, self).__init__(path)
        self._template = template
    def read_template(self):
        """ Reads lines from template macro. Returns lines in a list. """
        template_file = open(self._template)
        template = template_file.readlines()
        template_file.close()
        return template
    def edit_template(self, template):
        """ Edits template macro based on parameters from the Root filename
        used to initialise the class.
        """
        macro = []
        for line in template:
            search_result = re.search(line_to_edit, line)
            if search_result != None:
                line = line[(search_result.start()+1):]
                if line.find("/rat/procset") >= 0:
                    parts = line.split(" ")
                    new_line = parts[0] + " " + parts[1] + " \""
                    new_line += self._path
                    new_line += "\"\n"
                    line = new_line
                if line.find("/generator/add") >= 0:
                    parts = line.split(" ")
                    new_line = parts[0] + " " + parts[1] + " "
                    parts2 = parts[2].split(":")
                    new_line += self._generator.get_generator() + ":" + parts2[1]
                    new_line += "\n"
                    line = new_line
                if line.find("/generator/vtx/set") >= 0:
                    parts = line.split(" ")
                    new_line = parts[0] + " "
                    new_line += self._generator.get_type() + " "
                    new_line += self._generator.get_isotope().get_name() + " "
                    new_line += str(int(self._generator.get_level())) + " "
                    new_line += str(int(self._generator.get_mode())) + " "
                    new_line += str(self._generator.get_e_lo()) + " "
                    new_line += str(self._generator.get_e_hi())
                    new_line += "\n"
                    line = new_line
                    try:
                        assert(re.match\
                                   (r"^.*\s[0-9]{1,2}?\s[0-9]{1,2}?\s"
                                    "[0-9]+\.[0-9]+\s[0-9]+\.[0-9]+$", line)\
                                   != None),\
                                   "/generator/vtx/set line has incorrect format"
                    except AssertionError as detail:
                        print "generate_spectrum.edit_template: error,", detail
                        sys.exit(1)
                if line.find("/rat/run/start") >= 0:
                    parts = line.split(" ")
                    new_line = parts[0] + " " + str(self._n_events)
                    new_line += "\n"
                    line = new_line
            macro.append(line)
        return macro
    def write_macro(self, macro, mac_dir=""):
        """ Writes macro based on edited macro template """
        if (mac_dir == ""):
            mac_dir = os.getcwd()
        self._mac_dir = mac_dir+"/"
        self._mac_path = self._mac_dir+self._name+".mac"
        mac_file = open(self._mac_path, "w")
        for line in macro:
            mac_file.write(line)
        mac_file.close()
    def run_rat(self):
        """ Runs rat using the macro created in write_macro """
        try:
            assert (self._mac_path != None), \
                "method GenerateSpectrum.write_macro must be used before " \
                "GenerateSpectrum.run_rat"
            command = "rat "+self._mac_path
            os.environ["RATROOT"]
            process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
            output = process.communicate()[0]
        except AssertionError as detail:
            print "GenerateSpectrum.run_rat: error cannot locate macro,", detail
            sys.exit(1)
        except KeyError as detail:
            print "GenerateSpectrum.run_rat: error", detail, "not set"
            print " --> source correct environment scripts before running!"
            sys.exit(1)
    def clean_up(self):
        """ Move DQ outputs to their appropriate directory """
        try:
            data_dir = os.environ["DATA"]
            plots_dir = os.environ["PLOTS"]
            logs_dir = os.environ["LOGS"]
        except KeyError as detail:
            print "GenerateSpectrum.clean_up: error", detail, "not set"
            print " --> source analysis environment scripts before running!"
            sys.exit(1)
        for root, dirs, files in os.walk(os.getcwd()):
            for file in files:
                is_data = re.search(r".*\.root$", file)
                is_plot = re.search(r".*\.png$", file)
                hostname = socket.gethostname()
                is_log =  re.search(r"^rat\."+hostname+r"\.[0-9]+\.log$", file)
                if is_data:
                    try:
                        root_file = TFile(file)
                        tree = root_file.Get("T")
                        tree.ls()
                    except ReferenceError as detail:
                        "generate_spectrum.clean_up: error in TFile,", detail
                        sys.exit(1)
                    file_manips.copy_file(os.path.join(root, file), data_dir)
                elif is_plot:
                    file_manips.copy_file(os.path.join(root, file), plots_dir)
                elif is_log:
                    file_manips.copy_file(os.path.join(root, file), logs_dir)

if __name__=="__main__":
    import argparse
    import os
    import contextlib
    import tempfile
    import shutil

    parser = argparse.ArgumentParser()
    parser.add_argument("rat_mc_path", 
                        help="path rat generated root file you wish to create")
    parser.add_argument("template_path", 
                        help="supply the path a rat macro template")
    args = parser.parse_args()
    print args

    # set environment
    env=os.environ.copy()

    # make temporary directory to write macro files
    @contextlib.contextmanager
    def temporary_directory(*args, **kwargs):
        d = tempfile.mkdtemp(*args, **kwargs)
        try:
            yield d
        finally:
            shutil.rmtree(d)

    with temporary_directory() as temp_dir:
        generate_spectrum = GenerateSpectrum(args.rat_mc_path, 
                                         args.template_path)
        template = generate_spectrum.read_template()
        macro = generate_spectrum.edit_template(template)
        generate_spectrum.write_macro(macro, temp_dir)
        command = "diff " + generate_spectrum._mac_path + " " \
            + args.template_path
        os.system(command)
        generate_spectrum.run_rat()
        generate_spectrum.clean_up()
