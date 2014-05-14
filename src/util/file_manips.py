#!/usr/bin/env python
#
# file_manips.py
#
# A collection of useful methods for manipulating strings with respect to file
# paths and filenames
#
# Author A R Back - 31/01/2014 <ab571@sussex.ac.uk> : First revision
###############################################################################
import shutil

def cut_path(path):
    """ When supplied a filepath, returns a substring that is the just the 
    filename.
    """
    dir_ = path[:path.rfind("/")+1]
    filename = path[path.rfind("/")+1:]
    return (dir_, filename)

def strip_ext(filename):
    """ Returns a filename without its file extension e.g. myFile.cc would
    return myFile
    """
    s=filename
    sub=s[:s.find(".")]
    return sub

def split_ext(filename):
    """ As strip_ext, but returns the filename and extension in a tuple """
    s=filename
    name=s[:s.find(".")]
    ext=s[s.find("."):]
    return name, ext

def copy_file(source, destination):
    try:
        shutil.move(source, destination)
    except shutil.Error as detail:
        message = str(detail).split()
        if (message[2] == "'/home/ashley/Google"):
            true_dest = message[2]+message[3]
        else:
            true_dest = message[2]
        true_dest = true_dest[1:-1] # strip excess "'"
        file_name, ext = split_ext(true_dest)
        #print file_name, ext
        new_dest = file_name+"_#1"+ext
        #print new_dest
        print detail, "file not moved!"
        #print "writing to", new_dest
        #shutil.move(source, new_dest)
