#!/usr/bin/env python
#
# error_utils.py
#
# Collection of methods for common error handling constructions
#
# Author A R Back 
#
# 18/07/2014 <ab571@sussex.ac.uk> : First revision
###########################################################################
import sys

def check_exists_in_dict(dict_name, dict, key, location_text):
    """ Check the contents of a dict to see if it contains the key supplied

    :param dict_name: name of dictionary
    :type dict_name: str
    :param dict: dict to search
    :type dict: dict
    :param key: key to search for in dict
    :type key: str
    :param location_text: text to include in error message to identify location
    :type location_text: str
    """
    try:
        assert (dict.get(key) != None),\
            key + " was not found in dictionary " + dict_name
    except AssertionError as detail:
        print location_text + ": error - ", detail
        raise

