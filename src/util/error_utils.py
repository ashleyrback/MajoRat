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
def almost_equal(value1, value2, precision):
    """ Checks that two values are equal to a given precision.

    :param value1: first value to check
    :type value1: float
    :param value2: second value to check
    :type value2: float
    :param precision: precision to which the two values should be equal.
                      Specify the last significant figure that should be
                      equal, e.g. to two decimal places --> 1.0e-2, to the
                      nearest thousand --> 1.0e3 etc.
    :type precision: float
    """
    # Convert value 1
    value1_rounded = int((value1/precision)+0.5) * precision
    value2_rounded = int((value2/precision)+0.5) * precision
    precision_rounded = int((precision/precision)+0.5) * precision
    assert (value1_rounded == value2_rounded),\
        "error_utils.almost_equal: error - supplied values are not equal to a "\
        "precision of " + str(precision) + "\n --> " + str(value1_rounded) + \
        " != " + str(value2_rounded)
