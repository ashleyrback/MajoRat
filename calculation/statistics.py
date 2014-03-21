#!/usr/bin/env python
#
# statistics.py
#
# Some useful statistical methods
#
# Author A R Back - 14/02/2014 <ab571@sussex.ac.uk> : First revision
###############################################################################
import math

def log_likelihood(mc, data):
    """ Returns the value of the log likelihood function for a given
    number of MC events and a given number of data events
    """
    try:
        ll = -2*(mc - data + data*math.log(data/mc))
        print ll
    except ZeroDivisionError as detail:
        ll = 0.0
        print ("statistics.log_likelihood: runtime error:\n",
               detail)
    except ValueError as detail:
        ll = 0.0
        print ("statistics.log_likelihood: runtime error:\n",
               detail)
    return ll
