#!/usr/bin/env python
#
# detector_parameters.py
#
# Key (changeable) detector parameters - might move to a JSON file when I work
# out how to do this
#
# Author A R Back - 12/02/2014 <ab571@sussex.ac.uk> : First revision
###############################################################################
import math
from constants import density_LAB
from constants import Te_mass_fraction
from constants import n_avagadro
from constants import Te_atomic_weight

av_volume = 903.3 # m^3

def mass_Te(loading):
    """ Returns the mass of Te-130 in the AV, for a given natural Te loading 
    fraction
    """
    mass_LAB = density_LAB*av_volume
    mass = mass_LAB*Te_mass_fraction*loading
    print "mass = " + str(mass)
    return mass

def mass_Te_in_FV(loading, radius):
    """ Returns the mass of Te-130 in a FV - specified by radius - for a given
    natural Te loading fraction
    """
    volume = math.pow(radius, 3)
    volume *= math.pi
    volume *= (4.0/3.0)
    mass_LAB = density_LAB*volume
    mass = mass_LAB*Te_mass_fraction*loading
    print "mass = " + str(mass)
    return mass

def number_Te(mass_Te):
    """ Returns the number of Te-130 nuclei for a given mass (kg) of Te-130 """
    mass_Te *= 1e3 # convert mass to grams
    number = (mass_Te*n_avagadro)/Te_atomic_weight
    print "N_Te = " + str(number)
    return number

def n_decays(T_half, number):
    """ Returns the number of decays per year for a given half life and number
    of nuclei
    """
    decays = (math.log(2)*number)/T_half
    print "n_decays = " + str(decays)
    return decays
