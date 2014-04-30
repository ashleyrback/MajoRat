#!/usr/bin/env python
#
# isotope.py
#
# Isotope class contains useful methods for working with isotopes:
#   - number of decays
#   - number of nuclei
# SNOPlusTe is a specialised version of this class that can calculate these
# things in the SNO+ framework (e.g. for LAB in the SNO+ AV)
#
# Author A R Back 
# 
# 27/04/2014 <ab571@sussex.ac.uk> : First revision 
#                                   (previously detector_parameters.py)
###############################################################################
import constants

import math
import sys

class Isotope(object):
    """ Base isotope class, for dealing with numbers of isotopes/decays
    for any isotope
    """
    def __init__(self, isotope_name):
        """ Initialise class based on the name of the isotope """
        self._name = isotope_name
    def get_number_nuclei(self, isotope_mass):
        """ Returns the number of isotope for a given mass (kg) """
        try:
            assert constants.atomic_weights.get(self._name) != None, \
                "Atomic weight " + self._name + " is not available"
        except AssertionError as detail:
            print "Isotope:get_number_nuclei error:", detail
            sys.exit(1)
        atomic_weight = constants.atomic_weights.get(self._name)
        isotope_mass *= 1e3 # convert mass to grams
        number_nuclei = (isotope_mass*constants.n_avagadro)/atomic_weight
        return number_nuclei
    def get_number_decays(self, t_half, number_nuclei):
        """ Returns the number of decays (per year) for a given half life
        (in years) and number of nuclei
        """
        number_decays = (math.log(2)*number_nuclei)/t_half
        return number_decays

### SNO+ SPECIFIC #############################################################
class SNOPlusTe(Isotope):
    """ Special case for SNO+ Te-130, derived from isotope class. Provides
    additional useful methods for getting mass of Te inside AV, or a given
    fiducial volume.
    """
    def __init__(self, isotope_name="Te130"):
        """ Initialised with Te130 as the default choice for isotope """
        super(SNOPlusTe, self).__init__(isotope_name)
        # Constants defined here are taken from the Te Verification Report,
        # unless stated otherwise stated
        self._av_volume = 903.3 # m^3
        self._density_lab = 0.862e3 # kgm^-3
    def get_mass(self, loading = 0.003):
        """ Returns the mass of Te-130 in the AV, for a given natural Te
        loading fraction. Assumes pure LAB scintillator.
        """
        mass_lab = self._density_lab*self._av_volume
        mass = mass_lab*constants.mass_fractions.get(self._name)*loading
        return mass
    def get_mass_in_fv(self, radius, loading = 0.003):
        """ Returns the mass of Te-130 in a FV - specified by radius - 
        for a given natural Te loading fraction
        """
        volume = math.pow(radius, 3)
        volume *= math.pi
        volume *= (4.0/3.0)
        mass_lab = self._density_lab*volume
        mass = mass_lab*constants.mass_fractions.get(self._name)*loading
        return mass
    def get_number_nuclei(self, radius, loading = 0.003):
        """ For a given loading fraction (not percentage), and radius, 
        which defines a fiducial volume, the number of nuclei that could 
        decay is returned
        """
        mass_in_fv = self.get_mass_in_fv(radius, loading)
        number_nuclei = super(SNOPlusTe, self).get_number_nuclei(mass_in_fv)
        return number_nuclei
###############################################################################
