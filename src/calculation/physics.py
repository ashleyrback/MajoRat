#!/usr/bin/env python
#
# physics.py
#
# Useful physics calculations related to 0vBB and Majoron modes
#
# Author A R Back 
#
# 05/05/2014 <ab571@sussex.ac.uk> : First revision
# 09/05/2014 <ab571@sussex.ac.uk> : Modified error handling in conversion
#                                   methods and added some get methods
# 21/05/2014 <ab571@sussex.ac.uk> : Valid range for half life and 
#                                   corresponding effective mass.
###########################################################################
import constants

import math
import sys

class DoubleBeta(object):
    """ Base class designed as a utility to handle an isotope's double beta
    decay properties, e.g. lifetime.
    """
    def __init__(self, isotope_name="Te130",
                 t_half_min_scaling=1.0e-6,  # Assumed lower limit 10^-6 * 2nu
                 t_half_max_scaling=1.0e15): # Assumed upper limit 10^15 * 2nu
        """ Initialised for isotope name supplied (Te130 by default), with
        some useful nuclear factors defined.
        
        :param isotope_name: name of isotope
        :type isotope_name: float
        :param t_half_min_scaling: scaling factor for lower limit on t_half
        :type t_half_min_scaling: float
        :param t_half_max_scaling: scaling factor for upper limit on t_half
        :type t_half_max_scaling: float
        """
        self._phase_space = constants.phase_spaces.get(isotope_name).get("2nu")
        self._matrix_element = constants.matrix_elements.get(isotope_name).\
            get("2nu")
        self._t_half_min = self.get_half_life() * t_half_min_scaling
        self._t_half_max = self.get_half_life() * t_half_max_scaling
    def get_half_life(self):
        """ 
        :returns: SM double beta decay half life (years)
        :rtype: float
        """
        t_half = 1.0 / (self._phase_space*math.pow(self._matrix_element, 2))
        return t_half
    def get_t_half_min(self):
        """
        :returns: lower bound on double beta decay half life
        :rtype: float
        """
        return self._t_half_min
    def get_t_half_max(self):
        """
        :returns: lower bound on double beta decay half life
        :rtype: float
        """
        return self._t_half_max
class ZeroNuConverter(object):
    """ Base class designed as a utility to handle conversion between 
    different isotope properties, e.g. converting from lifetime to 
    effective double beta neutrino mass, or vice versa.
    """
    def __init__(self, isotope_name="Te130",
                 t_half_min_scaling=1.0e-6,  # Assumed lower limit 10^-6 * 2nu
                 t_half_max_scaling=1.0e15): # Assumed upper limit 10^15 * 2nu
        """ Initialised for isotope name supplied (Te130 by default), with
        some useful nuclear factors and constants defined.
        """
        two_nu = DoubleBeta(isotope_name,
                            t_half_min_scaling,
                            t_half_max_scaling)
        self._t_half_min = two_nu.get_t_half_min()
        self._t_half_max = two_nu.get_t_half_max()
        self._phase_space = constants.phase_spaces.get(isotope_name).get("0nu")
        self._matrix_element = constants.matrix_elements.get(isotope_name).\
            get("0nu")
        self._coupling_constant = constants.coupling_constant
        self._conversion_factor = None
    def set_conversion_factor(self):
        """ Calculation of conversion factor, which is then set as a class
        member, this way it only needs to be calculated once for each class
        instance.
        """
        electron_mass = constants.electron_mass * 1.0e6 # converted to eV
        # Conversion factor
        conversion_factor = math.pow(electron_mass, 2)
        conversion_factor /= self._phase_space
        conversion_factor /= math.pow(self._coupling_constant, 4)
        conversion_factor /= math.pow(self._matrix_element, 2)
        conversion_factor = math.sqrt(conversion_factor)
        self._conversion_factor = conversion_factor
    def get_conversion_factor(self):
        """ Method to get conversion factor. If conversion factor is None,
        sets conversion factor. Returns conversion factor.
        """
        if (self._conversion_factor == None):
            self.set_conversion_factor()
        return self._conversion_factor
    def get_t_half_min(self):
        """ Method to get minimum half life. Returns t_half_min """
        return self._t_half_min
    def get_t_half_max(self):
        """ Method to get maximum half life. Returns t_half_max """
        return self._t_half_max
    def get_mass_min(self):
        """ Method to get minimum effective mass. Returns mass_min """
        return self.half_life_to_mass(self._t_half_max)
    def get_mass_max(self):
        """ Method to get maximum effective mass. Returns mass_max """
        return self.half_life_to_mass(self._t_half_min)
    def half_life_to_mass(self, t_half):
        """ Method to convert half life (Yr) to OvBB effective mass. 
        Returns the effective neutrino mass in eV.
        """
        if (self._conversion_factor == None):
            self.set_conversion_factor()
        if (t_half < self.get_t_half_min()):
            t_half = self.get_t_half_min()
        if (t_half > self.get_t_half_max()):
            t_half = self.get_t_half_max()
        try:
            assert (self.get_t_half_min() <= t_half <= self.get_t_half_max()),\
                "half life does not fall within the accepted range"
        except AssertionError as detail:
            print "ZeroNuConverter.half_life_to_mass: ERROR", detail
            print " --> cannont convert to mass"
            sys.exit(1)
        effective_mass = self._conversion_factor * math.sqrt(1.0/t_half)
        return effective_mass
    def mass_to_half_life(self, mass):
        """ Method to convert effective neutrino mass in eV to half life.
        Returns the half life (Yr)
        """
        if (self._conversion_factor == None):
            self.set_conversion_factor()
        if (mass < self.get_mass_min()):
            mass = self.get_mass_min()
        if (mass > self.get_mass_max()):
            mass = self.get_mass_max()
        try:
            assert (self.get_mass_min() <= mass <= self.get_mass_max()),\
                "effective mass does not fall within the accepted range"
        except AssertionError as detail:
            print "ZeroNuConverter.mass_to_half_life: ERROR", detail
            print " --> cannont convert to half_life"
            sys.exit(1)
        t_half = math.pow(self._conversion_factor/mass, 2)
        return t_half

if __name__ == "__main__":
    te_double_beta = DoubleBeta("Te130")
    t_half_te_2nu = te_double_beta.get_half_life()
    t_half_xe = constants.half_life.get("Xe136").get(0)
    xe_converter = ZeroNuConverter("Xe136")
    xe_effective_mass = xe_converter.half_life_to_mass(t_half_xe)
    te_converter = ZeroNuConverter("Te130")
    t_half_te = te_converter.mass_to_half_life(xe_effective_mass)
    te_effective_mass = te_converter.half_life_to_mass(t_half_te)
    print "Xe136: half life = " + str(t_half_xe) + ", effective mass = "\
        + str(xe_effective_mass)
    print "Te130: half life = " + str(t_half_te) + ", effective mass = "\
        + str(te_effective_mass)
    print "Te130: 2nu half life = " + str(t_half_te_2nu)

    # Experiment with limits
    effective_mass_max = te_converter.half_life_to_mass(0.0) # when t_half=0
    t_half_min = te_converter.mass_to_half_life(effective_mass_max)
    t_half_max = te_converter.mass_to_half_life(0.0) # when mass=0
    effective_mass_min = te_converter.half_life_to_mass(t_half_max)
    print "Te130: half life (min) = " + str(t_half_min) + ", effective mass (max) = "\
        + str(effective_mass_max)
    print "Te130: half life (max) = " + str(t_half_max) + ", effective mass (min) = "\
        + str(effective_mass_min)
