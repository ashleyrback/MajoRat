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
###########################################################################
import constants

import math

class DoubleBeta(object):
    """ Base class designed as a utility to handle an isotope's double beta
    decay properties, e.g. lifetime.
    """
    def __init__(self, isotope="Te130"):
        """ Initialised for isotope name supplied (Te130 by default), with
        some useful nuclear factors defined.
        """
        self._phase_space = constants.phase_space.get(isotope).get("2nu")
        self._matrix_element = constants.matrix_element.get(isotope).\
            get("2nu")
    def get_half_life(self):
        """ Method to return the two neutrino half life (yrs) """
        t_half = 1.0 / (self._phase_space*math.pow(self._matrix_element, 2))
        return t_half
class ZeroNuConverter(object):
    """ Base class designed as a utility to handle conversion between 
    different isotope properties, e.g. converting from lifetime to 
    effective double beta neutrino mass, or vice versa.
    """
    def __init__(self, isotope="Te130"):
        """ Initialised for isotope name supplied (Te130 by default), with
        some useful nuclear factors and constants defined.
        """
        two_nu = DoubleBeta(isotope)
        self._two_nu_t_half = two_nu.get_half_life();
        self._t_half_min = self._two_nu_t_half*1.0e-6 # Assumed lower limit
        self._phase_space = constants.phase_space.get(isotope).get("0nu")
        self._matrix_element = constants.matrix_element.get(isotope).get("0nu")
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
    def half_life_to_mass(self, t_half):
        """ Method to convert half life (Yr) to OvBB effective mass. 
        Returns the effective neutrino mass in eV.
        """
        if (self._conversion_factor == None):
            self.set_conversion_factor()
        try:
            assert (t_half > self._t_half_min),\
                "half life is less than 1.0e-6*two_nu_t_half"
        except AssertionError as detail:
            print "ZeroNuConverter.half_life_to_mass: WARNING", detail
            t_half = self._t_half_min
            print " --> setting half life to t_half_min", self._t_half_min
        effective_mass = self._conversion_factor * math.sqrt(1.0/t_half)
        return effective_mass
    def mass_to_half_life(self, effective_mass):
        """ Method to convert effective neutrino mass in eV to half life.
        Returns the half life (Yr)
        """
        if (self._conversion_factor == None):
            self.set_conversion_factor()
        try:
            t_half = math.pow(self._conversion_factor/effective_mass, 2)
        except ZeroDivisionError as detail:
            print "ZeroNuConverter.half_life_to_mass: WARNING", detail
            t_half = self._t_half_min
            print " --> setting half life to t_half_min", self._t_half_min    
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

    # Experiment with zero limits
    te_effective_mass = te_converter.half_life_to_mass(0.0)
    t_half_te = te_converter.mass_to_half_life(0.0)
    print "Te130: half life = " + str(t_half_te) + ", effective mass = "\
        + str(te_effective_mass)
    
