#!/usr/bin/env python
#
# generators.py
#
# Classes to handle the different rat generators, when looking at 
# production data
#
# Author A R Back 
#
# 21/05/2014 <ab571@sussex.ac.uk> : First revision
#
###########################################################################
import isotope
import defaults

class Generator(object):
    """ Base class for handling generators """
    def __init__(self, isotope_name, e_lo="default", e_hi="default"):
        """ Initialise base class.

        :param isotope_name: name of isotope simulated
        :type isotope_name: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        if (e_lo == "default"):
            e_lo = defaults.spectrum.get("e_lo") 
        if (e_hi == "default"):
            e_hi = defaults.spectrum.get("e_hi")
        self.set_isotope(isotope_name)
        self._e_lo = e_lo
        self._e_hi = e_hi
        self._generator = None
        self._type = None
    def set_isotope(self, isotope_name):
        """ Create isotope.Isotope object for self._isotope

        :param isotope_name: name of isotope
        :type isotope_name: str
        """
        self._isotope = isotope.Isotope(isotope_name)
    def get_isotope(self):
        """ Gets isotope object
        
        :returns: self._isotope
        :rtype: isotope.Isotope (or derived object)
        """
        assert (self._isotope != None), "Generator.get_isotope: error - "\
            "isotope object does not exist" 
        return self._isotope
    def set_e_lo(self, e_lo):
        """ Sets lower KE limit for generated spectrum
        
        :param e_lo: lower KE limit
        :type e_lo: float
        """
        self._e_lo = e_lo
    def get_e_lo(self):
        """
        :returns: lower KE limit on spectrum
        :rtype: float
        """
        return self._e_lo
    def set_e_hi(self, e_hi):
        """ Sets lower KE limit for generated spectrum
        
        :param e_hi: upper KE limit
        :type e_hi: float
        """
        self._e_hi = e_hi
    def get_e_hi(self):
        """
        :returns: upper KE limit on spectrum
        :rtype: float
        """
        return self._e_hi
    def get_generator(self):
        """
        :returns: generator
        :rtype: str
        """
        return self._generator
    def get_type(self):
        """
        :returns: type
        :rtype: str
        """
        return self._type
class Decay0(Generator):
    """ Derived class for handling decay0 generator """
    def __init__(self, isotope_name, e_lo="default", e_hi="default" ):
        """ Initialise Decay0 class.

        :param isotope_name: name of isotope simulated
        :type isotope_name: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        super(Decay0, self).__init__(isotope_name, e_lo, e_hi)
        self._generator = "decay0"
class DecayChain(Generator):
    """ Derived class for handling decaychain generator """
    def __init__(self, isotope_name, e_lo="default", e_hi="default"):
        """ Initialise DecayChain class.

        :param isotope_name: name of isotope simulated
        :type isotope_name: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        super(DecayChain, self).__init__(isotope_name, e_lo, e_hi)
        self._generator = "decaychain"
    def set_isotope(self, isotope_name):
        """ Overloads set_isotope method for DecayChain
        
        :param isotope_name: name of isotope
        :type isotope_name: str
        """
        if (constants.half_lives.get("internal").get(isotope_name) != None):
            self._isotope = isotope.SNOPlusInternal(isotope_name)
        else:
            self._isotope = isotope.Isotope(isotope_name)
class Solar(Generator):
    """ Derived class for handling solar generator """
    def __init__(self, isotope_name, type_="nue", e_lo="default", e_hi="default"):
        """ Initialise Solar class.

        :param isotope_name: name of isotope simulated
        :type isotope_name: str
        :param type_: type of solar spectrum
        :type type_: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        super(Solar, self).__init__(isotope_name, e_lo, e_hi)
        self._generator = "solar"
        self._type = type_
class Backg(Decay0):
    """ Derived class for handling backg generator """
    def __init__(self, isotope_name, type_="backg", e_lo="default", e_hi="default"):
        """ Initialise Backg class.
        
        :param isotope_name: name of isotope simulated
        :type isotope_name: str
        :param type_: type of decay0 spectrum
        :type type_: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        super(Backg, self).__init__(isotope_name, e_lo, e_hi)
        self._type = type_
    def set_isotope(self, isotope_name):
        """ Overloads set_isotope method for Backg
        
        :param isotope_name: name of isotope
        :type isotope_name: str
        """
        if (constants.half_lives.get("internal").get(isotope_name) != None):
            self._isotope = isotope.SNOPlusInternal(isotope_name)
        else:
            self._isotope = isotope.Isotope(isotope_name)
class TwoBeta(Decay0):
    """ Derived class for handling backg generator """
    def __init__(self, isotope_name, mode, level=0,
                 type_="2beta", e_lo="default", e_hi="default"):
        """ Initialise TwoBeta class.
        
        :param isotope_name: isotope object
        :type isotope_name: isotope.Isotope or isotope.SNOPlusTe
        :param mode: decay0 2beta decay mode
        :type mode: float
        :param type_: type of decay0 spectrum
        :type type_: str
        :param level: decay0 2beta decay level
        :type level: float
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        self._type = type_
        self._level = level
        self._mode = mode
        super(TwoBeta, self).__init__(isotope_name, e_lo, e_hi)
    def set_isotope(self, isotope_name):
        """ Overloads set_isotope method for TwoBeta
        
        :param isotope_name: name of isotope
        :type isotope_name: str
        """
        if (self._mode == 1): # neutrinoless double beta decay
            self._isotope = isotope.ZeroNu(isotope_name)
        else:
            self._isotope = isotope.DoubleBetaIsotope(isotope_name)
        self._isotope.set_mode(self._mode)
    def get_level(self):
        """
        :returns: decay0 2beta level
        :rtype: float
        """
        return self._level
    def get_mode(self):
        """
        :returns: decay0 2beta mode
        :rtype: float
        """
        return self._mode

