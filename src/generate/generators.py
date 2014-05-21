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

class Generator(object):
    """ Base class for handling generators """
    def __init__(self, isotope_, e_lo=0.0, e_hi=3.5):
        """ Initialise base class.

        :param isotope_: name of isotope simulated
        :type isotopr_: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        self._isotope = isotope_
        self._e_lo = e_lo
        self._e_hi = e_hi
        self._generator = None
        self._type = None
    def get_isotope(self):
        """        
        :returns: isotope
        :rtype: str (or isotope.Isotope or isotope.SNOPlusTe)
        """
        return self._isotope
    def get_e_lo(self):
        """
        :returns: lower KE limit on spectrum
        :rtype: float
        """
        return self._e_lo
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
    def __init__(self, isotope_, e_lo=0.0, e_hi=3.5):
        """ Initialise Decay0 class.

        :param isotope_: name of isotope simulated
        :type isotopr_: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        super(Decay0, self).__init__(isotope_, e_lo, e_hi)
        self._generator = "decay0"
class DecayChain(Generator):
    """ Derived class for handling decaychain generator """
    def __init__(self, isotope_, e_lo=0.0, e_hi=3.5):
        """ Initialise DecayChain class.

        :param isotope_: name of isotope simulated
        :type isotopr_: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        super(DecayChain, self).__init__(isotope_, e_lo, e_hi)
        self._generator = "decaychain"
class Solar(Generator):
    """ Derived class for handling solar generator """
    def __init__(self, isotope_, type_="nue", e_lo=0.0, e_hi=3.5):
        """ Initialise Solar class.

        :param isotope_: name of isotope simulated
        :type isotopr_: str
        :param type_: type of solar spectrum
        :type type_: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        super(Solar, self).__init__(isotope_, e_lo, e_hi)
        self._generator = "solar"
        self._type = type_
class Backg(Decay0):
    """ Derived class for handling backg generator """
    def __init__(self, isotope_, type_="backg", e_lo=0.0, e_hi=3.5):
        """ Initialise Backg class.
        
        :param isotope_: name of isotope simulated
        :type isotopr_: str
        :param type_: type of decay0 spectrum
        :type type_: str
        :param e_lo: lower limit on energy spectrum, if applicable
        :type e_lo: float
        :param e_hi: upper limit on energy spectrum, if applicable
        :type e_hi: float
        """
        super(Backg, self).__init__(isotope_, e_lo, e_hi)
        self._type = type_
class TwoBeta(Decay0):
    """ Derived class for handling backg generator """
    def __init__(self, isotope_, mode, type_="2beta", 
                 level=0, e_lo=0.0, e_hi=3.5):
        """ Initialise TwoBeta class.
        
        :param isotope_: isotope object
        :type isotope_: isotope.Isotope or isotope.SNOPlusTe
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
        super(TwoBeta, self).__init__(isotope_, e_lo, e_hi)
        self._type = type_
        self._level = level
        self._mode = mode
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

