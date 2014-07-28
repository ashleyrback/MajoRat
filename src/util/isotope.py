#!/usr/bin/env python
#
# isotope.py
#
# Isotope class contains useful methods for working with isotopes:
#   - number of nuclei
#   - number of counts
# Internal class extends this so you can calculate easily calculate the
# number of counts/nuclei for an isotope within the scintillator
# DoubleBeta class extends this further for looking specifically at a
# double beta isotope loaded into the scintillator
# ZeroNu class is a special case of DoubleBeta that for neutrinoless
# double beta decays. The key extension provided here is to set either a
# half life or effective mass
#
# Author A R Back 
# 
# 27/04/2014 <ab571@sussex.ac.uk> : First revision 
#                                   (previously detector_parameters.py)
# 09/05/2014 <ab571@sussex.ac.uk> : Added get_decays_from_mass method and 
#                                   modified get_decays_from_t_half method
# 21/05/2014 <ab571@sussex.ac.uk> : Valid range for half life and 
#                                   corresponding effective mass.
# 05/06/2014 <ab571@sussex.ac.uk> : Improved defaults
# 23/07/2014 <ab571@sussex.ac.uk> : Completely re-designed Isotope class
#                                   and derived classes so that different
#                                   types of backgrounds can easily be 
#                                   added in and scaled from first 
#                                   principles (Background numbers)
###########################################################################
import constants
import defaults
import physics
import error_utils

import math
import sys
import re

def convert_to_years(time_string):
    """ Convert any expression of time into years

    :param time_string: expression of time with units e.g "24.08 d" or
                        "33.27 m" or "16.5 ns"
    :type time_string: str
    """
    time_list = time_string.split()
    time = float(time_list[0])
    unit = time_list[1]
    years_per_year = 1.0
    days_per_year = 365.25
    hours_per_year = 8766
    minutes_per_year = 525960
    seconds_per_year = 3.15576e7
    milliseconds_per_year = 3.15576e10
    microseconds_per_year = 3.15576e13
    nanoseconds_per_year = 3.15576e16
    if (unit == "y"):
        time_in_years = time / years_per_year
    elif (unit == "d"):
        time_in_years = time / days_per_year
    elif (unit == "h"):
        time_in_years = time / hours_per_year
    elif (unit == "m"):
        time_in_years = time / minutes_per_year
    elif (unit == "s"):
        time_in_years = time / seconds_per_year
    elif (unit == "ms"):
        time_in_years = time / milliseconds_per_year
    elif (unit == "us"):
        time_in_years = time / microseconds_per_year
    elif (unit == "ns" ):
        time_in_years = time / nanoseconds_per_year
    else:
        print "isotope.convert_to_years: error - unit", unit, "not recognised"
        sys.exit(1)
    return time_in_years

class Isotope(object):
    """ Base isotope class. Obtains and stores the key parameters that are
    relevant for any radioactive isotope. The base class is intended to 
    have one instance per isotope and per decay mode, i.e. if an isotope
    can decay through both alpha and beta modes, a separate Isotope
    instance should be created for each.
    """
    def __init__(self, isotope_name):
        """ Initialises the isotope class with the name of the isotope

        :param isotope_name: name of isotope
        :type isotope_name: str
        """
        self._name = isotope_name
        self._half_life = None
        self._branching_fraction = 1.0
        self._mode = None
        self._number_nuclei = None
        self._counts = None
    def get_name(self):
        """
        :returns: self._name
        :rtype: str
        """
        return self._name
    def set_number_nuclei(self, isotope_mass):
        """ Calculate number of nuclei for a given target mass.
        
        :param isotope_mass: mass of isotope (kg)
        :type isotope_mass: float
        """
        error_utils.check_exists_in_dict("constants.atomic_weights",
                                         constants.atomic_weights,
                                         self._name,
                                         "Isotope.set_number_nuclei")
        atomic_weight = constants.atomic_weights.get(self._name)
        isotope_mass *= 1e3 # convert mass to grams
        self._number_nuclei = (isotope_mass*constants.n_avagadro)/atomic_weight
    def get_number_nuclei(self):
        """ Get number of nuclei for isotope.

        :returns: number of nuclei
        :rtype: float
        """
        if (self._number_nuclei == None):
            print "Isotope.get_number_nuclei: error number_nuclei for",\
                self._name, "not set"
            sys.exit(1)
        else:
            return self._number_nuclei
    def set_mode(self, mode):
        """ Set half life and corresponding mode for decay
        
        :param mode: decay mode e.g. alpha or beta
        :type mode: str
        """
        self._mode = mode
    def get_mode(self):
        """ Get the half life for this Isotope instance
        
        :returns: self._mode
        :rtype: float
        """
        if (self._mode == None):
            print "Isotope.get_mode: error mode for", self._name,\
                "not set"
            sys.exit(1)
        else:
            return self._mode
    def set_half_life(self, half_life="default", mode="default"):
        """ Set half life and corresponding mode for decay
        
        :param half_life: half life for decay (in years) (default: get from
                          constants). Can also supply a half life and 
                          corresponding branching fraction in a tuple
                          --> (half_life, branching fraction)
        :type half_life: float, string (e.g "24.10 d") or tuple (half_life,
                         branching_fraction)
        :param mode: decay mode e.g. alpha or beta (default: already set)
        :type mode: str
        """
        if (mode != "default"):
            self.set_mode(mode)
        self.get_mode()
        if (half_life == "default"):
            error_utils.check_exists_in_dict("constants.half_lives",
                                             constants.half_lives,
                                             self._name,
                                             "Isotope.set_half_life")
            error_utils.check_exists_in_dict\
                ("constants.half_lives.get("+self._name+")",
                 constants.half_lives.get(self._name),
                 get_mode(),
                 "Isotope.set_half_life")
            half_life = constants.half_lives.get(self._name).get(get_mode())
        if isinstance(half_life, tuple):
            half_life, self._branching_fraction = half_life
        if isinstance(half_life, str):
            assert(re.search\
                       (r"^[0-9]+.?[0-9]*[eE]?\+?\-?[0-9]*\s[dhmnsuy]{1,2}$", \
                            half_life) != None), \
                            "SNOPlusInternal.set_half_life: error - string " \
                            + half_life + " not recognised"
            half_life = convert_to_years(half_life)
        self._half_life = half_life
    def get_half_life(self):
        """ Get the half life for this Isotope instance
        
        :returns: self._half_life
        :rtype: float
        """
        if (self._half_life == None):
            print "Isotope.get_half_life: error half life for", self._name,\
                "not set"
            sys.exit(1)
        else:
            return self._half_life
    def set_counts(self, livetime="default",
                   half_life="default",
                   mode="default",
                   isotope_mass="default"):
        """ Sets number of counts (per year - default). Default behaviour is 
        to have half life and isotope mass (for number of nuclei) already set,
        or set from constants dicts but can be set from here as well.

        :param livetime: experiment livetime (default: livetime_general 
                         (defaults) --> counts per year
        :type livetime: float
        :param half_life: half life for decay (in years) (default: already
                          set). Can also supply a half life and 
                          corresponding branching fraction in a tuple
                          --> (half_life, branching fraction)
        :type half_life: float (or tuple)
        :param mode: decay mode e.g. alpha or beta (default: already set)
        :type mode: str
        :param isotope_mass: mass of isotope (kg) (default: already set)
        :type isotope_mass: float
        """
        if (livetime == "default"):
            livetime = defaults.livetime_general
        if (half_life != "default"):
            assert (mode != "default"), "Isotope.set_counts:"\
                "error half_life must have corresponding mode"
            set_half_life(half_life, mode)
        if (isotope_mass != "default"):
            self.set_number_nuclei(isotope_mass)
        counts = (math.log(2)*self.get_number_nuclei())/self.get_half_life() 
        counts *= self._branching_fraction
        counts *= livetime
        self._counts = counts
    def get_counts(self):
        """ Gets number of counts

        :returns: self._counts (default: counts per year)
        :rtype: float
        """
        if (self._counts == None):
            print "Isotope.get_counts: error - counts for", self._name,\
                "not set"
            sys.exit(1)
        else:
            return self._counts
class SNOPlusInternal(Isotope):
    """ Special case of Isotope class for any decays originating from 
    within the AV. Includes methods calculate mass based on scintillator
    cocktail. If you require a fiducial volume (FV) this is set in the
    set_scintillator_masses method.
    """
    def __init__(self, isotope_name):
        """ Initialises the SNOPlusInternal class with name of the isotope

        :param isotope_name: name of isotope
        :type isotope_name: str
        """
        super(SNOPlusInternal, self).__init__(isotope_name)
        self._scintillator_masses = None
        self._isotope_masses = {}
    def set_scintillator_masses(self, scintillator_cocktail="default",
                                apply_fv_cut=False, fv_radius="default"):
        """ Set masses of scintillator cocktail components. A fiducial
        volume (FV) cut at the specified FV radius is applied by default
        
        :param scintillator_cocktail: fractional components of scintillator
                                      e.g. { "LAB" : 0.927, "Te_{nat}" : 
                                      0.003, "H_{2}O" : 0.02, "PRS" : 0.05 }
                                      (default: get from defaults)
        :type scintillator_cocktail: dict
        :param apply_fv_cut: should FV cut be applied (default: True, yes)
        :type apply_fv_cut: bool
        :param fv_radius: radius of FV (default: get from defaults)
        """
        if (scintillator_cocktail == "default"):
            scintillator_cocktail = defaults.scintillator.get("cocktail")
        self._scintillator_cocktail = scintillator_cocktail
        scintillator_mass = defaults.scintillator.get("mass") # kg
        if apply_fv_cut:
            if (fv_radius == "default"):
                fv_radius = defaults.analysis.get("fv_radius")
            fiducial_volume = (4.0/3.0) * math.pi * math.pow(fv_radius, 3)
            av_volume = constants.av_volume
            scintillator_mass *= fiducial_volume / av_volume
        scintillator_mass *= 1.0e3 # convert to grams
        fractions_sum = 0.0
        for component, fraction in scintillator_cocktail.iteritems():
            fractions_sum += fraction
        assert (fractions_sum == 1.0), "SNOPlusInternal.__init__:"\
            "error - sum of components fractions != 1.0"
        scintillator_masses = {}
        masses_sum = 0.0
        for component, fraction in scintillator_cocktail.iteritems():
            component_mass = scintillator_mass * fraction
            masses_sum += component_mass
            scintillator_masses[component] = component_mass
        assert (masses_sum == scintillator_mass), "SNOPlusInternal.__init__: "\
            "error - sum of component masses != mass of scintillator"
        self._scintillator_masses = scintillator_masses
    def get_cocktail_component(self, component):
        """ Get the cocktail fraction of a scintillator component

        :param component: scintillator component for which to get the fraction
        :type component: str
        :returns: fraction of component
        :rtype: float
        """
        try:
            assert (self._scintillator_cocktail != None),\
                "dict scintillator_cocktail = None,\n"\
                " --> call SNOPlusInternal.set_scintillator_masses first"
        except AssertionError as detail:
            print "SNOPlusInternal.get_cocktail_component: error -", detail
            raise
        error_utils.check_exists_in_dict\
            ("SNOPlusInternal._scintillator_cocktail",
             self._scintillator_cocktail,
             component,
             "SNOPlusInternal.get_cocktail_component")
        return self._scintillator_cocktail.get(component)
    def add_isotope_mass(self, isotope_gram_per_gram, component):
        """ Add to total isotope mass by supplying a gIsotope/gComponent 
        (gram per gram) value and corresponding scintillator component
        
        :param isotope_gram_per_gram: gram per gram value
        :type isotope_gram_per_gram: double
        :param component: scintillator component to which gram per gram
                          value corresponds
        :type component: str
        """
        try:
            assert (self._scintillator_masses != None),\
                "dict scintillator_masses = None,\n"\
                " --> call SNOPlusInternal.set_scintillator_masses first"
        except AssertionError as detail:
            print "SNOPlusInternal.add_isotope_mass: error -", detail
            raise
        error_utils.check_exists_in_dict\
            ("SNOPlusInternal._scintillator_masses",
             self._scintillator_masses,
             component,
             "SNOPlusInternal.add_isotope_mass")
        self._isotope_masses[component] = \
            self._scintillator_masses.get(component) * isotope_gram_per_gram
    def get_mass_of_component(self, component):
        """ Get the mass of a scintillator component

        :param component: scintillator component for which to get the mass
        :type component: str
        :returns: mass of component in kg
        :rtype: float
        """
        try:
            assert (self._scintillator_masses != None),\
                "dict scintillator_masses = None,\n"\
                " --> call SNOPlusInternal.set_scintillator_masses first"
        except AssertionError as detail:
            print "SNOPlusInternal.get_mass_of_component: error -", detail
            raise
        error_utils.check_exists_in_dict\
            ("SNOPlusInternal._scintillator_masses",
             self._scintillator_masses,
             component,
             "SNOPlusInternal.get_mass_of_component")
        return self._scintillator_masses.get(component) / 1.0e3 # convert to kg
    def get_mass_in_component(self, component):
        """ Get the mass of isotope in a scintillator component
        
        :param component: scintillator component for which to get mass of 
                          isotope
        :type component: str
        :returns: mass of isotope in scintillator component
        :rtype: float
        """
        try:
            assert (self._isotope_masses != {}),\
                "dict isotope_masses = None,\n"\
                " --> call SNOPlusInternal.add_isotope_masses first"
        except AssertionError as detail:
            print "SNOPlusInternal.get_mass_in_component: error -", detail
            raise
        error_utils.check_exists_in_dict\
            ("Internal._isotope_masses",
             self._isotope_masses,
             component,
             "SNOPlusInternal.get_mass_in_component")
        mass_in_grams = self._isotope_masses.get(component)
        return mass_in_grams / 1.0e3 # convert to kg
    def get_mass(self):
        """ Get total mass of isotope in scintillator
        
        :returns: total mass
        :rtype: float
        """
        try:
            assert (self._isotope_masses != {}),\
                "dict isotope_masses = None,\n"\
                " --> call SNOPlusInternal.add_isotope_masses first"
        except AssertionError as detail:
            print "SNOPlusInternal.get_mass: error -", detail
            raise
        total_mass_in_grams = 0.0
        for mass in self._isotope_masses.values():
            total_mass_in_grams += mass
        return total_mass_in_grams / 1.0e3 # convert to kg
    def set_half_life(self, half_life="default", mode="default"):
        """ Set half life and corresponding mode for decay
        
        :param half_life: half life for decay (in years) (default: get from
                          constants). Can also supply a half life and 
                          corresponding branching fraction in a tuple
                          --> (half_life, branching fraction)
        :type half_life: float, string (e.g "24.10 d") or tuple (half_life,
                         branching_fraction)
        :param mode: decay mode e.g. alpha or beta (default: already set)
        :type mode: str
        """
        if (mode != "default"):
            self.set_mode(mode)
        self.get_mode()
        if (half_life == "default"):
            error_utils.check_exists_in_dict("constants.half_lives.internal",
                                             constants.half_lives.get("internal"),
                                             self._name,
                                             "SNOPlusInternal.set_half_life")
            error_utils.check_exists_in_dict\
                ("constants.half_lives.internal.get("+self._name+")",
                 constants.half_lives.get("internal").get(self._name),
                 self.get_mode(),
                 "SNOPlusInternal.set_half_life")
            half_life = constants.half_lives.get("internal").get(self._name).\
                get(self.get_mode())
        if isinstance(half_life, tuple):
            half_life, self._branching_fraction = half_life
        if isinstance(half_life, str):
            assert(re.search\
                       (r"^[0-9]+.?[0-9]*[eE]?\+?\-?[0-9]*\s[dhmnsuy]{1,2}$", \
                            half_life) != None), \
                            "SNOPlusInternal.set_half_life: error - string " \
                            + half_life + " not recognised"
            half_life = convert_to_years(half_life)
        self._half_life = half_life
    def set_counts_from_expected_counts(self, livetime="default",
                                      expected_counts="default"):
        """ Sets number of counts (per year - default) based on expected 
        yearly counts for a detector filled with pure LAB (0.78 kTonnes LAB)
        
        :param livetime: experiment livetime (default: livetime_general 
                         (defaults) --> counts per year
        :type livetime: float
        :param expected_counts: expected counts for decay (in 
                                counts/yr/kTonneLAB) (default: get from
                                constants)
        :type expected_counts: float
        """
        if (livetime == "default"):
            livetime = defaults.livetime_general
        if (expected_counts == "default"):
            error_utils.check_exists_in_dict("constants.expected_counts",
                                             constants.expected_counts,
                                             self._name,
                                             "SNOPlusInternal.set_counts_from_expected_counts")
            expected_counts = constants.expected_counts.get(self._name)\
                # counts/yr/kTonne pure LAB
        try:
            assert (self._isotope_masses != None),\
                "dict isotope_masses = None,\n"\
                " --> call SNOPlusInternal.set_scintillator_masses first"
        except AssertionError as detail:
            print "SNOPlusInternal.set_counts_from_expected_counts: error -", detail
            raise
        scintillator_counts = {}
        summed_counts = 0.0
        lab_mass = self.get_mass_in_component("LAB") / \
            self.get_cocktail_component("LAB") # mass of isotope if detector
                                               # was pure LAB
        print "lab_mass:", lab_mass
        for component, mass in self._isotope_masses.iteritems():
            # First work out equivalent mass of LAB
            error_utils.check_exists_in_dict("SNOPlusInternal.isotope_masses",
                                             self._isotope_masses,
                                             "LAB",
                                             "SNOPlusInternal.set_counts_from_expected_counts")
            print "mass:", mass
            scintillator_counts[component] = (mass/1e3) / lab_mass
            scintillator_counts[component] *= 0.78 # kTonnes scintillator
            scintillator_counts[component] *= expected_counts
            scintillator_counts[component] *= livetime
            summed_counts += scintillator_counts[component]
        self._scintillator_counts = scintillator_counts
        self._counts = summed_counts
    def get_scintillator_counts(self, component):
        """ Gets the counts contribution from the specified scintillator
        component
        
        :param component: scintillator component for which to get counts
        :type component: str
        :returns: number of counts
        :rtype: float
        """
        try:
            assert (self._scintillator_counts != None),\
                "dict scintillator_counts = None,\n"\
                " --> call SNOPlusInternal.set_counts_from_expected_counts first"
        except AssertionError as detail:
            print "SNOPlusInternal.get_scintillator_counts: error -", detail
            raise
        error_utils.check_exists_in_dict\
            ("SNOPlusInternal._scitillator_counts",
             self._scintillator_counts,
             component,
             "SNOPlusInternal.get_scintillator_counts")
        return self._scintillator_counts.get(component)
class DoubleBetaIsotope(SNOPlusInternal):
    """ Special case of SNOPlusInternal class for scintillator loaded with
    double beta isotope
    """
    def __init__(self, isotope_name):
        """ Initialises the DoubleBetaIsotope class with name of the 
        isotope

        :param isotope_name: name of isotope
        :type isotope_name: str
        """
        super(DoubleBetaIsotope, self).__init__(isotope_name)
        self._two_nu = physics.DoubleBeta(isotope_name)
    def set_mode(self, mode):
        """ Set mode for decay (as defined in decay0)
        
        :param mode: decay0 decay mode e.g. (0nuBetaBeta = 1, 
                     2nuBetaBeta = 4) 
        :type mode: int
        """
        assert isinstance(mode, int), \
            "DoubleBetaIsotope.set_mode: error - mode is not and integer\n"\
            " --> DoubleBetaIsotope class uses integer modes defined in decay0"
        self._mode = mode
    def get_mass(self):
        """ *** SNO+ (Te130) ONLY METHOD *** 
        Get total mass of isotope loaded into scintillator
        
        :returns: total mass
        :rtype: float
        """
        try:
            assert (self.get_name() == "Te130"),\
                "SNO+ Te130 only method!"
            assert (self._scintillator_masses != None),\
                "dict scintillator_masses = None,\n"\
                " --> call DoubleBetaIsotope.set_scintillator_masses first\n\n"
        except AssertionError as detail:
            print "DoubleBetaIsotope.get_mass: error -", detail
            raise
        error_utils.check_exists_in_dict\
            ("DoubleBetaIsotope._scintillator_masses",
             self._scintillator_masses,
             "Te_{nat}",
             "DoubleBetaIsotope.get_mass")
        natural_te_mass = self._scintillator_masses.get("Te_{nat}")
        return self.get_mass_of_isotope(natural_te_mass)
    def get_mass_of_isotope(self, natural_isotope_mass,
                            mass_fraction="default"):
        """ Get mass of specific isotope based on mass of naturally
        occurring element and mass_fraction
        
        :param natural_isotope_mass: mass of naturally occurring element
        :type natural_isotope_mass: float
        :param mass_fraction: fraction of naturally occurring mass that
                              corresponds to isotope
        :returns: mass of isotope (default: get from constants)
        :rtype: float
        """
        if (mass_fraction == "default"):
            error_utils.check_exists_in_dict\
                ("constants.mass_fractions",
                 constants.mass_fractions,
                 self._name,
                 "DoubleBetaIsotope.get_mass_of_isotope")
            mass_fraction = constants.mass_fractions.get(self._name)
        return natural_isotope_mass * mass_fraction
    def get_two_nu(self):
        """ Get two nu DoubleBeta object. Useful method for setting the
        half life for two neutrino double beta decay (mode = 4)

        :returns: self._two_nu
        :rtype: physics.DoubleBeta object
        """
        return self._two_nu
    def set_half_life(self, half_life="default", mode="default"):
        """ Overrides base class version to provide bounds checking on half 
        life. Note mode and half life have swapped order as there is more 
        choice of modes for double beta decays, so it is more likely that
        you would want to get the default half life for a specific mode.

        Tip: If you want to easily set the half life for 2 nu (mode=4), try
            
            DoubleBetaIsotope.set_half_life\\n
                (4, DoubleBetaIsotope.get_two_nu().get_half_life())

        :param half_life: half life for decay (in years) (default: get from
                          constants).
        :type half_life: float or string (e.g "24.10 d")
        :param mode: decay0 mode e.g. 0nu=1, 2nu=4 (default: already set)
        :type mode: str       
        """
        if (mode != "default"):
            self.set_mode(mode)
        self.get_mode()
        if (half_life == "default"):
            error_utils.check_exists_in_dict\
                ("constants.zero_nu_lifetimes",
                 constants.zero_nu_lifetimes,
                 self._name,
                 "DoubleBetaIsotope.set_half_life")
            half_life = constants.zero_nu_lifetimes.get(self._name).\
                get(self.get_mode())
        if isinstance(half_life, str):
            assert(re.search\
                       (r"^[0-9]+.?[0-9]*[eE]?\+?\-?[0-9]*\s[dhmnsuy]{1,2}$", \
                            half_life) != None), \
                            "SNOPlusInternal.set_half_life: error - string " \
                            + half_life + " not recognised"
            half_life = convert_to_years(half_life)   
        t_half_min = self._two_nu.get_t_half_min()
        t_half_max = self._two_nu.get_t_half_max()
        if (half_life < t_half_min):
            half_life = t_half_min
            print "DoubleBetaIsotope.set_half_life: WARNING setting t_half = "\
                "t_half_min"
        if (half_life > t_half_max):
            half_life = t_half_max
            print "DoubleBetaIsotope.set_half_life: WARNING setting t_half = "\
                "t_half_max"
        try:
            assert (t_half_min <= half_life <= t_half_max),\
                "half life does not fall within the accepted range"
        except AssertionError as detail:
            print "DoubleBetaIsotope.set_half_life: ERROR", detail
            print " --> cannot calculate decays"
            raise
        self._half_life = half_life
class ZeroNu(DoubleBetaIsotope):
    """ Special case of DoubleBetaIsotope class for neutrinoless double
    beta decays
    """
    def __init__(self, isotope_name):
        """ Initialises the DoubleBetaIsotope class with name of the 
        isotope. Sets mode as 1 automatically, this class is only valid
        for mode=1.

        :param isotope_name: name of isotope
        :type isotope_name: str
        """
        super(ZeroNu, self).__init__(isotope_name)
        self._zero_nu = physics.ZeroNuConverter(self._name)
        self.set_mode(1)
        self._effective_mass = None
    def get_zero_nu(self):
        """ Get zero nu DoubleBeta object. Useful method for converting to
        effective mass

        :returns: self._zero_nu
        :rtype: physics.ZeroNuConverter object
        """
        return self._zero_nu
    def set_effective_mass(self, effective_mass):
        """ Sets the effective mass for a neutrinoless double beta decay
        (mode = 1 ONLY) double beta isotpe
        
        :param effective_mass: effective mass for decay in eV
        :type effective_mass: float
        """
        assert (self.get_mode() == 1), \
            "ZeroNu.set_effective_mass: error - invalid mode\n"\
            " --> ZeroNu class is only valid for mode=1!"
        effective_mass_min = self._zero_nu.get_mass_min()
        effective_mass_max = self._zero_nu.get_mass_max()
        if (effective_mass < effective_mass_min):
            effective_mass = effective_mass_min
            print "ZeroNu.set_effective_mass: WARNING setting effective_mass"\
                " = effective_mass_min"
        if (effective_mass > effective_mass_max):
            effective_mass = effective_mass_max
            print "ZeroNu.set_effective_mass: WARNING setting effective_mass"\
                " = effective_mass_max"
        try:
            assert (effective_mass_min <= effective_mass <= effective_mass_max),\
                "half life does not fall within the accepted range"
        except AssertionError as detail:
            print "ZeroNu.set_effective_mass: ERROR", detail
            print " --> cannont calculate decays"
            raise
        self._effective_mass = effective_mass
    def get_effective_mass(self):
        """ Get the effective mass for this ZeroNu instance
        
        :returns: self._effective_mass
        :rtype: float
        """
        if (self._effective_mass == None):
            print "ZeroNu.get_effective_mass: error - effective mass for", self._name,\
                "not set"
            sys.exit(1)
        else:
            return self._effective_mass
    def set_counts_by_mass(self, livetime="default",
                           effective_mass="default",
                           isotope_mass="default"):
        """ Sets number of counts (per year - default). Default behaviour is 
        to have effective mass and isotope mass (for number of nuclei) 
        already set, or set from constants dicts but can be set from here 
        as well.

        :param livetime: experiment livetime (default: livetime_general 
                         (defaults) --> counts per year
        :type livetime: float
        :param effective_mass: effective mass for decay (in eV) (default:
                               already set)
        :type effective_mass: float
        :param isotope_mass: mass of isotope (kg) (default: already set)
        :type isotope_mass: float
        """
        assert (self.get_mode() == 1), \
            "ZeroNu.set_counts_by_mass: error - invalid mode\n"\
            " --> ZeroNu class is only valid for mode=1!"
        if (livetime == "default"):
            livetime = defaults.livetime_general
        if (effective_mass != "default"):
            set_effective_mass(effective_mass)
        if (isotope_mass != "default"):
            set_number_nuclei(isotope_mass)
        counts = math.log(2)*math.pow\
            (self.get_effective_mass()/self._zero_nu.get_conversion_factor(), 2)\
            * self.get_number_nuclei()
        counts *= livetime
        self._counts = counts

############################################################################
if (__name__ == "__main__"):
    U238 = SNOPlusInternal("U238")
    U238.set_mode("alpha")
    U238.set_half_life()
    U238.set_scintillator_masses()
    U238.add_isotope_mass(1.0e-15, "Te_{nat}")
    U238.add_isotope_mass(1.60e-17, "LAB")
    U238.add_isotope_mass(3.5e-14, "H_{2}O")
    U238.add_isotope_mass(3.5e-14, "PRS")
    U238.set_number_nuclei(U238.get_mass())
    U238.set_counts()
    print U238.get_name(), U238.get_counts(), "events/yr"
    U238.set_counts_from_expected_counts()
    print U238.get_name() + " " + str(U238.get_scintillator_counts("LAB")) + \
        " + " + str(U238.get_scintillator_counts("Te_{nat}")) + \
        " + " + str(U238.get_scintillator_counts("H_{2}O")) + \
        " + " + str(U238.get_scintillator_counts("PRS")) + \
        " = " + str(U238.get_counts())

