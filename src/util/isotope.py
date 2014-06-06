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
# 09/05/2014 <ab571@sussex.ac.uk> : Added get_decays_from_mass method and 
#                                   modified get_decays_from_t_half method
# 21/05/2014 <ab571@sussex.ac.uk> : Valid range for half life and 
#                                   corresponding effective mass.
# 05/06/2014 <ab571@sussex.ac.uk> : Improved defaults
###########################################################################
import constants
import defaults
import physics

import math
import sys

class Isotope(object):
    """ Base isotope class, for dealing with numbers of isotopes/decays
    for any isotope
    """
    def __init__(self, isotope_name):
        """ Initialise class based on the name of the isotope """
        self._name = isotope_name
        self._zero_nu_converter = physics.ZeroNuConverter(self._name)
        self._number_nuclei = None
    def get_name(self):
        """
        :return: isotope name
        :rtype: str
        """
        return self._name
    def get_number_nuclei(self, isotope_mass="default"): 
        """ Number of nuclei for a given target mass.
        
        :param isotope_mass: mass of isotope/target (kg)
        :type isotope_mass: float
        :returns: number of nuclei
        :rtype: float
        """
        if (isotope_mass == "default"):
            isotope_mass = defaults.target_mass
        try:
            assert constants.atomic_weights.get(self._name) != None, \
                "Atomic weight " + self._name + " is not available"
        except AssertionError as detail:
            print "Isotope:get_number_nuclei error:", detail
            sys.exit(1)
        atomic_weight = constants.atomic_weights.get(self._name)
        isotope_mass *= 1e3 # convert mass to grams
        self._number_nuclei = (isotope_mass*constants.n_avagadro)/atomic_weight
        return self._number_nuclei
    def get_decays_from_t_half(self, t_half, 
                               livetime="default", 
                               number_nuclei="default"):
        """ Expected number of decays.

        :param t_half: half life in years (default = self._t_half)
        :type t_half: float
        :param livetime: livetime in years
        :type livetime: float
        :param number_nuclei: number of target nuclei in target mass
        :type number_nuclei: float
        :returns: expected number of decays
        :rtype: float
        """
        if (livetime == "default"):
            livetime = defaults.livetime_general
        if (number_nuclei == "default"):
            if (self._number_nuclei == None):
                self.get_number_nuclei()
            number_nuclei = self._number_nuclei
        t_half_min = self._zero_nu_converter.get_t_half_min()
        t_half_max = self._zero_nu_converter.get_t_half_max()
        if (t_half < t_half_min):
            t_half = t_half_min
            print "Isotope.get_number_of_decays: WARNING setting t_half = "\
                "t_half_min"
        if (t_half > t_half_max):
            t_half = t_half_max
            print "Isotope.get_number_of_decays: WARNING setting t_half = "\
                "t_half_max"
        try:
            assert (t_half_min <= t_half <= t_half_max),\
                "half life does not fall within the accepted range"
        except AssertionError as detail:
            print "Isotope.get_decays_from_t_half: ERROR", detail
            print " --> cannot calculate decays"
            sys.exit(1)
        decay_rate = math.log(2)/t_half
        number_decays = decay_rate*number_nuclei*livetime
        return number_decays
    def get_decays_from_mass(self, effective_mass,
                             livetime="default",
                             number_nuclei="default"):
        """ Expected number of decays.

        :param effective_mass: effective double beta mass in eV (no default)
        :type effective_mass: float
        :param livetime: livetime in years
        :type livetime: float
        :param number_nuclei: number of target nuclei in target mass
        :type number_nuclei: float
        :returns: expected number of decays
        :rtype: float
        """
        if (livetime == "default"):
            livetime = defaults.livetime_general
        if (number_nuclei == "default"):
            if (self._number_nuclei == None):
                self.get_number_nuclei()
            number_nuclei = self._number_nuclei
        effective_mass_min = self._zero_nu_converter.get_mass_min()
        effective_mass_max = self._zero_nu_converter.get_mass_max()
        if (effective_mass < effective_mass_min):
            effective_mass = effective_mass_min
            print "Isotope.get_number_of_decays: WARNING setting effective_mass"\
                " = effective_mass_min"
        if (effective_mass > effective_mass_max):
            effective_mass = effective_mass_max
            print "Isotope.get_number_of_decays: WARNING setting effective_mass"\
                " = effective_mass_max"
        try:
            assert (effective_mass_min <= effective_mass <= effective_mass_max),\
                "half life does not fall within the accepted range"
        except AssertionError as detail:
            print "Isotope.get_decays_from_mass: ERROR", detail
            print " --> cannont calculate decays"
            sys.exit(1)
        decay_rate = math.log(2)*math.pow\
            (effective_mass/self._zero_nu_converter.get_conversion_factor(), 2)
        number_decays = decay_rate*number_nuclei*livetime
        return number_decays

### SNO+ SPECIFIC #############################################################
class SNOPlusTe(Isotope):
    """ Special case for SNO+ Te-130, derived from isotope class. Provides
    additional useful methods for getting mass of Te inside AV, or a given
    fiducial volume.
    """
    def __init__(self, isotope_name="Te130"):
        """ Initialised with choice of isotope.

        :param isotope_name: name of isotope, Te130 default
        :type isotope_name: str
        """
        super(SNOPlusTe, self).__init__(isotope_name)
        self._av_volume = constants.av_volume
        self._density_lab = constants.density_lab
    def get_mass(self, loading="default"):
        """ Mass of Te-130 in AV, for given natural Te loading fraction.

        :param loading: loading fraction of Te-nat
        :type loading: float
        :returns: mass Te130, in given AV
        :rtype: float
        """
        if (loading == "default"):
            loading = defaults.snoplus.get("te_loading")
        mass_lab = self._density_lab*self._av_volume
        mass = mass_lab*constants.mass_fractions.get(self._name)*loading
        return mass
    def get_mass_in_fv(self, radius="default", loading="default"):
        """ Mass of Te-130 in a FV, for given natural Te loading fraction.

        :param radius: specify FV radius (m)
        :type radius: float
        :param loading: loading fraction of Te-nat
        :type loading: float
        :returns: mass Te130, in given FV
        :rtype: float
        """
        if (radius == "default"):
            radius = defaults.snoplus.get("fv_radius")
        if (loading == "default"):
            loading = defaults.snoplus.get("te_loading")
        volume = math.pow(radius, 3)
        volume *= math.pi
        volume *= (4.0/3.0)
        mass_lab = self._density_lab*volume
        mass = mass_lab*constants.mass_fractions.get(self._name)*loading
        return mass
    def get_number_nuclei(self, fv_cut=True, 
                          radius="default", 
                          loading="default"):
        """ Number of nuclei that could decay.

        :param fv_cut: cut on fiducical volume (FV)?
        :type fv_cut: bool
        :param radius: specify FV radius (m)
        :type radius: float
        :param loading: loading fraction of Te-nat
        :type loading: float
        :returns: number of nuclei
        :rtype: float
        """
        if (radius == "default"):
            radius = defaults.snoplus.get("fv_radius")
        if (loading == "default"):
            loading = defaults.snoplus.get("te_loading")
        if (fv_cut):
            mass = self.get_mass_in_fv(radius, loading)
        else:
            mass = self.get_mass(loading)
        self._number_nuclei = super(SNOPlusTe, self).get_number_nuclei(mass)
        return self._number_nuclei

###############################################################################
if __name__ == "__main__":
    te_converter = physics.ZeroNuConverter("Te130")
    spectrum_isotope = SNOPlusTe()

    t_half = 5.0e25 # years # KLZ limit from Gando et al.
    decays = spectrum_isotope.get_decays_from_t_half(t_half)
    print "Te130: t_half =", t_half, "decays =", decays
    effective_mass = te_converter.half_life_to_mass(t_half)
    decays = spectrum_isotope.get_decays_from_mass(effective_mass)
    print "Te130: effective mass =", effective_mass, "decays =", decays

    # Experiment with limits
    effective_mass_max = te_converter.half_life_to_mass(0.0) # when t_half=0
    t_half_min = te_converter.mass_to_half_life(effective_mass_max)
    t_half_max = te_converter.mass_to_half_life(0.0) # when mass=0
    effective_mass_min = te_converter.half_life_to_mass(t_half_max)

    max_decays = spectrum_isotope.get_decays_from_t_half(t_half_min)
    min_decays = spectrum_isotope.get_decays_from_t_half(t_half_max)
    print "Max decays =", max_decays, "min decays =", min_decays
    max_decays = spectrum_isotope.get_decays_from_mass(effective_mass_max)
    min_decays = spectrum_isotope.get_decays_from_mass(effective_mass_min)
    print "Max decays =", max_decays, "min decays =", min_decays

