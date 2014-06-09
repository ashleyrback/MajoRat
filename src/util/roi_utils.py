#!/usr/bin/env python
#
# roi_utils.py
#
# Class for handling regions of interest
#
# Author A R Back 
#
# 04/06/2014 <ab571@sussex.ac.uk> : First revision
#
###########################################################################
import constants

class RoIUtils(object):
    """ Class for handling and creating regions of interest

    """
    def __init__(self, roi_name="default"):
        """ Initialise class with the name of a ROI
         - reconstructed_energy
         - gaussian_smeared_mc_energy
         - nhit_per_180

        :param roi_name: name of ROI to initialise
        :type roi_name: str
        """
        self._name = roi_name
        if (self._name == "default"):
            self._e_lo = None
            self._e_hi = None
            self._mu = None
            self._sigma = None
        else:
            try:
                self._e_lo = constants.roi.get(self._name).get("e_lo")
                self._e_hi = constants.roi.get(self._name).get("e_hi")
                self._mu = float((self._e_lo+self._e_hi)/2.0)
                self._sigma = self._e_hi - self._mu
            except AttributeError:
                print "RoIUtils.__init__: ERROR, could not find ROI", roi_name
    def get_lower_limit(self):
        """
        :returns: self._e_lo
        :rtype: float
        """
        return self._e_lo
    def get_upper_limit(self):
        """
        :returns self._e_hi
        :rtype: float
        """
        return self._e_hi
    def set_lower_by_sigma(self, offset):
        """
        :param offset: sets lower limit as self._mu + (offset*self._sigma)
        :type offset: float
        """
        self._e_lo = self._mu + offset*self._sigma
    def set_upper_by_sigma(self, offset):
        """
        :param offset: sets lower limit as self._mu + (offset*self._sigma)
        :type offset: float
        """
        self._e_hi = self._mu + offset*self._sigma
        
