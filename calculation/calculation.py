#!/usr/bin/env python
#
# signal_esitmates.py
#
# Back-of-the-envelope-style calculations of expected signals
#
# Author A R Back - 12/02/2014 <ab571@sussex.ac.uk> : First revision
###############################################################################
import os
from spectrum_data import SpectrumData
from detector_parameters import mass_Te_in_FV
from detector_parameters import number_Te
from detector_parameters import n_decays

# Calculate number of nuclei
number_Te = number_Te(mass_Te_in_FV(0.003, 4.0))

half_life = {0 : 5.0e25,
             1 : 3.0e24,
             2 : 1.0e24,
             3 : 5.0e23,
             5 : 1.0e21,
             7 : 1.0e22}

class SignalEstimates(SpectrumData):
    """ Derived class, adapts spectrum data class for looking at estimated 
    signals and backgrounds.
    """
    def __init__(self):
        """ Initialise the container """
        super(SinalEstimates, self).__init__(filepath)
        if (self._mode is "1") or (self._mode is "4"):
            self._signal = False
        else:
            self._signal = True
        self._t_half = half_life.get(self._spectral_index)
        self._events
        self._fraction_of_events = self._get_fraction_events(2.0)
        

    def get_fraction_events(self, elo=0.0, ehi=3.0):
        """ Integrates the normalised energy spectrum between the limits E_lo
        and E_hi. Returns the fraction of events within these limits.
        """
        spectrum = self.GetHistogram()
        n_bins = spectrum.GetNbinsX()
        bin_width = spectrum.GetBinWidth()
        integrate_from = elo/bin_width
        integrate_to = ehi/bin_width
        spectrum.Scale(1/spectrum.Integral("width"))
        fraction = spectrum.Integral(integrate_from, integrate_to)
        return fraction

print __name__

if __name__ == "__main__":
    files = ["RAT4-5_1k_decay0_2beta_Te130_0_1_ELT3-5.root",
             "RAT4-5_1k_decay0_2beta_Te130_0_4_ELT3-5.root",
             "RAT4-5_1k_decay0_2beta_Te130_0_5_ELT3-5.root",
             "RAT4-5_1k_decay0_2beta_Te130_0_6_ELT3-5.root",
             "RAT4-5_1k_decay0_2beta_Te130_0_7_ELT3-5.root",
             "RAT4-5_1k_decay0_2beta_Te130_0_8_ELT3-5.root"]
    print files

