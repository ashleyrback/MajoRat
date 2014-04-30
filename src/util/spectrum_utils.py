#!/usr/bin/env python
#
# spectrum_utils.py
#
# A collection of methods relating to spectra
#
# Author A R Back 
#
# 29/04/2014 <ab571@sussex.ac.uk> : First revision
###############################################################################
def get_spectral_index(decay0_mode):
    """ For a supplied decay0_mode, returns the spectral index """
    mode_to_index = {1:0, 4:5, 5:1, 6:2, 7:3, 8:7}
    spectral_index = mode_to_index.get(decay0_mode)
    return spectral_index
def get_label(spectral_index):
    """ Returns the correct (Root formatted) label, for a given spectral
    index
    """
    index_to_label = {0 : "0#nu#beta#beta",
                      1 : "0#nu#beta#beta#chi^{0} (n=1)",
                      2 : "0#nu#beta#beta#chi^{0} 'bulk' (n=2)",
                      3 : "0#nu#beta#beta#chi^{0}(#chi^{0}) (n=3)",
                      5 : "2#nu#beta#beta (n=5)",
                      7 : "0#nu#beta#beta#chi^{0}#chi^{0} (n=7)"}
    label = index_to_label.get(spectral_index)
    return label
