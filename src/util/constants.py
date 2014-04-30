#!/usr/bin/env python
#
# constants.py
#
# Main constants associated with double beta decay
#
# Author A R Back 
#
# 12/02/2014 <ab571@sussex.ac.uk> : First revision
# 29/04/2014 <ab571@sussex.ac.uk> : Major re-facoring for log-likelihood
#                                   calculations moved some isotopic info 
#                                   here in isotope dicts. Moved SNO+ 
#                                   specific info to SNOPlus(isotope) class
###############################################################################

### Physical constants ###
n_avagadro = 6.02214129e23 # mol^-1 # PDG 2012

### Isotope dictionaries ###
mass_fractions = { "Te130" : 0.34696, } # Te Verification Report
atomic_weights = { "Te130" : 129.906229, } # gmol^-1 # Te Verification Report
isotope_mass = { "Xe136" : 125.0 } # +/- 7 kg # Gando et al. 2012
half_life = { "Xe136" : { 0 : 5.0e25, # KLZ limits Gando et al. 2012
                          1 : 3.0e24,
                          2 : 1.0e24,
                          3 : 5.0e23,
                          5 : 1.0e21,
                          7 : 1.0e22 } }

### SNO+ detector ###
snoplus = { "fv_radius" : 4.0 } # Assumong 4.0 m FV

