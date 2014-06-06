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
# 06/05/2014 <ab571@sussex.ac.uk> : Added matrix elements and phase space 
#                                   factors
# 04/06/2014 <ab571@sussex.ac.uk> : Added ROI info
# 05/06/2014 <ab571@sussex.ac.uk> : Moved snoplus experiment constants to 
#                                   defaults.
###########################################################################

## Physical constants ##
n_avagadro = 6.02214129e23 # mol^-1 # PDG 2012
electron_mass = 0.510998928 # MeV PDG 2012

## SNO+ specific ##
# Constants defined here are taken from the Te Verification Report,
# unless stated otherwise stated
density_lab = 0.862e3 # kgm^-3
av_volume = 903.3 # m^3

## ROI ##
# Taken from "External Background Simulation Studies for NewNd" 
# SNO+-doc-1720-v7
roi = {
    "reconstructed_energy" : {
        "e_lo" : 2.40, # MeV
        "e_hi" : 2.64
        },
    "gaussian_smeared_mc_energy" : {
        "e_lo" : 2.450,
        "e_hi" : 2.635
        },
    "nhit_over_180" : {
        "e_lo" : 2.53,
        "e_hi" : 2.89
        },
    }
    
## Isotope dictionaries ##

### Mass, weight and abundance ###
mass_fractions = { "Te130" : 0.34696, } # Te Verification Report
atomic_weights = { "Te130" : 129.906229, } # gmol^-1 # Te Verification 
                                           # Report (TVR)
isotope_mass = { "Xe136" : 125.0 } # +/- 7 kg # Gando et al. 2012

### Matrix elements, phase space and lifetime ###
""" g_A Phys. Rev. C 87, 014315 (2013), Table IV, as used in TVR """
coupling_constant = 1.269
""" |M| - nuclear matrix elements, assumes:
 * decay to GS (0_1^+)
 * IBM-2
source: Phys. Rev. C 87, 014315 (2013), Table IV (0nu), Table XV (2nu)
"""
matrix_element = {
    "Te130" : {
        "0nu" : 4.03,
        "2nu" : 3.31,
        },
    "Xe136" : {
        "0nu" : 3.33, 
        "2nu" : 2.76,
        },
    "Nd150" : {
        "0nu" : 2.32,
        "2nu" : 1.54,
        },
    }
""" G^(0) - phase space (yr^-1)
 * screened exact finite-sized Coulomb wave functions
source: Phys. Rev. C 85, 034316 (2012), (0nu) Table III, (2nu) Table I
"""
phase_space = { 
    "Te130" : {
        "0nu" : 14.22e-15,
        "2nu" : 1529e-21,
        },
    "Xe136" : {
        "0nu" : 14.58e-15,
        "2nu" : 1433e-21,
        },
    "Nd150" : {
        "0nu" : 63.03e-15,
        "2nu" : 36430e-21,
        },
    }
half_life = { 
    "Xe136" : {
        0 : 6.2e24, # KLZ limits Gando et al. 2012
        1 : 2.6e24,
        2 : 1.0e24,
        3 : 4.5e23,
        5 : 2.3e21,
        7 : 1.1e22, 
        },
    }
