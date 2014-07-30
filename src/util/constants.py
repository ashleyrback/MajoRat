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

# Physical constants #
n_avagadro = 6.02214129e23 # mol^-1 # PDG 2012
electron_mass = 0.510998928 # MeV PDG 2012

# SNO+ specific #
# Constants defined here are taken from the Te Verification Report,
# unless stated otherwise stated
density_lab = 0.862e3 # kgm^-3
av_volume = 903.3 # m^3

# ROI #
# Taken from "External Background Simulation Studies for NewNd" 
# SNO+-doc-1720-v7
rois = {
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
    
# Isotopes #

## Relating to any Isotope ##

### Mass, weight and abundance ###
atomic_weights = { # from http://physics.nist.gov/cgi-bin/Compositions/stand_al
                   # one.pl unless otherwise stated
    "Te130" : 129.906229,  # Te Verification Report (TVR)
    "Xe136" : 135.907219,
    "U238" : 238.0507882,
    } # gmol^-1
isotope_masses = { "Xe136" : 125.0 } # +/- 7 kg # Gando et al. 2012

## Double beta isotopes ONLY ##
mass_fractions = { "Te130" : 0.34696, } # Te Verification Report

### Matrix elements, phase space and lifetime ###
""" g_A Phys. Rev. C 87, 014315 (2013), Table IV, as used in TVR """
coupling_constant = 1.269
""" |M| - nuclear matrix elements, assumes:
 * decay to GS (0_1^+)
 * IBM-2
source: Phys. Rev. C 87, 014315 (2013), Table IV (0nu), Table XV (2nu)
"""
matrix_elements = {
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
phase_spaces = { 
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
zero_nu_lifetimes = { # 
    "Xe136" : {
        0 : "6.2e24 y", # KLZ limits Gando et al. 2012
        1 : "2.6e24 y",
        2 : "1.0e24 y",
        3 : "4.5e23 y",
        5 : "2.3e21 y",
        7 : "1.1e22 y", 
        },
    }

## Backgrounds ONLY ##
half_lives = { # years
    "internal" : {
        ### Uranium chain ###
        "U238" : {
            "alpha" : "4.468E+9 y"
            },
        "Th234" : {
            "beta" : "24.10 d"
            },
        "Pa234m" : {
            "beta" : "1.159 m"
            },
        "U234" : {
            "alpha" : "2.455E+5 y"
            },
        "Th230" : {
            "alpha" : "7.538E+4 y"
            },
        "Ra226" : {
            "alpha" : "1600 y"
            },
        "Rn222" : {
            "alpha" : "3.8235 d"
            },
        "Po218" : {
            "alpha" : "3.098 m"
            },
        "Pb214" : {
            "beta" : "26.8 m"
            },
        "Bi214" : {
            "alpha" : ("19.9 m", 0.00021),
            "beta" : ("19.9 m" , 0.99979)
            },
        "Tl210" : {
            "beta" : ("1.30 m", 0.00021)
            },
        "Po214" : {
            "alpha" : ("164.3 us", 0.99979)
            },
        "Pb210" : {
            "beta" : "22.20 y"
            },
        "Bi210" : {
            "beta" : "5.012 d"
            },
        "Po210" : {
            "alpha" : "138.376 d"
            },
        #####################
        ### Thorium Chain ###
        "Th232" : {
            "alpha" : "14.0E+9 y"
            },
        "Ra228" : {
            "beta" : "5.75 y"
            },
        "Ac228" : {
            "beta" : "6.15 h"
            },
        "Th228" : {
            "alpha" : "1.912 y"
            },
        "Ra224" : {
            "alpha" : "3.66 d"
            },
        "Rn220" : {
            "alpha" : "55.6 s"
            },
        "Po216" : {
            "alpha" : "0.145 s"
            },
        "Pb212" : {
            "beta" : "10.64 h"
            },
        "Bi212" : {
            "alpha" : ("60.55 m" , 0.3594),
            "beta" : ("60.55 m" , 0.6406)
            },
        "Tl208" : {
            "beta" : ("3.053 m", 0.3594)
            },
        "Po212" : {
            "alpha" : ("0.299 us", 0.6404)
            },
        #####################
        "Ar39" : {
            "beta" : "269 y"
            },
        "C14" : {
            "beta" : "5700 y"
            },
        "Kr85" : {
            "beta" : "10.739 y"
            },
        "K40" : {
            "beta" : ("1.248E+9 y", 0.8928),
            "beta+" : ("1.248E+9 y", 0.1072)
            }
        },
    }

expected_counts = { # from SNO+-docdb-507v18, Table 2 (counts/year/kTonneLAB)
    ### Uranium chain ###
    "U238" : 6278,
    "Th234" : 1.9872e-10,
    "Pa234m" : 1.9872e-10,
    "U234" : 1.9872e-10,
    "Th230" : 1.9872e-10,
    "Ra226" : 1.9872e-10,
    "Rn222" : 1.9872e-10,
    "Po218" : 1.9872e-10,
    "Pb214" : 1.9872e-10,
    "Bi214" : 1.9872e-10,
    "Tl210" : 4.1795e-14,
    "Po214" : 1.9872e-10,
    "Pb210" : 1.7308e-9,
    "Bi210" : 1.7308e-9,
    "Po210" : 6.92731e-7,
    #####################
    ### Thorium Chain ###
    #"Th232" : ,
    #"Ra228" : ,
    #"Ac228" : ,
    #"Th228" : ,
    #"Ra224" : ,
    #"Rn220" : ,
    #"Po216" : ,
    #"Pb212" : ,
    #"Bi212" : ,
    #"Tl208" : ,
    #"Po212" : ,
    #####################
    #"Ar39" : ,
    #"C14" : ,
    #"Kr85" : ,
    #"K40" : 
    }
