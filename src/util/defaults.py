#!/usr/bin/env python
#
# defaults.py
#
# Default values for parameters used by MajoRat
#
# Author A R Back 
#
# 05/06/2014 <ab571@sussex.ac.uk> : First revision
#
###########################################################################

# Spectra #
spectrum = {
    "bin_width" : 0.02,
    "e_lo" : 0.0,
    "e_hi" : 3.5
    }

# Analysis #
analysis = {
    "fv_radius" : 3.5, # Assumed FV radius (m)
    "production_label" : {
        "RAT4.5" : {
            "energy_resolution" : 200.0 # NHit per MeV
            },
        "RAT4.4" : {
            "energy_resolution" : 180.0 # NHit per MeV
            }
        }
    }

# Scaling #
## Scintillator ##
scintillator = {
    "cocktail" : {
        "PRS" : 0.05,
        "H_{2}O" : 0.02,
        "Te_{nat}" : 0.003,
        "LAB" : 0.927
        },
    "mass" : 780.00e3 # kg (780 Tonnes)
    }

## Spectral plots ##
spectral_plot = {
    "livetime" : 2.0, # Livetime used for spectral plots SNO+-doc-1975-v1
    "effective_mass" : 0.27 # eV
    }

# Log-likelihood analysis #
ll_analysis = {
    "delta_chi_squared" : 2.71, # corresponds to 90% CL PDG 2012
    "livetime" : 5.0, # Yr
    "outer_mass_lo" : 0.0, # eV
    "outer_mass_hi" : 0.27, # eV
    "outer_mass_step" : 0.01, # eV
    "inner_mass_lo" : 0.0, # eV
    "inner_mass_hi" : 0.27, # eV
    "inner_mass_step" : 0.01, # eV
    "outer_t_half_lo" : 0.0, # Yr
    "outer_t_half_hi" : 1.0e25, # Yr
    "outer_t_half_test" : 0.5e24, # Yr
    "inner_t_half_lo" : 0.0, # Yr
    "inner_t_half_hi" : 1.0e25, # Yr
    "inner_t_half_test" : 0.5e24 # Yr
    }

## Other defaults ##
livetime_general = 1.0 # years
target_mass=240.534021855 # Expected SNO+ ^{nat}Te mass
