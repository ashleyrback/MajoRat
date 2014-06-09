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

### SNO+ experiment defaults ###
snoplus = { 
    "fv_radius" : 3.5, # Assumed FV radius (m)
    "te_loading" : 0.003, # 0.3% loading
    "energy_resolution" : 180.0, # NHit per MeV
    "livetime" : 2.0 # Livetime used for spectral plots SNO+-doc-1975-v1
    }

## Spectra ##
spectrum = {
    "bin_width" : 0.02,
    "e_lo" : 0.0,
    "e_hi" : 3.5
    }

## Other defaults ##
livetime_general = 1.0 # years
target_mass=240.534021855 # Expected SNO+ ^{nat}Te mass
