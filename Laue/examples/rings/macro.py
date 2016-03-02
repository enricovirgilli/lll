#!/usr/bin/env python
from LLL import Lenses
# How to use the Generic lens object
# Lens definition
lens=Lenses.Rings(name = "rings", # prefix for file names
                  profile = "spherical", # lens profile
                  thickness_profile = "fixed", # crystal thickness policy
                  thickness_threshold = (0, 10), # crystal thickness limits
                  thickness_factor = 1., # thickness reduction for saving mass
                  datadir="data", # directory where files are stored
                  Focal = 2000., # focal length (cm)
                  Emin=200., # minimum of the nominal energy working band (keV)
                  Emax=205., # maximum of the nominal energy working band (keV)
                  dim=[3., 1., 0.4], # crystal tile (radial, tangential, thickness in cm)
                  tan_framewidth = 0.10, # space around each tile (cm)
                  rad_framewidth = 0.1,
                  fwhm=1., # mosaic crystal FWHM (arcmin)
                  hkl=[1,1,1], # Miller indices
                  Z=32, # atomic number
                  microthick=0., # microblock thickness (A)
                  pixels=(50,50), # detector pixels (x,y)
                  pitches=(0.2,0.2), # detector pitch (x,y in cm)
                  error=0., # orientation error
                  )
# Running the program
lens.macro(INFO=True, # generic lens info
           XTALINFO=True, # crystal position and orientation
           ERANGE=(200,205), # energy range to be investigated
           AREA=True, # effective area
           MC=True, # Monte Carlo
           SENS=True, # sensitivity (needs MC)
           S=True, # equivalent surface (needs MC)
           G=True, # focusing factor (needs MC)
           PHOTONS=1e4, # number of photons (MC)
           OFFSET=0.0, # source offset
           FRACTION=0.5, # considered photon fraction
           Tobs=1.e6, # observation time
          )
