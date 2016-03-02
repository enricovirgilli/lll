#!/usr/bin/env python
from LLL import Lenses
# How to use the Generic lens object
# Lens definition
lens=Lenses.Spiral(name = "spiral", # prefix for file names
                   profile = "spherical", # lens profile
                   thickness_profile = "fixed", # crystal thickness policy
                   thickness_threshold = (0, 10), # crystal thickness limits
                   thickness_factor = 1., # thickness reduction for saving mass
                   datadir="data", # directory where files are stored
                   Focal = 2000., # focal length (cm)
                   Emin=100., # minimum of the nominal energy working band (keV)
                   Emax=600., # maximum of the nominal energy working band (keV)
                   dim=[1., 1., 0.4], # crystal tile (radial, tangential, thickness in cm)
                   framewidth = 0.20, # space around each tile (cm)
                   fwhm=1., # mosaic crystal FWHM (arcmin)
                   hkl=[1,1,1], # Miller indices
                   Z=29, # atomic number
                   microthick=0., # microblock thickness (A)
                   pixels=(50,50), # detector pixels (x,y)
                   pitches=(0.2,0.2), # detector pitch (x,y in cm)
                   error=0., # orientation error
                   )
# Running the program
lens.macro(INFO=True, # generic lens info
           XTALINFO=True, # crystal position and orientation
           ERANGE=(200,400), # energy range to be investigated
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