#!/usr/bin/env python
from LLL import Lenses
from math import *
# How to use the Generic lens object
# Lens definition
class Custom(Lenses.Generic):
    """This will be a class  to deal with a petal composed of rings"""
    def __init__(self,
                 name = "custom",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 Emin=100.,
                 Emax=600.,
                 dim=[1., 1., 0.4],
                 framewidth = 0.20,
                 fwhm=1.,
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        """Init of the class. The same of the parent Lens."""
        Lenses.Generic.__init__(self, name=name, geo="spiral", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, Emin=Emin, Emax=Emax,
                        dim=dim, framewidth=framewidth, fwhm=fwhm,
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)

    def xtals(self, COUNT=False):
        '''Custom crystal generator'''
        for xtal in self.gen_spiral_fixed(COUNT=COUNT):
            if 0<=self.r[1]*self.r[0]<=1000: yield None


lens=Custom()
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