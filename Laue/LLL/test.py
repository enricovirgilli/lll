from LLL import Lenses
from pylab import *
from math import cos, sin
from pylab import plot, show, ylim, yticks
from scipy import integrate
import os, Physics

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import spline
import pylab,numpy
from numpy.fft import fft2, ifft2
from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import
from pylab import plot, show, ylim, yticks
from LLL import Lenses
from LLL import __path__ as LLLpath
from pylab import *
import pylab
from math import cos, sin, asin
from pylab import plot, show, ylim, yticks, xlim
from scipy import integrate
import scipy
from scipy.special import jn, erf
from numpy import arcsin
from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import
import os, Physics
import sys
from Numeric import *
import Numeric
import Gnuplot, Gnuplot.funcutils
from myvec import *
math=Physics
Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi
from numpy import arange, where
import numpy


def EA_ring_at_keVs(self, xtal_in_ring=1., NUM_OF_SIGMA=4):
    tB = math.pi/2-self.gtheta # always on axis
    max_order = int(self.keVlist[-1]/self.Bragg2keV(tB)+1)
    ring_cross_section=self.cross_section * xtal_in_ring
    for order in range(1,max_order):
        if not self.sfs[order]==0.:
            Emin, Emax = self.Erange(order, OFFSET=0., NUM_OF_SIGMA=NUM_OF_SIGMA)
            keVs=[(i,keV) for i,keV in enumerate(self.keVlist) if Emin<=keV<=Emax]
            poss=[(i, self.reflectivity(keV, 0., order)) for i, keV in keVs]
            for i, pos in poss:
                if pos!=0.: self.arealist[i] = self.arealist[i] + pos * ring_cross_section



def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4):
        self.arealist = zeros(len(self.keVlist), 'f')
        self.EA_ring_at_keVs(self.noofxtals, NUM_OF_SIGMA)



class Rings(Generic):
    """This will be a class  to deal with systems of rings"""
    def __init__(self,
                 name = "rings",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 Emin=100.,
                 Emax=600.,
                 dim=[1.5, 1.5, 0.4],
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
        Generic.__init__(self, name=name, geo="rings", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, Emin=Emin, Emax=Emax,
                        dim=dim, framewidth=framewidth, fwhm=fwhm,
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)
        self.noofrings=1+max([x[0] for x in enumerate(self.rings(COUNT=True))])

    def __str__(self, filename=False):
        """Displays ottica values"""
        dim=(self.dim[0], self.dim[1], self.averagethickness)
        string=("Geometry:                 %s" % self.geo,
                "Profile:                  %s" % self.profile,
                "Thickness profile:        %s" % self.thickness_profile,
                "Thickness threshold:      %s" % repr(self.thickness_threshold),
                "Thickness factor:         %s" % self.thickness_factor,
                "Focal Length:             %s m" % (self.Focal/100),
                "Inner/Outer radius:       %s/%s m" % (self.rmin/100, self.rmax/100),
                'Z (Atomic number):        %s' % self.Z,
                "Energy range:             [%g, %g] keV" % (self.Emin, self.Emax),
                'hkl:                      (%d,%d,%d)' % tuple(self.hkl),
                "d_hkl:                    %s A" % Booklet.d_hkl(self.Z, self.hkl),
                'Cell volume:              %s A^3' % Booklet.volume[self.Z],
                'Structure factor:         %s' % self.sf,
                'Microcrystal size:        %s micron' % (self.microthick/10000),
                'A_0 constant:             %s keV' % self.A_0const,
                'Material factor constant: %s cm^-1' % self.mat_fac,
                'Crystal tile size:        %gx%gx%g cm^3' % tuple(dim),
                'Number of crystals:       %s' % self.noofxtals,
                'Filling factor:           %s' % self.filling_factor,
                'Number of rings:          %s' % self.noofrings,
                'Crystals volume:          %s cm^3' % self.volume,
                'Crystals weight:          %s kg' % self.weight,
                'mosaicity (FWHM):         %s arcmin' % self.fwhm,
                'Frame width:              %s cm' % self.framewidth,)
        return "\n".join(string)

    def rings(self, COUNT=False):
        '''Generator for returning crystals that are placed on concentric rings
        to diffract energies in the nominal range [Emin, Emax]'''
        keV=float(self.Emax)
        while(keV >= self.Emin):
            # yield xtals of ONE ring
            xtal_in_ring = max(x[0] for x in enumerate(self.gen_one_ring(keV=keV, COUNT=COUNT)))+1
            # Next ring: radius and energy calculation (PhD thesis)
            part1 = (self.rho+self.dim[0]/2.+self.framewidth)**2
            part2 = (self.dim[1]/2.+self.framewidth)**2
            sqrt = math.sqrt(part1+part2)
            next_rho = self.dim[0]/2. + self.framewidth + sqrt
            # next energy
            keV = self.radius2keV(next_rho)
            yield xtal_in_ring

    def rings_info(self):
        """returns info about rings structure"""
        return [(ring, self.rho, self.radius2keV(self.rho), self.dim[2]) for ring in self.rings()]

    def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4):
        """Calculates the effective area for a single xtal and then multiplies
        by the number of xtals for each ring."""
        self.arealist = zeros(len(self.keVlist), 'f')
        # monitoring tools
        countercondition=(self.noofxtals>1000)
        if countercondition:
            counter=Progressbar(update_time=1.,
                               steps=self.noofxtals,
                               title=self.name,
                               name='Effective area',
                               )
        # effective calculation
        if OFFSET==0. and self.error==0.:
            for xtal_in_ring in self.rings():
                self.EA_ring_at_keVs(xtal_in_ring, NUM_OF_SIGMA)
                if countercondition: counter.index+=xtal_in_ring
        else:
            for xtal in self.xtals():
                self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)
                if countercondition: counter.index+=1
        if countercondition: counter.close()
