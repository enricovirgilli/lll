# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python
from LLL import Lenses
from pylab import *
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from numpy.fft import fft2, ifft2
from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import
from LLL import __path__ as LLLpath
from scipy.special import jn, erf
import scipy.io.array_import
import sys, Numeric, os, Physics, numpy
import Gnuplot, Gnuplot.funcutils
from myvec import *
import Booklet
from  Numeric  import zeros

Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi

rem_angstrom = 2.81794092e-5 # [electron radius in Angstrom]
hc=Booklet.hc

#Z='GaAs'                       # atomic number
Z = 14
hkl=(1,1,1)                    # miller indices
#hkl=(2,2,0)# miller indices
#hkl=(2,2,0)


external_radius = 40.  # meters
internal_radius = 2.46 * external_radius  # meters

curvature_s = 1 / internal_radius

#fwhm = 1. #in arcmin.

dim1=3.0    #
dim2=1.0    # dimensions of the cristal [cm]
dim3=0.1     # thickness

density = Booklet.density[Z]

lattice = Booklet.lattice[Z]
volume = Booklet.volume[Z]
TDebye = Booklet.TDebye[Z]
atomic_mass = Booklet.atomic_mass[Z]
sf = Booklet.sf(Z, hkl)

def d_hkl(Z, hkl):
    sum_square = sum(i**2 for i in hkl)
    return lattice / math.sqrt(sum_square)

d_hkl=d_hkl(Z, hkl)

print d_hkl

def fwhm2eta(x):
    return x/math.sqrt(8.*math.log(2.))

def mu(Z,keV):
    return Booklet.mu(Z, keV)


def keV2Bragg(keV):
    return math.asin(hc / (2. * d_hkl * keV))  #is given in radians

def Bragg2keV(angle):
    return hc / (2. * d_hkl *math.sin(angle))

def lambda0(keVB, keV):
    tB = keV2Bragg(keVB)
    up = math.pi  * math.cos(keV2Bragg(keVB))
    down =  (hc/keV) * rem_angstrom * (1 + abs(math.cos(2*tB)))
    v_by_f = volume / abs(sf)
    r = v_by_f * (up / down)
    return r

def peak_refl_curved(Z, keVB, keV, t):
    tB=keV2Bragg(keVB)
    cp = curvature_s/10000000000.0
    expo = (math.pi**2 * d_hkl)/(cp*lambda0(keVB,keV)**2)
    #diffraction = ( 1 - ( Numeric.exp(-(math.pi**2 * d_hkl) / (cp * lambda0(keVB,keV)**2) )) )
    diffraction = ( 1 - math.exp(-expo))
    absorption = Numeric.exp(- (mu(Z,keV)) * (t) / cos(tB))
    return diffraction * absorption

elist=arange(6, 600, 1)
for i, e in enumerate(elist):
    print e, peak_refl_curved(Z,e,e, dim3)


#########################################
#PSF Calculation:
##########################################



