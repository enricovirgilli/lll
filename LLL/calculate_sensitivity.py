# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python

from LLL import Lenses
from pylab import *
from pylab import plot, show, ylim, yticks, xlim
from math import cos, sin, asin
from scipy import integrate
from scipy.interpolate import spline
from scipy import *
from scipy.optimize import leastsq
from scipy.special import jn, erf
from numpy.fft import fft2, ifft2
from numpy import arcsin, arange, where

from LLL import __path__ as LLLpath

from Numeric import *
from myvec import *

import os, Physics
import matplotlib.pyplot as plt
import numpy as np
import scipy.io.array_import
import Numeric
import Gnuplot, Gnuplot.funcutils
import scipy.io.array_import
import pylab, numpy, scipy, sys
import numpy

math=Physics
Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi

#rem_angstrom = 2.81794092e-5 # [electron radius in Angstrom]

#external_curvature = 40.0  # meters
#Z=14                  # atomic number
#microthick_micron=9   # microthickness [micron]
#fwhm=1.            # mosaic spread [arcmin]
#hkl=(1,1,1)            # miller indices

#dim1=1.5     #
#dim2=1.5     # dimensions of the cristal [cm]
#dim3=0.5   #  thickness

#EB = 300.   # peak energy of diffraction

#alpha=dim1/(external_curvature * 100)    
    
#theta_0=0

Z=14
density = Booklet.density[Z]
lattice = Booklet.lattice[Z]
volume = Booklet.volume[Z]
TDebye = Booklet.TDebye[Z]
atomic_mass = Booklet.atomic_mass[Z]
hc=Booklet.hc


filename ="EA_lens.dat" 

e=[float(x.split()[0]) for x in file("EA_lens.dat")]
elow=e[0]
eup = e[len(e) - 1]

enebottom=arange(0, elow-1, 1)

ee = []
aa = []
for i in range(0,int(elow)):
    ee.append(i)
    aa.append(100)

e=[float(x.split()[0]) for x in file("EA_lens.dat")]
a=[float(x.split()[1]) for x in file("EA_lens.dat")]

eee = ee + e
aaa = aa + a

def area(keV):
    return aaa[int(keV)]


def Sfactor_at_keV(keV, fraction=0.8):
    radius = 0.1
    Ad=math.pi*radius**2
    return (fraction*area(keV))**2/Ad


def sensitivity_at_keV(keV, fraction=0.5, Tobs=1.e5, nsigma=3., effdet=0.9, BG="default"):
    eqS=Sfactor_at_keV(keV, fraction=fraction)+2.05641071127e-10
    #DeltaE=0.2 * keV + 2.05641071127e-15 #0.5*keV+2.05641071127e-10
    DeltaE=0.5 * keV + 2.05641071127e-15
    BG = Booklet.bg(keV)/2 
    return (nsigma/effdet)*math.sqrt(2. * BG / (Tobs*DeltaE*eqS))


enelist=arange(0, eup, 1)
for i, e in enumerate(enelist):
    print e, e*e*sensitivity_at_keV(e, fraction=0.5, Tobs=1.e5, nsigma=3., effdet=0.9, BG="default")



sys.exit(0)

#print sensitivity_at_keV(250.0, fraction=0.5, Tobs=1.e6, nsigma=3., effdet=1., BG="default")


    #if eqS==0: return 1.
    #if BG=="default":
    #    BG = Booklet.bg(keV)
    #elif BG=="ezio":
    #    BG = Booklet.spibg.read(keV)/Booklet.mu("CdTe", keV)
    #    return nsigma/effdet*math.sqrt(2.*BG / (Tobs*DeltaE*eqS))


#print sensitivity_at_keV(self, keVs, areaeff=None, fraction=0.5, radius=None, Tobs=1.e6, nsigma=3., effdet=1., BG="default")



def fwhm2eta(x):
    return x/math.sqrt(8.*math.log(2.))

eta = Booklet.arcmin2rad(Physics.fwhm2eta(fwhm)) # radians

def mu(Z,keV):
    """Wrapper function"""
    return Booklet.mu(Z, keV)

sf = Booklet.sf(Z, hkl) # []

d_hkl = Booklet.d_hkl(Z, hkl) # [Angstrom]



def keV2Bragg(keV):
    return math.asin(Booklet.hc / (2. * d_hkl*keV))  #is given in radians

def Bragg2keV(theta):
    return Booklet.hc / (2. * d_hkl*math.sin(theta))

tB = keV2Bragg(EB)
Emin = Bragg2keV(tB + alpha/2)   # min and max energy diffracted by this curved
Emax = Bragg2keV(tB - alpha/2)   # crystal due to the dim[1] and to curvature radius

#print keV2Bragg(Emin)/math.pi*180*60*60-keV2Bragg(Emax)/math.pi*180*60*60
#print keV2Bragg(Emax)/math.pi*180*60*60

#print tB + alpha/2
#print tB - alpha/2
#print tB
def best_thickness(keV):
    """Wrapper for the Physics module function"""
    return Physics.best_thickness(keV, Z, hkl, eta)


def extinction_lenght(keV, Z, hkl, th0=0.):
    a = ( volume[Z] * keV * cos(th0) ) / rem_angstrom * hc * sf(Z, hkl)
    return a

def extinction_factor(keVB, theta_0):
    return Physics.extinction_factor(keVB, Z, hkl, microthick, theta_0)

mat_facs = Booklet.mat_fac(Z, hkl)  


def gaussian(x, eta, xB):
    """Normalized Gaussian distribution function centered in 0. by default."""
    delta=((x-xB)/eta)
    exponential=math.exp(-delta**2/2.)
    normalization = math.sqrt(2.*pi) * eta
    return exponential/normalization

def sigma(keVB, keV):
    t = keV2Bragg(keV)
    tB = keV2Bragg(keVB)
    weight=gaussian(t, eta, tB)
    ext_fa=extinction_factor(keVB, theta_0)
    return Q(keV, keVB)* weight * ext_fa


def peak_refl_curved(keVB, keV, t):
    up_for_lambda = math.pi * volume * math.cos(keV2Bragg(keVB))
    down_for_lambda =  (hc/keVB) * rem_angstrom * abs(sf) * (1 + abs(math.cos(2*keV2Bragg(keVB))))
    lambda0 = up_for_lambda / down_for_lambda 
    omega = (fwhm / 60) * (math.pi / 180) # conversion from arcmin to radians
    T = t * 100000000.0
    cp = omega / T
    muT = mu(Z,keV)
    mu_angstrom = muT/100000000
    diffraction = 1 - math.exp(-(math.pi**2 * d_hkl) / (cp * lambda0**2) )
    absorption = math.exp(-(mu_angstrom * omega ) / (cp * math.cos(keV2Bragg(keVB))))
    return diffraction #* absorption

peak_gauss = peak_refl_curved(EB,EB, dim3)



##############################################################
# 1) for Energy vs Peak reflectivity (at Tbest)
##############################################################
a=[]
elist=arange(100, 500, 1)
tlist=arange(0.1, 10, 0.1)
for i, e in enumerate(elist):
    for i, x in enumerate(tlist):
        a.append(peak_refl_curved(e, e, x))
#    print e, max(a)
##############################################################


###############################################################
# 2) for thickness vs peak reflectivity
###############################################################
tlist=arange(0.1, 10, 0.1)
#for i, t in enumerate(tlist):
#    print t, peak_refl_curved(EB, EB, t)
###############################################################



###############################################################
# for Energy vs best thickness
###############################################################
b=[]
c=[]
elist=arange(100, 500, 1)
tlist=arange(0.1, 5, 0.005)
#for i, e in enumerate(elist):
#    for i, t in enumerate(tlist):
#        b.append(peak_refl_curved(e, e, t))
#        c.append(t)
#    mass=b.index(max(b))
#    tbest=c[mass]
#    print e, tbest
####################################################################

enelist=arange(100, 500, 1)
#for i, e in enumerate(enelist):
#    print e, peak_refl_curved(e,e, dim3)
################################################################


    
def keVgen(keVB, NUM_OF_SIGMA):
    DeltaEnergy = NUM_OF_SIGMA * (Emax-Emin)
    Einf = Emin - DeltaEnergy/2
    Esup = Emax + DeltaEnergy/2
    s=arange(Einf*200,Esup*200)
    return s/200

def rectfun(x):
    if  Emin < x < Emax : return peak_refl_curved(EB, x, dim3)
    else:
        return 0

def smooth(x,window_len=11,window='hanning'):
    s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat':
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]


rectarray = array([rectfun(e) for e in keVgen(EB, 2)])

f = keVgen(EB, NUM_OF_SIGMA=2)

reflecurved = smooth(rectarray,window_len=EB**1.09-85,window='hanning')


plot (f,reflecurved, label='line 1')
plot (f,rectarray)
plt.show()


