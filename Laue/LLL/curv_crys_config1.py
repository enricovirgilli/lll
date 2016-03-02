# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python
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


math=Physics
Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi
from numpy import arange

rem_angstrom = 2.81794092e-5 # [electron radius in Angstrom]
#fwhm = 0.5 


curvature = 1.45e-2  # (1/m)
curv_radius = 100 #1 / curvature


Z=32                # atomic number
#microthick_micron=9   # microthickness [micron]

hkl=(1,1,1)            # miller indices

dim1=1.0     #
dim2=1.0     # dimensions of the cristal [cm]
dim3=1.0   # thickness

#fwhm = (dim3 * 180 * 60 ) / (internal_curv_radius * 100 * pi)  # mosaic spread [arcmin]

#omega = dim3 / (internal_curvature * 100)

fwhm = (dim3 * 180 * 60 ) / (curv_radius * 100 * pi)  # arcmin

EB = 400.   # peak energy of diffraction

alpha=dim3/(curv_radius * 100)    
    
theta_0=0
density = Booklet.density[Z]
lattice = Booklet.lattice[Z]
volume = Booklet.volume[Z]
TDebye = Booklet.TDebye[Z]
atomic_mass = Booklet.atomic_mass[Z]
hc=Booklet.hc

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
    up_for_lambda = volume * math.pi  #* math.cos(keV2Bragg(keVB))
    down_for_lambda =  (hc/keVB) * rem_angstrom * abs(sf) #* (1 + abs(math.cos(2*keV2Bragg(keVB))))
    lambda0 = up_for_lambda / ( down_for_lambda)
    omega = (fwhm * math.pi) / (60. * 180.) # conversion from arcmin to radians
    T = t * 100000000.0
    cp = omega / T
    muT = mu(Z,keV)
    mu_angstrom = muT/100000000.0
    diffraction = 1 - ( Numeric.exp(-(math.pi**2 * d_hkl * T) / (omega * lambda0**2) ))
    absorption = math.exp(-(mu_angstrom * T ) / (math.cos(keV2Bragg(keVB))))
    return diffraction * absorption


def lambda0(keVB, keV):
    up = volume * math.pi  * math.cos(keV2Bragg(keVB))
    down =  (hc/keV) * rem_angstrom * abs(sf) *  (1 + abs(math.cos(2*tB)))
    r = up / down
    return r

##########################################################################
def tbest(keVB, keV):
    omega = (fwhm / 60) * (math.pi / 180) # conversion from arcmin to radians
    #up_for_lambda = math.pi * volume * math.cos(keV2Bragg(keVB))
    #down_for_lambda =  (hc/keVB) * rem_angstrom * abs(sf) * (1 + abs(math.cos(2*keV2Bragg(keVB))))
    #lambda0 = up_for_lambda / down_for_lambda
    muT = mu(Z,keV)
    mu_angstrom = muT/100000000
    M=(math.pi**2 * d_hkl / lambda0(keVB,keV)**2)
    factor = (math.pi**2 * d_hkl * math.cos(keV2Bragg(keVB))) / (lambda0(keVB,keV)**2 * mu_angstrom * omega) 
    logarithm = log ( 1 + factor)
    return (omega * logarithm)/M
##########################################################################

#print tbest(EB, EB)/100000000.0

peak_gauss = peak_refl_curved(EB,EB, dim3)

#print math.pi**2 *d_hkl * dim3 * 100000000 /lambda0(100,100)**2

#print mu(Z,100)/(curvature * math.cos(keV2Bragg(100)))


    
    #f = (-A/x**2)*math.exp(-A/x) *math.exp(-B*x) - B* math.exp(-B*x) + B*math.exp(-A/x) * math.exp(-B*x)
    #f=math.exp(-A/x) * (A/x**2+B) - B
    #return f

#xlist=arange(0.0000001, 0.01, 0.0000001)
#for i, e in enumerate(xlist):
#    print e, function(e)

##############################################################
# 1) for Energy vs Peak reflectivity (at Tbest)
##############################################################
#elist=arange(10, 1000, 1)
#for i, e in enumerate(elist):
    #print e, tbest(e, e)/100000000, peak_refl_curved(e, e, tbest(e, e)/100000000)
    #t_best=tbest(e,e)/100000000
    #print e, peak_refl_curved(e,e, t_best)
##############################################################


###############################################################
# 2) for thickness vs peak reflectivity
###############################################################
#tlist=arange(0.1, 10, 0.1)
#for i, t in enumerate(tlist):
    #print t, peak_refl_curved(EB, EB, t)
###############################################################



###############################################################
# for Energy vs best thickness
###############################################################
#b=[]
#c=[]
#elist=arange(100, 500, 1)
#tlist=arange(0.1, 5, 0.005)
#for i, e in enumerate(elist):
#    for i, t in enumerate(tlist):
#        b.append(peak_refl_curved(e, e, t))
#        c.append(t)
#    mass=b.index(max(b))
#    tbest=c[mass]
#    print e, tbest
####################################################################

#enelist=arange(100, 400, 1)
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


