# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python
from LLL import Lenses
from pylab import *
from math import cos, sin
from pylab import plot, show, ylim, yticks
from scipy import integrate
import os, Physics
math=Physics
Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi
from numpy import arange

rem_angstrom = 2.81794092e-5 # [electron radius in Angstrom]

Z=29                   # atomic number
microthick_micron=15   # microthickness [micron]
fwhm=3.0              # mosaic spread [arcmin]
hkl=(1,1,1)            # miller indices

D=2000       # distance source-cristal [cm]

dim1=1.5     #
dim2=1.5     # dimensions of the cristal [cm]
dim3=0.3     #

EB = 104.0   # peak energy of diffraction

theta_0=0
density = Booklet.density[Z]
lattice = Booklet.lattice[Z]
volume = Booklet.volume[Z]
TDebye = Booklet.TDebye[Z]
atomic_mass = Booklet.atomic_mass[Z]
hc=Booklet.hc

microthick=microthick_micron*10000
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

#def keVB2Bragg(keVB):
#    return math.asin(Booklet.hc / (2. * d_hkl*keVB))  #is given in radians

def Bragg2keV(theta):
    return Booklet.hc / (2. * d_hkl*math.sin(theta))


def best_thickness(keV):
    """Wrapper for the Physics module function"""
    return Physics.best_thickness(keV, Z, hkl, eta)


def extinction_lenght(keV, Z, hkl, th0=0.):
    a = ( volume[Z] * keV * cos(th0) ) / rem_angstrom * hc * sf(Z, hkl)
    return a

def extinction_factor(keVB, theta_0):
    return Physics.extinction_factor(keVB, Z, hkl, microthick, theta_0)

mat_facs = Booklet.mat_fac(Z, hkl)  

def angular_factor(keV, keVB):
    tetaB = keV2Bragg(keVB)
    teta =  keV2Bragg(keV)
    return sin(tetaB)**2  *  ( 1. + (cos(2.*tetaB) )**2. ) / cos(tetaB)
    #return 2 * sin(tetaB)**3 * (1 + (cos(2*tetaB))**2  ) / ( 2 * sin(2*tetaB))
     
    #return 2 * math.sin(keV2Bragg(keV))**3 * (1 + (math.cos(2 * keV2Bragg(keVB))**2 ) /  (2 * math.sin( 2 * keV2Bragg(keVB)))
   
def Q(keV, keVB):
    return mat_facs * angular_factor(keV, keVB)


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

def reflectivity(keVB, keV):
    T = dim3
    muT= mu(Z,keV) * T
    sigmaT = abs(sigma(keV, keVB) * T)
    r = 0.5 * (1 - math.exp(-2 * sigmaT)) * math.exp(-muT)
    return r


def peak_refl_curved(keVB, keV):

    up_for_lambda = math.pi * volume * math.cos(keV2Bragg(keVB))
    down_for_lambda =  (hc/EB) * rem_angstrom * sf * (1 + math.cos(2*keV2Bragg(keVB)))
    
    lambda0 = up_for_lambda / down_for_lambda 
    
    sigma = (fwhm / 60) * (math.pi / 180) # conversion from arcmin to radians
    
    cp = sigma / (dim2 * 100000000)

    T = dim3 / 100000000

    muT = mu(Z,keV) * T

    diffraction = 1 - math.exp(-(math.pi**2 * d_hkl) / (cp * lambda0**2) )
    
    absorption = math.exp(-(muT * sigma) / cp * math.cos(keV2Bragg(keVB)))

    #return math.exp(-(math.pi**2 * d_hkl) / (cp * lambda0**2) )
    
    return lambda0 /10000

    
print peak_refl_curved(EB,100)

def f(keVB, keV):
    T = dim3
    muT= mu(Z,keV) * T
    sigmaT = abs(sigma(keV, keVB) * T)
    ratio=hc/(2*d_hkl*keVB)
    radice=1/math.sqrt(1-ratio**2)
    g = 0.5 * (1 - math.exp(-2 * sigmaT)) * math.exp(-muT) * hc/(2*d_hkl) * radice / keVB**2
    return g


def integrefle(keV):

    n=100
    
    tBmin= keV2Bragg(EB)-dim1/(2*D)
    tBmax= keV2Bragg(EB)+dim1/(2*D)
    EBmin=Bragg2keV(tBmax)
    EBmax=Bragg2keV(tBmin)    

    h = float (EBmax - EBmin)/n

    #s = 0.5 * ( reflectivity(EBmin, keV) +  reflectivity(EBmax, keV))
    s = 0.5 * ( f(EBmin, keV) +  f(EBmax, keV))
    for i in range (1,n):
        #s=s+reflectivity(EBmin+i*h, keV) * h
        s=s+f(EBmin+i*h, keV) * h
    return s * D / dim1 


def keVgen(keVB, NUM_OF_SIGMA=8):
    """Generates an energy list"""
    Dthetamax = NUM_OF_SIGMA*eta
    tB = keV2Bragg(keVB)
    Emin = Bragg2keV(tB + Dthetamax)
    Emax = Bragg2keV(tB - Dthetamax)
    
    s=arange(Emin*20,Emax*20)
    return s/20


inte = array([integrefle(e) for e in keVgen(EB, NUM_OF_SIGMA=8)])
refle = array([reflectivity(EB, e) for e in keVgen(EB, NUM_OF_SIGMA=8)])

f = keVgen(EB, NUM_OF_SIGMA=8)

#plot(f, inte, f, refle)
plot(f, refle)
show()


