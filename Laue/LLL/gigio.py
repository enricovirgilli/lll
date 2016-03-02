# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python
from LLL import Lenses
from pylab import *
from scipy import integrate
import os, Physics
math=Physics
Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi
from numpy import arange

Z=29
microthick_micron=20

microthick=microthick_micron*10000

fwhm=2     #mosaic spread [arcmin]
hkl=(1,1,1)
density = Booklet.density[Z]
lattice = Booklet.lattice[Z]
volume = Booklet.volume[Z]
TDebye = Booklet.TDebye[Z]
atomic_mass = Booklet.atomic_mass[Z]


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

g=20
theta_0=0

D=2000
dim1=1.5
dim2=1.5
dim3=0.3

g=0
gnormalized = Physics.normalize(g)
gtheta = 0.0


EB = 100.0

def extinction_factor(keVB, theta_0):
    return Physics.extinction_factor(keVB, Z, hkl, microthick, theta_0)


def getDelta(self, tB, OFFSET=0.):
    return pi/2-gtheta-tB


mat_facs = Booklet.mat_fac(Z, hkl)


def angular_factor(keV, keVB):
   return 2 * math.sin(keV2Bragg(keV))**3 * (1 + math.cos(2 * keV2Bragg(keVB))**2 ) /  (2 * math.sin( 2 * keV2Bragg(keVB)))


def Q(keV, keVB):
    return mat_facs * angular_factor(keV, keVB)

#def squarefn()

def gaussian(x, eta, xB):
    """Normalized Gaussian distribution function centered in 0. by default."""
    delta=((x-xB)/eta)
    exponential=math.exp(-delta**2/2.)
    normalization = math.sqrt(2.*pi) * eta
    return exponential/normalization

def sigma(keVB, keV):
    t = keV2Bragg(keV)
    tB = keV2Bragg(keVB)
    #weight=gaussian(t, eta, tB)
    weight=fwhm
    ext_fa=extinction_factor(keVB, theta_0)
    return Q(keV, keVB)* weight*ext_fa

def reflectivity(keVB, keV):
    T = dim3
    muT= mu(Z,keV) * T
    sigmaT = abs(sigma(keV, keVB) * T)
    r = 0.5 * (1 - math.exp(-2 * sigmaT)) * math.exp(-muT)
    return r


#print reflecdiv(100)

def keVgen(keVB, NUM_OF_SIGMA=4):
    """Generates an energy list"""
    Dthetamax = NUM_OF_SIGMA*eta
    tB = keV2Bragg(keVB)
    Emin = Bragg2keV(tB + Dthetamax)
    Emax = Bragg2keV(tB - Dthetamax)
    s=arange(Emin*100,Emax*100)
    return s/100


refl = array([reflectivity(e, EB) for e in keVgen(EB, NUM_OF_SIGMA=4)])


#def reflecdiv(keV):

s = keVgen(EB, NUM_OF_SIGMA=4)

#size_s = size(s)
#divrefl = [None]*size_s
#for j in range (0,size_s):
#    EB = s[j]
#    tBmin= keV2Bragg(EB)-dim1/(2*D)
#    tBmax= keV2Bragg(EB)+dim1/(2*D)
#    EBmin=Bragg2keV(tBmax)
#    EBmax=Bragg2keV(tBmin)

#    x=arange(EBmin*100,EBmax*100)/100
#    size_x = size(x)
#    result = [None]*size_x
#    for i in range (0,size_x):
#        result[i] = reflectivity(x[i], EB)

    #divrefl = [None]*size_x
#    divrefl[j] = integrate.trapz(result)
   
#print divrefl    

#divrefl = array([reflecdiv(e) for e in keVgen(EB, NUM_OF_SIGMA=4)])

#s = keVgen(EB, NUM_OF_SIGMA=4)

plot(s, refl)

show()


