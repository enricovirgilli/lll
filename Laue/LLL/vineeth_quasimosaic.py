# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python

from pylab import *
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from numpy.fft import fft2, ifft2
from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import
from scipy.special import jn, erf
import scipy.io.array_import
import sys, Numeric, os, Physics, numpy
import Gnuplot, Gnuplot.funcutils
from myvec import *
import Booklet
from  Numeric  import zeros
from math import *
from Booklet import *
from scipy.special import jn, erf
from numpy import array
from random import uniform
import sys
from myvec import *
import Booklet
import Xtal
import Numeric

hc=12.39842 # [keV x �]
rem_angstrom = 2.81794092e-5 # [�]
keVB = 300.
hkl=[2,2,0]
Z = 14
d_hkl = Booklet.d_hkl(Z, hkl)
sf = Booklet.sf(Z, hkl)

volume = Booklet.volume[Z]
tB = Physics.BraggAngle(keVB, d_hkl)
omega_arcsec = 40.0
omega = omega_arcsec / 3600 * (math.pi/180)
R_ext = 40. # in meters
R_int = R_ext * 2.5 # in meters
base_thickness_t0 = 0.1 # in micron
Total_thickness = 4. #iim mm

basic_thickness = base_thickness_t0 * 10.e4 #in angstrom
t_n = Total_thickness * 10.e7 #in angstrom
R_int_angstr =  R_int * 10.e10 #in angstrom
n = 1

angle = math.pi/2 - tB

ang_alpha = angle  #in radiance
ang_beta = angle  #in radiance
gamma_0 = sin(ang_alpha)
gamma_H = sin(ang_beta) #in radiance, and in Laue case, gamma_0 = gamma_H 


def mu0(Z, keV):
    return Booklet.mu(Z, keV)

def gamma0_n(n):
    return sin( ang_alpha + (n * basic_thickness / R_int_angstr ) )

def gammaH_n(n):
    return sin( ang_beta + (n * basic_thickness / R_int_angstr ) )

def psi_H_first(keV, Z, hkl, th0=0.):
    num = abs(sf) * (hc**2) 
    den= (keV**2) * math.pi * volume
    g = -(num/den)
    return g

def y(keV, Z, hkl, th0=0.):
    #num = 4.0 * keVB * (sin(tB))**2.0 * math.pi * volume * (-keV)**2
    #den = (keV-keVB) * 1.0 * abs(1+cos(2.*tB)) * abs(sf) * hc**2
    num = alpha(keV, Z, hkl, th0=0.)/2
    den = abs((1+cos(2.*tB))/2) * abs(psi_H_first(keV, Z, hkl, th0=0.))
    a = (num/den)
    return a


#print y(101, Z, hkl, th0=0.)
def r(keV, Z, hkl, th0=0.):
    m =  abs(y(keV, Z, hkl, th0=0.)**2 - 1) +  y(keV, Z, hkl, th0=0.)**2
    return m - sqrt(m**2-1)

def qpz2(keV, Z, hkl, th0=0.):
    temp = (((1 + cos(2*tB))/2)**2) * (((((hc/keV)**2) * sf)/(math.pi * volume))**2) * (1 + y(keV, Z, hkl, th0=0.)**2)
    return abs(temp)

def Ag(keV, Z, hkl, n, th0=0.):
    num = basic_thickness * (-mu0(Z, keV))*(hc/keV)
    den = d_hkl * sin(tB) * math.sqrt(gamma0_n(n)*gammaH_n(n))
    return num/den

def A(keV, Z, hkl, n, th0=0.):
    num = (((1 + cos(2*tB))/2)) * (hc**2) * sf * basic_thickness
    den = d_hkl * sin(tB) * (keV**2) * volume * math.sqrt(gamma0_n(n)*gammaH_n(n))
    return -num/den

def q(keV, Z, hkl, th0=0.):
    temp = ((((hc/keV)**2) * sf/(math.pi*volume) )**2) * rem_angstrom**2
    return abs(temp)

def alpha(keV, Z, hkl, th0=0.):
    temp = (4 * keVB * sin(tB)**2) / ((keV - keVB))
    return temp


def D(keV, Z, hkl, n, th0=0.):
    one = abs(qpz2(keV, Z, hkl, th0=0.))
    two = (one - (alpha(keV, Z, hkl, th0=0.)/2)**2)*sin(abs((A(keV, Z, hkl, n, th0=0.)* math.sqrt(1 + y(keV, Z, hkl, th0=0.)**2))))
    three = 0.5 * (math.sqrt(abs((one - (alpha(keV, Z, hkl, th0=0.)/2)**2) -  q(keV, Z, hkl, th0=0.)**2)) * sin(abs(2 * A(keV, Z, hkl, n, th0=0.) * math.sqrt(1 + y(keV, Z, hkl, th0=0.)**2))))
    return one - two + three

def tk(keV, Z, hkl, n, th0=0.):
    num = qpz2(keV, Z, hkl, th0=0.) * math.exp(Ag(keV, Z, hkl, n, th0=0.))
    return num / D(keV, Z, hkl, n, th0=0.)

def exp_rn(Z, keV, n):
    return math.exp( mu0(Z, keV) * basic_thickness * (10e-8) / gammaH_n(n) )




def refl_n(keV, Z, hkl, th0=0.):
    n = 50#int( math.floor( t_n / basic_thickness))
    r_n = 0
    temp_r = r(keV, Z, hkl, th0=0.)
    #print temp_r
    for i in xrange(1,n):
        r_n = temp_r / ((1 - tk(keV, Z, hkl, n, th0=0.)**i ) * exp_rn(Z, keV, n))
        temp_r = r_n
    refl = r_n
    return refl



for g in range(100, 2999):
    i = g/10.
    #print y(i, Z, hkl, th0=0.), r(i, Z, hkl, th0=0.)
    print i, abs(refl_n(i, Z, hkl, th0=0.))

for g in range(3001,8000 ):
    i = g/10.
    #print y(i, Z, hkl, th0=0.), r(i, Z, hkl, th0=0.)
    print i, abs(refl_n(i, Z, hkl, th0=0.))
    #print i, abs(r(i, Z, hkl, th0=0.))
