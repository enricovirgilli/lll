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

keV = 100.
hkl=[2,2,0]
Z = "GaAs"
d_hkl = Booklet.d_hkl(Z, hkl)
sf = Booklet.sf(Z, hkl) 
volume = Booklet.volume[Z]
tB = Physics.BraggAngle(keV, d_hkl)
omega_arcsec = 80.0
omega = omega_arcsec / 3600 * (math.pi/180)
etaarcsecondi=30.
eta= 4.848*0.000001*etaarcsecondi
mu=Booklet.mu(Z, keV)

def extinction_constant(Z, hkl, th0=0.):
    """Material dependent constant useful for primary extinction calculation. Unit is keV/A"""
    return rem_angstrom * hc * abs(sf) / ( math.pi * volume * cos(tB) )

def extinction_lenght(keV, Z, hkl, th0=0.):
    """Material and energy dependent constant useful for primary extinction"""
    return keV/ (extinction_constant(Z, hkl, th0=0.))


def lambda_zero(keV, Z, hkl, th0=0.):
    #tB = BraggAngle(keV, d_hkl(Z, hkl))
    ex_le =  extinction_lenght(keV, Z, hkl, th0=0.)
    polarization = (1 + abs(cos(2*tB)))/2
    return ex_le/polarization

def lambda_curvo(keV, Z, hkl, th0=0.):
    ex_le =  extinction_lenght(keV, Z, hkl, th0=0.)
    polarization = (1 + abs(cos(2*tB)))/2
    l = ex_le/polarization
    r = (l*l*omega)/(math.pi*math.pi*d_hkl)
    return r

#print omega
#print lambda_zero(keV, Z, hkl, th0=0.)/10000000
#print lambda_curvo(keV, Z, hkl, th0=0.)/10000000

lc=lambda_curvo(100, Z, hkl, th0=0.)/100000000
#print lc, 2*lc, 3*lc
for i in range(int(keV-keV/2),int(keV+keV/2)):
    t0=i/100.0
    microthick=50
    th0=0
    n=3.0
    t = n * lc
    mu=Booklet.mu(Z, i)
    T=2.5
    lc1=lambda_curvo(i, Z, hkl, th0=0.)/100000000
    #Physics.best_thickness(i, Z, hkl, eta)
    #print  t0, exp(-t0*Booklet.mu(Z, 300.0))
    #print  i,lambda_curvo(i, Z, hkl, th0=0.)/10000000
    ######print i, 1-math.exp(-t/lc1), math.exp(-Booklet.mu(Z, i)*t) 
#print  i,lambda_zero(i, Z, hkl, th0=0.)/10000000
    #print i, Physics.best_thickness(i, Z, hkl, eta)
    #print i, Physics.rockingcurve(i, keV, Z, hkl, microthick, th0, eta, T, Tamorph=0.)


for i in range(50, 300):
    
    print  i,4*(lambda_zero(i, Z, hkl, th0=0.)/10000000000)/(math.pi*4.8e-6)
    
