# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python
from LLL import Lenses
import Physics, Xtal
from Physics import *
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
import Booklet, Lenses
from Lenses import *
from  Numeric  import zeros

cryst_mat = Xtal.Material
cryst = Xtal.Xtal
Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi
volumes = Booklet.volume

hc=Booklet.hc
dim1=3.0    #
dim2=1.0    # dimensions of the cristal [cm]
dim3=01.0

Z = "GaAs"
hkl = [2,2,0]


density = Booklet.density[Z]

lattice = Booklet.lattice[Z]
volume = Booklet.volume[Z]
TDebye = Booklet.TDebye[Z]
atomic_mass = Booklet.atomic_mass[Z]
sf = Booklet.sf(Z, hkl)

def dhkl(Z, hkl):
    sum_square = sum(i**2 for i in hkl)
    print lattice, sum_square
    return lattice / math.sqrt(sum_square)

#d_hkl=dhkl(Z, hkl)

def mu(Z,keV):
    return Booklet.mu(Z, keV)




microthick = 12 #in microns


#def extinctionfactor(keV, order=1):
#    theta_0 = keV2Bragg(keV)
#    return Physics.extinction_factor(keV, Z, hkl[order], microthick, theta_0)

#print extinctionfactor(100, order=1)
#energy_keV = 30


keVrange = numpy.arange(50,600,1)
for i, energy_keV in enumerate(keVrange):
    theta_0 = 0#keV2Bragg(energy_keV)
    
    #print energy_keV, Physics.rockingcurve(energy_keV, energy_keV, Z, hkl, microthick, theta_0, eta, dim3, Tamorph=0.)

    #---Peak reflectivity----
    #------------------------
    #print energy_keV, reflectivity2(energy_keV, theta_0, Z, hkl, microthick, theta_0, eta, dim3, Tamorph=0.)


best_try = Lenses.EA_bendring_at_keVs(self, xtal_in_ring=1., NUM_OF_SIGMA=4)





#elist=arange(6, 600, 1)
#for i, e in enumerate(elist):
#    print e, peak_refl_curved(Z,e,e, dim3)


#########################################
#PSF Calculation:
##########################################



