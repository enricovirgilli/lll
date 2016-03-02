from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import
from pylab import plot, show, ylim, yticks
from LLL import Lenses
from LLL import __path__ as LLLpath
from pylab import *
from math import cos, sin
from pylab import plot, show, ylim, yticks, xlim
from scipy import integrate
from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import
import os, Physics
import sys

from Numeric import *

import Gnuplot, Gnuplot.funcutils
math=Physics
Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi
from numpy import arange, where
import numpy

Z=29
density={6 : 2.267,
         8 : 1.e-4,
        13 : 2.700,
        14 : 2.330,
        26 : 7.874,
        29 : 8.920,
        32 : 5.323,
        42 : 10.280,
        47 : 10.490,
        73 : 16.650,
        74 : 19.250,
        79 : 19.30,
        82 : 11.34,
        "air"  : 1.23e-3, # missing reference
        "SiO2" : 2.651, # missing reference
        "CdTe" : 6.0, # Ezio
        "GaAs" : 5.2, # Filippo
########################################################
         28 : 8.908,
         "InAs" : 5.680,
        }

EB=100
fwhm=3
microthick_micron=15

p_start = array([EB , fwhm, microthick_micron])


def muxxx(Z, keV):
    file = ('/home/ale/code/python/lib/LLL/29keV_all.dat')
    data = scipy.io.array_import.read_array(file)
    energy = data[:,0]
    ass = data[:,2]
    item=numpy.where(energy==int(keV))
    return ass[item]*density

######################################################################
# fitting of mu tables 
######################################################################
def poly(x,p):
    a = p[0]
    b = p[1]*x**(-1)
    c = p[2]*x**(-2)
    d = p[3]*x**(-3)
    e = p[4]*x**(-4)
    f = p[5]*x**(1)
    return a + b + c + d + e + f

def poliresid(p, y, x):
    err = y-poly(x,p) 
    return err

par_start=array([1, 2, 5000, 2e6, -8e6, 4])

file=('/home/ale/code/python/lib/LLL/data/cross_sections/29.dat')
data = scipy.io.array_import.read_array(file)

y = density[Z] * data[:,2]
#y1= y[10:len(y)]

x = 1000 * data[:,0]
#x1= x[10:len(x)]

plsq = leastsq(poliresid, par_start, args=(y, x), maxfev=2000)

A = plsq[0][0]
B = plsq[0][1]
C = plsq[0][2]
D = plsq[0][3]
E = plsq[0][4]
F = plsq[0][5]

def mu(Z,x):
    a = A
    b = B*x**(-1)
    c = C*x**(-2)
    d = D*x**(-3)
    e = E*x**(-4)
    f = F*x**(1)
    return a + b + c + d + e + f
##########################################################################






def reflectivity(Z, x, p): 
    return p[0] + p[1] * x + p[2] + mu(Z, x)

def residuals(p, y, x):
    err = y-reflectivity(Z,x,p) 
    return err


filename=('refletest.dat')
data = scipy.io.array_import.read_array(filename)

y = data[:,1]
keV = data[:,0]


EB=0.001
fwhm=0.3
microthick_micron=1




#p_start = array([EB , fwhm, microthick_micron])

plsq = leastsq(residuals, p_start, args=(y, keV), maxfev=2000)

print plsq

plot(keV,y,keV,reflectivity(Z,keV,plsq[0]))

#plot(keV,y,keV,reflectivity(Z,keV,p_start))

show()
