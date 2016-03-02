from scipy import *
from scipy.optimize import leastsq
import scipy.io.array_import
from pylab import plot, show, ylim, yticks
from LLL import Lenses
from LLL import __path__ as LLLpath
from pylab import *
from math import cos, sin, asin
from pylab import plot, show, ylim, yticks, xlim
from scipy import integrate
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
        "air"  : 1.23e-3,
        "SiO2" : 2.651,
        "CdTe" : 6.0,
        "GaAs" : 5.2,
         28 : 8.908,
         "InAs" : 5.680,
        }

fwhm=4
EB=99.2
microthick_micron=40   # microthickness [micron]

rem_angstrom = 2.81794092e-5 # [electron radius in Angstrom]

Z=29                   # atomic number


hkl=(1,1,1)            # miller indices

D=2000       # distance source-cristal [cm]

dim1=1.5     #
dim2=1.5     # dimensions of the cristal [cm]
dim3=0.3     #

theta_0=0

lattice = Booklet.lattice[Z]
volume = Booklet.volume[Z]
TDebye = Booklet.TDebye[Z]
atomic_mass = Booklet.atomic_mass[Z]
hc=Booklet.hc

microthick=microthick_micron*10000

p_start = array([fwhm, microthick])
#def arcmin2rad(theta):
#    """arcminute to radian conversion"""
#    return theta*math.pi/60./180.

#def fwhm2eta(x):
#    return x/math.sqrt(8.*math.log(2.))

#def eta(p):
#    return arcmin2rad(fwhm2eta(p[0]))


def  eta(p):
    return  ( p[0] * math.pi / (180 * 60) ) / (math.sqrt(8.*math.log(2.)) )

######################################################################
# fitting of mu tables 
######################################################################
def poly(x,v):
    a = v[0]
    b = v[1]*x**(-1)
    c = v[2]*x**(-2)
    d = v[3]*x**(-3)
    e = v[4]*x**(-4)
    f = v[5]*x**(1)
    return a + b + c + d + e + f

def poliresid(v, y, x):
    err = y-poly(x,v) 
    return err

var_start=array([1, 1, 1, 1, 1, 1])

file=('/home/ale/code/python/lib/LLL/data/cross_sections/29.dat')
data = scipy.io.array_import.read_array(file)

ymu = density[Z] * data[:,2]
y1= ymu[10:4000]

xmu = 1000 * data[:,0]
x1= xmu[10:4000]

best = leastsq(poliresid, var_start, args=(y1, x1), maxfev=2000)

A = best[0][0]
B = best[0][1]
C = best[0][2]
D = best[0][3]
E = best[0][4]
F = best[0][5]

def mu(Z,x):
    a = A
    b = B*x**(-1)
    c = C*x**(-2)
    d = D*x**(-3)
    e = E*x**(-4)
    f = F*x**(1)
    return a + b + c + d + e + f

filename=('reflectivity_dx.dat')
data = scipy.io.array_import.read_array(filename)
absorp = data[:,1]
energy = data[:,0]

##########################################################################


sf = Booklet.sf(Z, hkl) # []


d_hkl = Booklet.d_hkl(Z, hkl) # [Angstrom]


def keV2Bragg(keV):
    return arcsin(Booklet.hc / (2. * d_hkl*keV))  #is given in radians


def Bragg2keV(theta):
    return Booklet.hc / (2. * d_hkl*math.sin(theta))


#############################################################
def BraggAngle(keV, d_hkl):
    return arcsin(hc / (2*d_hkl*keV) )

def sumoddbess(a, lower_limit=1.e-10):
    tmp, i =0., 0
    ctrl1, ctrl2 = 1., 1.
    while (abs(ctrl1)+abs(ctrl2) >= lower_limit):
        ctrl2=ctrl1
        ctrl1 = jn(2*i+1, a)
        tmp+=ctrl1
        i+=1
    return tmp

def extinction_constant(Z, hkl, th0=0.):
    stru_fac = abs(sf) 
    return rem_angstrom * hc * stru_fac / ( volume * cos(th0) )

def extinction_lenght(keVB, Z, hkl, th0=0.):
    return keVB/extinction_constant(Z, hkl, th0=0.)

def extinction_parameter(keVB, Z, hkl, p, th0=0.):
    #tB = BraggAngle(keV, d_hkl)
    #polarization = 1 + abs(cos(2*tB))
    return p[1] / extinction_lenght(keVB, Z, hkl, th0)

def extinction_factor(keVB, Z, hkl, p, th0=0.): 
    tB = BraggAngle(keVB, d_hkl)
    A = 2 * extinction_parameter(keVB, Z, hkl, p, th0=0.)
    cos2theta = cos(2.*tB)
    arg0 = 2 * A
    arg1 = arg0 * cos2theta
    numerator = 2 * sumoddbess(arg0) + cos2theta * sumoddbess(arg1)
    denominator =  A  *  ( 1. + cos2theta**2. )
    return numerator/denominator

#############################################################

def mat_fac(Z, hkl):
    #return ((sf(Z, hkl) * rem_angstrom / volume )**2. * d_hkl**3 *2)*1.e8
    return (( sf * rem_angstrom / volume )**2. * d_hkl**3 *8)*1.e8


def angular_factor(keV, keVB):
    tetaB = keV2Bragg(keVB)
    teta =  keV2Bragg(keV)
    return  numpy.sin(teta)**3  *  ( 1. + (numpy.cos(2.*tetaB) )**2. ) / ( 2 * numpy.sin(2. * tetaB))

def Q(keV, keVB):
    return mat_fac(Z, hkl) * angular_factor(keV, keVB)

def gaussian(x, xB, p):
    delta=(x-xB) /eta(p)
    exponential=numpy.exp(-delta**2/2)
    normalization = math.sqrt(2.*pi) * eta(p)
    return exponential/normalization


def sigma(keVB, keV, p):
    t = keV2Bragg(keV)
    tB = keV2Bragg(keVB)
    weight=gaussian(t, tB, p)
    ext_fa=extinction_factor(keVB, Z, hkl, p, th0=0.)
    return Q(keV, keVB) * weight * ext_fa
    

def reflectivity(keVB, keV, p):
    T = dim3 
    muT= mu(Z,keV) * T
    sigmaT = abs(sigma(keVB, keV, p) * T)
    r = 0.5 * (1 - numpy.exp(-2 * sigmaT)) * Numeric.exp(-muT)
    return r


#def ref(Z,keV):
#    T = dim3
#    muT= mu(Z,keV) * T
#    r = Numeric.exp(-muT)
#    return r

def residuals(p, y, keV):
    keVB=EB
    err = y-reflectivity(keVB, keV, p) 
    return err


filename=('reflectivity_dx.dat')
data = scipy.io.array_import.read_array(filename)

y = data[:,1]
x = data[:,0]


plsq = leastsq(residuals, p_start, args=(y, x), maxfev=2000)

plot(x,y, x, reflectivity(EB, x , plsq[0]), x, reflectivity(EB, x , p_start) )

#plot (x,y, x, reflectivity(EB, x , p_start))

print plsq[0][0]
print plsq[0][1]/10000

show()

