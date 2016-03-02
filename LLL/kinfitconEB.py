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

fwhm=0.9
EB=99.2
microthick_micron=25   # microthickness [micron]
fondo = 0.
rem_angstrom = 2.81794092e-5 # [electron radius in Angstrom]

Z=29                   # atomic number


hkl=(1,1,1)            # miller indices

D=200    # distance source-cristal [cm]

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

p_start = array([fwhm, microthick, EB, fondo])


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

def extinction_lenght(p, Z, hkl, th0=0.):
    return p[2]/extinction_constant(Z, hkl, th0=0.)

def extinction_parameter(Z, hkl, p, th0=0.):
    #tB = BraggAngle(keV, d_hkl)
    #polarization = 1 + abs(cos(2*tB))
    return p[1] / extinction_lenght(p, Z, hkl, th0)

def extinction_factor(Z, hkl, p, th0=0.): 
    tB = BraggAngle(p[2], d_hkl)
    A = 2 * extinction_parameter(Z, hkl, p, th0=0.)
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


def angular_factor(keV, p):
    tetaB = keV2Bragg(p[2])
    teta =  keV2Bragg(keV)
    return  numpy.sin(teta)**3  *  ( 1. + (numpy.cos(2.*tetaB) )**2. ) / ( 2 * numpy.sin(2. * tetaB))
    

def Q(keV, p):
    return mat_fac(Z, hkl) * angular_factor(keV, p)

def gaussian(x, xB, p):
    delta=(x-xB) /eta(p)
    exponential=numpy.exp(-delta**2/2)
    normalization = math.sqrt(2.*pi) * eta(p)
    return exponential/normalization


def sigma(keV, p):
    t = keV2Bragg(keV)
    tB = keV2Bragg(p[2])
    weight=gaussian(t, tB, p)
    ext_fa=extinction_factor(Z, hkl, p, th0=0.)
    return Q(keV, p) * weight * ext_fa
    

def reflectivity(keV, p):
    T = dim3 
    muT= mu(Z,keV) * T
    sigmaT = abs(sigma(keV, p) * T)
    r = p[3] + 0.5 * (1 - numpy.exp(-2 * sigmaT)) * Numeric.exp(-muT)
    return r

def residuals(p, y, keV):
    err = y-reflectivity(keV, p) 
    return err


filename=('reflectivity_dx.dat')
data = scipy.io.array_import.read_array(filename)

y = data[:,1]
x = data[:,0]

#plot(x, y, 'o')

plsq = leastsq(residuals, p_start, args=(y, x), maxfev=2000)
microth = plsq[0][1]/10000


####################################################################################################
# Divergence effect
####################################################################################################
EB = 99
back= 0.003

c_start = array([fwhm, microthick])

def  eta(c):
    return  ( c[0] * math.pi / (180 * 60) ) / (math.sqrt(8.*math.log(2.)) )

def extinction_lenght(EB, Z, hkl, th0=0.):
    return EB/extinction_constant(Z, hkl, th0=0.)

def extinction_parameter(Z, hkl, c, th0=0.):
    return c[1] / extinction_lenght(EB, Z, hkl, th0)

def extinction_factor(Z, hkl, c, th0=0.): 
    tB = BraggAngle(EB, d_hkl)
    A = 2 * extinction_parameter(Z, hkl, c, th0=0.)
    cos2theta = cos(2.*tB)
    arg0 = 2 * A
    arg1 = arg0 * cos2theta
    numerator = 2 * sumoddbess(arg0) + cos2theta * sumoddbess(arg1)
    denominator =  A  *  ( 1. + cos2theta**2. )
    return numerator/denominator

def angular_factor(keV, EB):
    tetaB = keV2Bragg(EB)
    teta =  keV2Bragg(keV)
    return  numpy.sin(teta)**3  *  ( 1. + (numpy.cos(2.*tetaB) )**2. ) / ( 2 * numpy.sin(2. * tetaB))
    #return sin(tetaB)**2  *  ( 1. + (cos(2.*tetaB) )**2. ) / cos(tetaB)

def Q(keV, EB):
    return mat_fac(Z, hkl) * angular_factor(keV, EB)

def gaussian(x, xB, c):
    delta=(x-xB) /eta(c)
    exponential=numpy.exp(-delta**2/2)
    normalization = math.sqrt(2.*pi) * eta(c)
    return exponential/normalization

def sigma(keV, EB, c):
    t = keV2Bragg(keV)
    tB = keV2Bragg(EB)
    weight=gaussian(t, tB, c)
    ext_fa=extinction_factor(Z, hkl, c, th0=0.)
    return Q(keV, EB) * weight * ext_fa
    
def reflectivity(keV, c):
    T = dim3 
    muT= mu(Z,keV) * T
    sigmaT = abs(sigma(keV, EB, c) * T)
    r = 0.5 * (1 - numpy.exp(-2 * sigmaT)) * Numeric.exp(-muT)
    return r

def f(c, EB, keV):
    T = dim3
    muT= mu(Z,keV) * T
    sigmaT = abs(sigma(keV, EB, c) * T)
    ratio=hc/(2*d_hkl*EB)
    radice=1/math.sqrt(1-ratio**2)
    g =  0.5 * (1 - numpy.exp(-2 * sigmaT)) * Numeric.exp(-muT) * hc/(2*d_hkl) * radice / EB**2
    return g


def integrefle(c, EB, keV):
    n=30
    tBmin= keV2Bragg(EB)-dim1/(D)
    tBmax= keV2Bragg(EB)+dim1/(D)
    EBmin=Bragg2keV(tBmax)
    EBmax=Bragg2keV(tBmin)    
    h = 0.5 #float (EBmax - EBmin)/n
    s = 0.5 * ( f(c, EBmin, keV) +  f(c, EBmax, keV))
    for i in range (1,n):
        s=s + f(c, EBmin+i*h, keV) * h
    return s * D / (dim1*dim1)


###################################################################################################################

#legend( ("mu = %.2f" % plsq[0][0], '$t_{0}$=  %.2f' % microth, '$E_{p}$= %.2f' % plsq[0][2],  'bkg= %.4f' % plsq[0][3] ), loc=0 )
#z=x/x - 1

plot(x, integrefle(c_start,EB, x) )
show()
sys.exit(0)


#pylab.subplot(212)
#plot(x, y-reflectivity(x , plsq[0]), 'o' )
#plot(x,z, '--')

print "Mu: %.2f" % plsq[0][0]
print "t0: %.2f" % microth
print "Ep: %.2f" % plsq[0][2]
print "Background: %.4f" % plsq[0][3]

