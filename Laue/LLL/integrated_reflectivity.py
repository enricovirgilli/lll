# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python
from LLL import Lenses
import spibg
import numpy
import math, cmath, glob, sys
from pylab import load
from xraylib import *
from LLL import __path__ as LLLpath
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
import Booklet
from  Numeric  import zeros

Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi
rem_angstrom = 2.81794092e-5 # [electron radius in Angstrom]
hc=Booklet.hc

eta_arcsec = 20.0
Z= "GaAs" #14                      # atomic number
#Z = 32
#hkl=(1,1,1)                    # miller indices
hkl=(2,2,0)# miller indices
#hkl=(1,1,1)
#keV = 59.2
microthick = 100
d_hkl = Booklet.d_hkl(Z, hkl) # [Angstrom]
thickness = 10.0 #[cm]


eta = (eta_arcsec / 3600) * (math.pi/180)
fwhm = 20 #Physics.eta2fwhm(eta)

name={6  : "HOPG",
      8  : "Oxygen",
      13 : "Alluminum",
      14 : "Silicon",
      29 : "Copper",
      32 : "Germanium",
      42 : "Molybdenum",
      47 : "Silver",
      79 : "Gold",
      82 : "Lead",
      "air" : "Air",
      "Si02": "Quartz",
###########################
      "GaAs": "GaAs",
      28 : "Nikel",
      "InAs" : "InAs",
      "CdTe" : "CdTe",
      73 : "Tantalum",
      74 : "Tungsten",
      }



# SYMBOLS
symbols={6  : "C",
         8  : "O",
         13 : "Al",
         14 : "Si",
         29 : "Cu",
         32 : "Ge",
         42 : "Mo",
         47 : "Ag",
         79 : "Au",
         82 : "Pb",
         ##############################
         "GaAs": "GaAs",
         28 : "Ni", 
         "InAs" : "InAs" ,
         "CdTe" : "CdTe",
         73 : "Ta",
         74 : "W",  
         }


# LATTICE PARAMETER [�]
# http://www.webelements.com
lattice={6  : 6.711, # hexagon side = 2.464
         13 : 4.0495,
         14 : 5.4309,
         26 : 2.8665,
         29 : 3.6149,
         32 : 5.6575,
         42 : 3.147,
         47 : 4.0853,
         73 : 3.3013,
         74 : 3.1652,
         79 : 4.0782,
         "GaAs" : 5.6535, # Filippo says
#############################################
         28 : 3.524,
         "InAs" : 6.0583,
         "CdTe" : 6.482,
         73 : 3.3013,
         74 : 3.1652,          
         }



# CELL VOLUMES [�]^3
# Except for the HOPG (Z=6) the cell is cubic...
volume={6  : 35.1660,
        13 : lattice[13]**3,
        14 : lattice[14]**3,
        26 : lattice[26]**3,
        29 : lattice[29]**3,
        32 : lattice[32]**3,
        42 : lattice[42]**3,
        47 : lattice[47]**3,
        73 : lattice[73]**3,
        74 : lattice[74]**3,
        79 : lattice[79]**3,
        "GaAs" : lattice["GaAs"]**3,
        "SiO2" : 113.01,
###################################################
        28 : lattice[28]**3,
        "InAs" : lattice["InAs"]**3,
        "CdTe" : lattice["CdTe"]**3, 
        }



# DENSITIES (g/cc)
# http://www.webelements.com
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

# DEBYE TEMPERATURES [K]
# http://www.dl.ac.uk/MEIS/periodictable/text/index.htm
TDebye={14: 645.,
        26: 470.,
        29: 343.,
        32: 374.,
        42: 450.,
        47: 225.,
        79: 225.,
        "GaAs": 360., # Valore medio... calcolato da me
############################################
        "InAs": 350.,
        "CdTe": 350.,  # inventato da me, cercarlo in qualche tabella
        28: 450.,
        73: 240.,       
        74: 400.,


        }

# ATOMIC MASSES [amu]
# http://www.webelements.com
atomic_mass={6 : 12.0107,
             13 : 26.981538,
             14 : 28.0855,
             26 : 55.845,
             29 : 63.546,
             32 : 72.64,
             42 : 95.94,
             47 : 107.8682,
             73 : 180.9479,
             74 : 183.84,
             79 : 196.96655,
             "GaAs" : 144.63, # Valore medio... calcolato da me
###########################################################
             "InAs" : 94.87,
             "CdTe" : 100.00,  # inventato da me!!!!!
             28 : 58.6934,
             50 : 144.33,
             }


#######################################################################
def ff(Z, hkl):
    """Atomic form factor for material Z and Miller indices hkl through xraylib"""
    if Z=="GaAs":
        ff_Ga=FF_Rayl(31, 0.5/d_hkl)
        ff_As=FF_Rayl(32, 0.5/d_hkl)
        return ff_Ga+ff_As*cmath.exp(1.j*cmath.pi/2*sum(hkl))
    if Z=="CdTe":
        ff_Cd=FF_Rayl(48, 0.5/d_hkl)
        ff_Te=FF_Rayl(52, 0.5/d_hkl)
        return ff_Cd+ff_Te*cmath.exp(1.j*cmath.pi/2*sum(hkl))
    return FF_Rayl(Z, 0.5/d_hkl)
###########################################################################

def reduced_sf_HOPG(Z, hkl):
    """Returns the reduced structure factor for HOPG.
    Works only for planes like (0 ,0, 2X)"""
    if Z==6: # HOPG only for (0,0, X) planes
        if not (hkl[0], hkl[1])==(0,0):
            sys.stderr.write("Not implemented diffraction by planes different from (0 ,0, 2X) for HOPG")
            sys.exit()
        if hkl[2]%2: return 4.
        else: return 0.

def reduced_sf_copper(Z, hkl):
    """Returns the reduced structure factor for copper structure."""
    tmp_reduced_sf=1.
    for i in range(3):
        tmp_reduced_sf+=cmath.exp(cmath.pi*1.j*(hkl[i%3] + hkl[(i+1)%3]))
    return abs(tmp_reduced_sf)

def reduced_sf_diamond(Z, hkl):
    """Returns the reduced structure factor for diamond structure.
    Based on copper structure factor with the inclusion of the basis."""
    return abs(reduced_sf_copper(Z, hkl)*(1.+cmath.exp(1.j*cmath.pi/2*sum(hkl))))

def reduced_sf(Z, hkl):
    """Returns the reduced structure factor"""
    if Z==6: reduced_sf_HOPG(Z, hkl)
    elif Z in (29, 47, 79, "GaAs", "CdTe"):
        return reduced_sf_copper(Z, hkl)
    elif Z in (14, 32):
        return reduced_sf_diamond(Z, hkl)
    elif Z in (42,):
        if sum(hkl)%2==0: return 2.
        else: return 0.
    else: return 0.

def sf(Z, hkl):
    """Structure factor for Z element and hkl Miller indices."""
    return ff(Z, hkl)*reduced_sf(Z, hkl)
    

def mat_fac(Z, hkl):
    return ((abs(sf(Z, hkl)) * rem_angstrom / volume[Z] )**2. * d_hkl**3 *2)*1.e8

def BraggAngle(keV, d_hkl, order=1):
    """Returns the scattering angle from the Bragg law for a given order"""
    return math.asin(order * hc / (2*d_hkl*keV) )

def BraggEnergy(angle, d_hkl, order=1):
    """Returns the Energy associated with the Bragg angle for a given order"""
    return order * hc / ( 2*d_hkl*sin(angle) )

def angular_factor(tB):
    """Returns the angular dependency of the reflectivity given by Zachariasen."""
    return sin(tB)**2  *  ( 1. + (cos(2.*tB) )**2. ) / cos(tB)

def R_hkl(Z, keV, th0, t_micro, hkl):
    tB = BraggAngle(keV, d_hkl)
    ang_fac = angular_factor(tB)

    # MATERIAL FACTOR times angular factor [1/cm]
    Q=mat_fac(Z, hkl)*ang_fac

    # Now I mix all this parts togheter and divide for cos(th0)
    # return Q*weight/cos(th0)
    t_cm=t_micro/1.e4 # Converts from micron to cm
    return Q* t_cm/cos(th0)

def hat(x, fwhm=1., x0=0.):
    """Hat distribution function normalized to one and centered in 0."""
    if abs(2.*(x-x0)) < fwhm: return 1./fwhm
    else: return 0.

def gaussian(x, eta, x0=0.):
    """Normalized Gaussian distribution function centered in 0. by default."""
    delta=((x-x0)/eta)
    exponential=exp(-delta**2/2.)
    normalization = sqrt(2.*pi) * eta
    return exponential/normalization

def sigma(Z, keV, th0, eta, hkl, distribution="gaussian"):
    th0=0
    tB=BraggAngle(keV, d_hkl)
    Q=mat_fac(Z, hkl)*angular_factor(tB)
    DeltaTheta = tB - th0
    weight=gaussian(DeltaTheta, eta)
    return Q*weight/cos(th0)

def sumoddbess(x, lower_limit=1.e-10):
    """Returns the series sum(i=1, inf) jn(2i+1, x).
    Stops the series when the sum of the last two encountered terms become lower
    than lower_limit.
    [Zach, p. 169]
    """
    tmp, i =0., 0
    # Control variables
    ctrl1, ctrl2 = 1., 1.
    while (abs(ctrl1)+abs(ctrl2) >= lower_limit):
        ctrl2=ctrl1
        ctrl1 = jn(2*i+1, x)
        tmp+=ctrl1
        i+=1
    return tmp

def extinction_constant(Z, hkl, th0=0.):
    """Material dependent constant useful for primary extinction calculation. Unit is keV/A"""
    return rem_angstrom * hc * abs(sf(Z, hkl)) / ( volume[Z] * cos(th0) )

def extinction_lenght(keV, Z, hkl, th0=0.):
    """Material and energy dependent constant useful for primary extinction calculation. Unit is Angstrom"""
    return keV/extinction_constant(Z, hkl, th0=0.)

def extinction_parameter(keV, Z, hkl, microthick, th0=0.):
    
    tB = BraggAngle(keV, d_hkl)
    polarization = 1 + abs(cos(2*tB))
    microthickA=microthick#*1.e4 # from micron to Angstrom
    return microthickA * polarization /extinction_lenght(keV, Z, hkl, th0)


def extinction_factor(keV, Z, hkl, microthick, th0=0.):
    if microthick==0.: return 1.

    tB = BraggAngle(keV, d_hkl)
    A = extinction_parameter(keV, Z, hkl, microthick, th0)
    cos2theta = cos(2.*tB)

    arg0 = 2. * A
    arg1 = arg0 * cos2theta

    numerator = sumoddbess(arg0) + cos2theta * sumoddbess(arg1)
    denominator =  A  *  ( 1. + cos2theta**2. )
    return numerator/denominator



def int_refl(Z, keV, th0, eta, hkl):
    extinction = extinction_factor(keV, Z, hkl, microthick, th0)
    tB=BraggAngle(keV, d_hkl)
    Q=mat_fac(Z, hkl)*angular_factor(tB)
    A = Q * thickness * extinction / cos(tB)
    muTgamma=Booklet.mu(Z, keV) * thickness/cos(tB)
    return A# *  exp(-muTgamma)

kevRange = arange(50, 600, 1)
for i, energy_ in enumerate(kevRange):
    print int_refl(Z, energy_, 0, eta, hkl)

