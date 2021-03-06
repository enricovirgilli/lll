# -*- coding: iso-8859-15 -*-
"""Contains the physical constants needed for the calculations"""
import math, cmath, glob, sys
from os import environ
import spibg
from pylab import load
from xraylib import *
XRayInit()

# COSTANTS
# from the Particle Physics Booklet
# Speed of light
# c=299792458. # [m/s]
# Planck
# h=6.6260755e-34 # [Js]
# hbar=1.05457266e-34 # [Js]

# classical electron radius [fm] and [�] and rest mass [keV]
rem = 2.81794092 # [fm]
rem_angstrom = 2.81794092e-5 # [�]
m_e = 510.999 # [keV]

# Conversion constant
hc=12.39842 # [keV x �]
# Boltzmann constant
kB=8.617385e-5 # [eV/K]
# Atomic mass unit
amu=931494.32 # [keV/c2]


# some useful conversions about angles
# unit is radian
# degree = 0.01745329
# arcmin = degree/60.
# arcsec = arcmin/60.

# NAMES
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
        "GaAs": 300., # Valore medio... calcolato da me
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
             "GaAs" : 70., # Valore medio... calcolato da me
             }

# LINEAR ABSORPTION COEFFICIENT [cm^-1]
# reads data from xcom data
mudata={}

def lnbg_parser():
    """Default BG since 07/06/18"""
    fname="%s/code/python/lib/data/lnbg/lnbg.dat" % environ["HOME"]
    return load(fname)[:,1]/4.5*2 # Takes into account the conversion from cts/(s bin) to cts/(s cm^2 keV)
lnbgdata=lnbg_parser()

def bg(keV):
    """Gets bg for the specified energy"""
    if 100.<=keV<1400.:
        i=(float(keV)-100.)*2 # conversion to bin
        ii=int(i)
        deltai=ii-i
        q=lnbgdata[i]
        m=(lnbgdata[i+1]-q)
        return m*deltai+q # linear interpolation. There is better around...
    else:
        sys.stderr.write("Warning: background information unavailable!\n")
        return 1.

def arcmin2rad(theta):
    """arcminute to radian conversion"""
    return theta*math.pi/60./180.

def rad2arcmin(theta):
    """radian to arcminute conversion"""
    return theta*60.*180./math.pi

def mu(Z, keV):
    """Linear absorption coefficient the Z element at energy keV.
    Data are automatically loaded if needed using the exception control.

    TODO: add interpolation such that non integer energy values can
    be better evaluated.
    """
    try:
        keV_int=int(keV)
        keV_float=keV-int(keV)
        y0=mudata[Z][keV_int]
        y1=mudata[Z][keV_int+1]
        y=y0+(y1-y0)*keV_float
        return y
    except KeyError:
        DATADIR = "%s/code/python/lib/Xtal/data/cross_sections" % environ["HOME"]
        mudata[Z]=[float(x.split()[2])*density[Z] for x in file("%s/%s.dat" % (DATADIR, Z))]
        return mudata[Z][int(keV)]

def absorb(Z, keV, thickness):
    """Photon fraction absorbed by a layer of the Z element of thickness T [cm]."""
    return math.exp(-mu(Z, keV)*thickness)

def d_hkl(Z, hkl):
    """Lattice plane distance for Z element and hkl Miller indices [Angstrom].
    The structure is assumed to be cubic but works also for HOPG (0, 0, 2X)."""
    sum_square = sum(i**2 for i in hkl)
    return lattice[Z] / math.sqrt(sum_square)

def ff(Z, hkl):
    """Atomic form factor for material Z and Miller indices hkl through xraylib"""
    if Z=="GaAs":
        ff_Ga=FF_Rayl(31, 0.5/d_hkl("GaAs", hkl))
        ff_As=FF_Rayl(32, 0.5/d_hkl("GaAs", hkl))
        return ff_Ga+ff_As*cmath.exp(1.j*cmath.pi/2*sum(hkl))
    return FF_Rayl(Z, 0.5/d_hkl(Z,hkl))

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
    elif Z in (29, 47, 79, "GaAs"):
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

def phi_hkl(Z, E, hkl):
    """Fourier component of the polarizability as in Zachariasen p. 114 but
    taken in absolute value and slightly modified"""
    lambdasquare=(hc/E)**2
    return rem_angstrom*lambdasquare*sf(Z, hkl)/(pi*volume[Z])

def mat_fac(Z, hkl):
    """Frequently used constant for element and hkl Miller indices [cm]."""
    # the final factor is needed for the conversion from 1/Angstrom to 1/cm
    return ((sf(Z, hkl) * rem_angstrom / volume[Z] )**2. * d_hkl(Z, hkl)**3 *2)*1.e8

if __name__=='__main__':
    for keV in range(100, 1000): print keV, bg(keV)