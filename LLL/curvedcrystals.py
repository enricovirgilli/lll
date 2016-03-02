# -*- coding: iso-8859-15 -*-
#!/usr/bin/env python
from LLL import Lenses
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

#Z='GaAs'                       # atomic number
Z = 14
hkl=(1,1,1)                    # miller indices
#hkl=(2,2,0)# miller indices
#hkl=(2,2,0)


external_radius = 40.  # meters
internal_radius = 2.56 * external_radius  # meters

curvature_s = 1 / internal_radius

fwhm = 1. #in arcmin.

dim1=3.0    #
dim2=1.0    # dimensions of the cristal [cm]
dim3=0.2     # thickness

EB = 200.0   # peak energy of diffraction

density = Booklet.density[Z]

lattice = Booklet.lattice[Z]
volume = Booklet.volume[Z]
TDebye = Booklet.TDebye[Z]
atomic_mass = Booklet.atomic_mass[Z]
sf = Booklet.sf(Z, hkl)
#v_by_f = 0
#v_by_f = volume/abs(sf)
#print v_by_f

def d_hkl(Z, hkl):
    sum_square = sum(i**2 for i in hkl)
    return lattice / math.sqrt(sum_square)

d_hkl=d_hkl(Z, hkl)

def fwhm2eta(x):
    return x/math.sqrt(8.*math.log(2.))

def mu(Z,keV):
    return Booklet.mu(Z, keV)


def keV2Bragg(keV):
    return math.asin(hc / (2. * d_hkl * keV))  #is given in radians

def Bragg2keV(angle):
    return hc / (2. * d_hkl *math.sin(angle))

tB = keV2Bragg(EB)
alpha=2*math.atan((dim1/2)/(external_radius * 100))    
Emin = Bragg2keV(tB + alpha/2)   # min and max energy diffracted by this curved
Emax = Bragg2keV(tB - alpha/2)   # crystal due to the dim[1] and to curvature radius


def lambda0(keVB, keV):
    tB = keV2Bragg(keVB)
    up = math.pi  * math.cos(keV2Bragg(keVB))
    down =  (hc/keV) * rem_angstrom * (1 + abs(math.cos(2*tB)))
    v_by_f = volume / abs(sf)
    r = v_by_f * (up / down)
    return r


def function(x):
    muT = mu(Z,x)
    mu_a = muT/100000000
    cp = curvature_s / 10000000000
    A=math.pi**2 *d_hkl /lambda0(x,x)**2
    B=mu_a/(cp * math.cos(keV2Bragg(x)))
    exponent = math.exp(A/cp) - 1
    f = (2*A) / (B*exponent)
    down =  2 * cp**2 * B * (1 - math.exp(A/cp))
    g = -A/down
    #return f/cp *1e-7
    return g *1e-7

#alist=arange(50, 900, 1)
#for i, e in enumerate(alist):
#    print e, function(e)


def peak_refl_curved(Z, keVB, keV, t):
    tB=keV2Bragg(keVB)
    #omega = 2 * math.atan( (t/2) / (internal_radius * 100)) # Internal Radius

    #T = t * 100000000.0
    #cp =  omega_s(keV) /(t*10000000.0)# (tbest(keVB, keV)*  10000000.0)
    cp = curvature_s/10000000000.0
    expo = (math.pi**2 * d_hkl)/(cp*lambda0(keVB,keV)**2)
    #diffraction = ( 1 - ( Numeric.exp(-(math.pi**2 * d_hkl) / (cp * lambda0(keVB,keV)**2) )) )
    diffraction = ( 1 - math.exp(-expo))
    absorption = Numeric.exp(- (mu(Z,keV)) * (t) / cos(tB))#(tbest(keVB, keV)* 10000000.0) / cos(tB))
    return diffraction * absorption #* omega/alpha

#print peak_refl_curved(Z, 200, 200, 0)

def tbest_T(keVB, keV):
    #omega = 2 * math.atan( (t/2) / (internal_radius * 100)) # Internal Radius
    omega = (fwhm / 60) * (math.pi / 180) # conversion from arcmin to radians
    muT = mu(Z,keV)
    tB=keV2Bragg(keVB)
    mu_angstrom = muT/100000000
    factor = ((math.pi**2) * d_hkl * math.cos(tB)) / (lambda0(keVB,keV)**2 * mu_angstrom * omega) 
    logarithm = log ( 1 + factor)
    return ((omega * logarithm)/((math.pi**2) * d_hkl / lambda0(keVB,keV)**2))/10000000

##############################
def omega_s(keV):
    muT = mu(Z,keV)
    mu_angstrom = muT/100000000
    t_B = keV2Bragg(keV)
    M_factor = ((math.pi**2) * d_hkl) / (lambda0(keV,keV)**2)
    Om = ((math.cos(t_B)/mu_angstrom) * M_factor) / ((math.exp(M_factor/curvature_s)) - 1 )
    return Om

def tbest(keVB, keV):
    #omega = 2 * math.atan( (t/2) / (internal_radius * 100)) # Internal Radius
    #omega = (fwhm / 60) * (math.pi / 180) # conversion from arcmin to radians
    omega = omega_s(keV)
    muT = mu(Z,keV)
    tB=keV2Bragg(keVB)
    mu_angstrom = muT/100000000
    factor = ((math.pi**2) * d_hkl * math.cos(tB)) / (lambda0(keVB,keV)**2 * mu_angstrom * omega) 
    logarithm = log ( 1 + factor)
    return ((omega * logarithm)/((math.pi**2) * d_hkl / lambda0(keVB,keV)**2))/10000000.0 #tBest in mm

def tbest_1(keVB,keV):
    muT = mu(Z,keV)
    tB=keV2Bragg(keVB)
    ab = math.cos(tB) /  muT
    #ab = omega_s(keV) / curvature_s
    return ab

def M_by_N(kevB, keV):
    muT = mu(Z,keV)
    mu_angstrom = muT/100000000
    t_B = keV2Bragg(keV)
    M_factor = ((math.pi**2) * d_hkl) / (lambda0(keV,keV)**2)
    N_factor = (mu_angstrom * omega_s(keV))/math.cos(t_B)
    return M_factor/N_factor

elist=arange(50, 500, 1)
    
#C=[]
#for i, e in enumerate(elist):
    #C.append(M_by_N(e, e))
    #print e, tbest_1(e,e)*10

#plot (elist,C)
#plt.show()

    

######################################################################
# 1) for Energy vs Peak reflectivity fixing the value of thickness
#    which is defined by required quasi-mosaic fwhm
######################################################################
elist=arange(50, 600, 1)
for i, e in enumerate(elist):
    print e, peak_refl_curved(Z,e,e, dim3)
##############################################################


###############################################################
# 2) for thickness vs peak reflectivity
###############################################################

def peak_refl_curved_forThickness(Z, keVB, keV, t):
    tB=keV2Bragg(keVB)
    omega = 2 * math.atan( (t/2) / (internal_radius * 1000)) # Internal Radius

    T = t * 10000000.0
    cp = omega/T
    #cp =  omega_s(keV) / (tbest(keVB, keV)*  10000000.0)
    #cp = curvature_s/10000000000.0
    #muT = mu(Z,keV)
    #mu_angstrom = muT / 100000000.0
    diffraction = ( 1 - ( Numeric.exp(-(math.pi**2 * d_hkl) / (cp * lambda0(keVB,keV)**2) )) )
    absorption = Numeric.exp(- (mu(Z,keV)/100000000.0) * T / cos(tB))
    return diffraction * absorption

def peak_diffraction_curved_forThickness(Z, keVB, keV, t):
    tB=keV2Bragg(keVB)
    omega = 2 * math.atan( (t/2) / (internal_radius * 1000)) # Internal Radius

    T = t * 10000000.0
    cp = omega/T
    #cp =  omega_s(keV) / (tbest(keVB, keV)*  10000000.0)
    #cp = curvature_s/10000000000.0
    #muT = mu(Z,keV)
    #mu_angstrom = muT / 100000000.0
    #diffraction = ( 1 - ( Numeric.exp(-(math.pi**2 * d_hkl) / (cp * lambda0(keVB,keV)**2) )) )
    interaction = 1 - (Numeric.exp(- (mu(Z,keV)/100000000.0) * T / cos(tB)))
    return interaction


elist=arange(100, 500, 100)#in keV
tlist = arange(0,50,.2) #in mm
Curvature_list = arange(10, 120, 10) #in 1/m 
aR = []
aT = []
aI = []

#########################################################################
#Plotting Thickness vs Omega in arc seconds...!!!
########################################################################
#aTheta = []
#for i, c in enumerate(Curvature_list):
    #aT.append(tbest(e, e))
    #aTheta = []
    #for i, t in enumerate(tlist):
        #aTheta.append(3600. * 180. * (1 / math.pi) * (t/1000.) /c  )
    #plot (tlist, aTheta)
    #plt.hold(True)
#plt.axis([0, 20, 0, 80])
#plt.xlabel('Thickness in mm')
#plt.ylabel('Omega in arcSec')
#plt.title('Thicness vs omega')
#plt.grid(True)
#plt.savefig('T_vs_omega_in_arcSec.ps')
#plt.show()

#---------------------------------------------------------------------------
    
        
for i, e in enumerate(elist):
    #aT.append(tbest(e, e))
    aR = []
    aI = []
    aM = []
    for i, t in enumerate(tlist):
        #aT.append(t)
        #t = tbest(e, e)
        #print e, peak_refl_curved_forThickness(Z,EB, EB, e)
        aR.append(peak_refl_curved_forThickness(Z, e, e, t))
        aI.append(peak_diffraction_curved_forThickness(Z, e, e, t))
        aM.append(peak_refl_curved_forThickness(Z, e, e, t)*peak_diffraction_curved_forThickness(Z, e, e, t))
        #aR.append(peak_refl_curved_forThickness(Z, e, e, t))

    
    #matplotlib.lines.Line2D(tlist, aR)
    if e == 100:
        plot (tlist, aR, label='100')
        plt.hold(True)
        plot (tlist, aI, label='100')
        plt.hold(True)
        plot (tlist, aM, label='100')
        plt.hold(True)
       
    elif e == 200:
        plot (tlist, aR, label='200')
        plt.hold(True)
        plot (tlist, aI, label='200')
        plt.hold(True)
        plot (tlist, aM, label='200')
        plt.hold(True)
        
    elif e == 300:
        plot (tlist, aR, label='300')
        plt.hold(True)
        plot (tlist, aI, label='300')
        plt.hold(True)
        plot (tlist, aM, label='300')
        plt.hold(True)
       
    elif e == 400:
        plot (tlist, aR, label='400')
        plt.hold(True)
        plot (tlist, aI, label='400')
        plt.hold(True)
        plot (tlist, aM, label='400')
        plt.hold(True)
       
    elif e == 500:
        plot (tlist, aR, label='500keV')
        plt.hold(True)
        plot (tlist, aI, label='500keV')
        plt.hold(True)
        plot (tlist, aM, label='500keV')
        plt.hold(True)
        
        
    #plot (tlist, aR)  #, label='energy')
    #plt.hold(True)

plt.grid(True)
plt.legend()
plt.xlabel('Thickness in mm')
plt.ylabel('Reflectivity')
plt.title('GaAs(333)')
plt.savefig('RwithIvsThick_mult_GaAs333.ps')
plt.show()


    
###############################################################

######################################################################
# 1) for Thickness vs Energy 
#    
######################################################################
#elist=arange(50, 500, 1)
#B=[]
#for i, e in enumerate(elist):
    #B.append(tbest(e, e))
    #print e, tbest_1(e,e)

#plot (elist,B)
#plt.show()
##############################################################

#omega = 2 * math.atan( (dim3/2) / (internal_radius * 100)) 
#print "alpha=", alpha/3.14*180*60*60, "omega=", omega/3.14*180*60*60


#enelist=arange(50, 400, 1)
#for i, e in enumerate(enelist):
#    print e, peak_refl_curved(Z,e,e, dim3)
################################################################

#sys.exit(0)

    
def keVgen(keVB, NUM_OF_SIGMA):
    DeltaEnergy = NUM_OF_SIGMA * (Emax-Emin)
    Einf = Emin - DeltaEnergy/2
    Esup = Emax + DeltaEnergy/2
    s=arange(Einf*200,Esup*200)
    return s/200

def rectfun(xB, x, t):
    T = tbest(xB, x) * 10000000.0
    #omega = 2 * math.atan( (t/2) / (internal_radius * 100))
    cp = omega_s(xB) / T
    muT = mu(Z,x)/100000000.0
    #absorption = Numeric.exp(-(mu(Z,x)/100000000.0 * omega ) / (cp * math.cos(keV2Bragg(x)))  )
    if  Emin < x < Emax : return peak_refl_curved(Z, xB, x, t) #* 1-0.01*abs(x-EB)**15/EB
    else:
        return 0


###############################################################
# for Energy vs best thickness
###############################################################
elow = 50.
eup = 600.
estep=1.
eerange= int(round((eup - elow )/estep))

tlow = 0.1 
tup =  5.0
tstep=0.1
ttrange= int(round((tup - tlow )/tstep))


elist=arange(elow, eup, estep)
tlist=arange(tlow, tup, tstep)

b = [[]*eerange]  # initialize M to have one row with 6 zeros, as first row in matirx

#print 250, peak_refl_curved(Z, 250, 250, 1.45)
a = []

#for i, e in enumerate(elist):
    #print peak_refl_curved(Z, EB, EB, 0)
    #a.append(peak_refl_curved(Z, e, e, 0))
    #a.append(lambda0(EB,EB))
    
#plot (elist,a)
#plt.show()


#for i, e in enumerate(elist):
    #print peak_refl_curved(Z, EB, EB, 0)
    #a.append(peak_refl_curved(Z, e, e, 0))
    
    
#plot (elist,a)
#plt.show()


    #ee=int(round(e-elow))
    #for i, t in enumerate(tlist):     
        #tt=int(round(t*10-tlow))        
        #a.append(tbest(e, e))
     #   a.append(peak_refl_curved(Z, e, e, t))

e = 100.
#rr = peak_refl_curved(Z, e, e, 0)
#print rr


#plot (f,rectarray)
#plt.show()

#################################################################


#def smooth(x,window_len=11,window='hanning'):
#    s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
#    if window == 'flat':
#        w=numpy.ones(window_len,'d')
#    else:
#        w=eval('numpy.'+window+'(window_len)')
#    y=numpy.convolve(w/w.sum(),s,mode='same')
#    return y[window_len:-window_len+1]


#rectarray = array([rectfun(EB, e, dim3) for e in keVgen(EB, 2)])
#f = keVgen(EB, NUM_OF_SIGMA=2)


#reflecurved = smooth(rectarray,window_len=EB**1.09-85,window='hanning')

#plot (f,reflecurved, label='line 1')
#plot (f,rectarray)
#plt.show()


#####################################
#PLOTTING REFLECTIVITY AS A FUNCTION OF (VOLUME BY STRUCTURE FACTOR)
#####################################

def lambda_0(keVB, keV, v_f):
    tB = keV2Bragg(keVB)
    up = math.pi  * math.cos(keV2Bragg(keVB))
    down =  (hc/keV) * rem_angstrom * (1 + abs(math.cos(2*tB)))
    #v_by_f = volume / abs(sf)
    #v_by_f = v_f
    r = v_f * (up / down)
    return r


def peak_refl_curved_test(Z, keVB, keV, lmbd0, t=0):
    tB=keV2Bragg(keVB)
    #omega = 2 * math.atan( (t/2) / (internal_radius * 100)) # Internal Radius

    #T = t * 100000000.0
    #cp =  omega_s(keV) / (tbest(keVB, keV)*  10000000.0)
    cp = curvature_s/10000000000.0
    #muT = mu(Z,keV)
    #mu_angstrom = muT / 100000000.0
    diffraction = ( 1 - ( math.exp(-(math.pi**2 * d_hkl) / (cp * lmbd0**2) )) )
    absorption = Numeric.exp(- (mu(Z,keV)/100000000.0) * (tbest(keVB, keV)* 10000000.0) / cos(tB))
    return diffraction * absorption


#ab = []
#v_by_F = arange(0.01, 6., .01)
#for i, e in enumerate(v_by_F):
    #print peak_refl_curved(Z, EB, EB, 0)
    #lmbd0 = lambda_0(EB, EB, e)
    #ab.append(peak_refl_curved_test(Z, EB, EB, lmbd0, 0))
    #print peak_refl_curved_test(Z, EB, EB, lmbd0, 0)
        
#plot (v_by_F,ab)
#plt.show()
#plt.grid




#####################################
#PLOTTING AN EXPONENTIAL FUNCTION
#####################################

#exp_arr = arange(0, 10, .1)

#a = []
#for i, e in enumerate(exp_arr):
    #print peak_refl_curved(Z, e, e, 0)
    #a.append (1- math.exp(-(1/e)))
    
#exp_reslt = array( math.exp(-e) for e in arange(1, 100, 1))
#plot (exp_arr, a)
#plt.show()

