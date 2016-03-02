"""Contains all the needed physical laws and mathematical relations.
Divided in sections:
- vectors
- math
- real physics

"""
from math import *
from Booklet import *
from scipy.special import jn, erf
from numpy import array
from random import uniform
import sys
from myvec import *

def eta2fwhm(eta):
    return sqrt(8.*math.log(2.))*eta

def fwhm2eta(fwhm):
    return fwhm/sqrt(8.*math.log(2.))

################################################################################
# Math
#
def int_erf(a, r):
    '''Returns the integral of an error function with integration extremes (r+a), (r-a)'''
    if r=='+inf': return 4. * a
    S,D = r+a, r-a
    erf_part = S*erf(S) - D*erf(D)
    exp_part = 1./(math.sqrt(math.pi)) * ( math.exp(-S**2) - math.exp(-D**2) )
    return 2. * ( erf_part + exp_part )
#
def convolution(a, A, B, C, x):
    """Result of a convolution between a rect function and a Gaussian function.
    At the moment A and C are unusefull parameters (amplitude of the rect and the gaussian).
    """
    return (C*A)/B * (math.sqrt(math.pi)/2.) * ( erf(B*(x+a)) - erf(B*(x-a)) )
#
def Rfinder(a, fraction=0.5):
    '''Calculates the value of r that results in a fraction (default 0.5) of the maximum of the
    function.
    At the moment the function is hardcoded and is int_erf (see above)'''
    rguess = a
    halflim = int_erf(a, '+inf')*fraction
    diff = int_erf(a, rguess)/halflim - 1.
    if diff > 0.:
        lower, upper = 0., rguess
    else:
        lower, upper = rguess, rguess*1.e4
    while not abs(diff) < 1.e-6:
        if diff > 0.:
            upper = rguess
        else:
            lower = rguess
        rguess = (lower + upper) /2
        diff = int_erf(a, rguess)/halflim - 1.
    return rguess, diff
#

################################################################################
# Physics
#
def BraggAngle(keV, d_hkl, order=1):
    """Returns the scattering angle from the Bragg law for a given order"""
    return asin(order * hc / (2*d_hkl*keV) )

def BraggEnergy(angle, d_hkl, order=1):
    """Returns the Energy associated with the Bragg angle for a given order"""
    return order * hc / ( 2*d_hkl*sin(angle) )

def angular_factor(tB):
    """Returns the angular dependency of the reflectivity given by Zachariasen."""
    return sin(tB)**2  *  ( 1. + (cos(2.*tB) )**2. ) / cos(tB)

def R_hkl(Z, keV, th0, t_micro, hkl):
    """Returns the integrated reflecting power
    th0 is the incidence angle and t_micro is the thickness of the crystallites
    in micron."""
    # ANGULAR FACTOR
    tB = BraggAngle(keV, d_hkl(Z, hkl))
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
    """Returns the probability per unit of length of a photon of given
    energy to be scattered. (Measure unit cm^-1)"""
    tB=BraggAngle(keV, d_hkl(Z, hkl))
    Q=mat_fac(Z, hkl)*angular_factor(tB)
    DeltaTheta = tB - th0
    if distribution=="hat": weight=hat(DeltaTheta, eta2fwhm(eta))
    elif distribution=="gaussian": weight=gaussian(DeltaTheta, eta)
    else: return 0.
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
    return rem_angstrom * hc * sf(Z, hkl) / ( volume[Z] * cos(th0) )

def extinction_lenght(keV, Z, hkl, th0=0.):
    """Material and energy dependent constant useful for primary extinction calculation. Unit is Angstrom"""
    return keV/extinction_constant(Z, hkl, th0=0.)

def extinction_parameter(keV, Z, hkl, microthick, th0=0.):
    """Returns the adimensional parameter A_0 related to the secondary absorption
    [Zach, p. 169]
    """
    microthickA=microthick*1.e4 # from micron to Angstrom
    return microthickA/extinction_lenght(keV, Z, hkl, th0)

def extinction_factor(keV, Z, hkl, microthick, th0=0.):
    """Returns secondary extinction factor for Laue diffraction
    [Zach, p. 169]
    """
    # Workaround for no secondary extinction
    if microthick==0.: return 1.

    # Start here...
    tB = BraggAngle(keV, d_hkl(Z, hkl))
    A = extinction_parameter(keV, Z, hkl, microthick, th0)
    cos2theta = cos(2.*tB)

    # Useful abbreviations
    arg0 = 2. * A
    arg1 = arg0 * cos2theta

    numerator = sumoddbess(arg0) + cos2theta * sumoddbess(arg1)
    denominator = A  *  ( 1. + cos2theta**2. )
    return numerator/denominator

def reflectivity(keV, Z, hkl, microthick, th0, eta, T, Tamorph=0.):
    """Returns the fraction of photons diffracted in that direction"""
    extinction_factor_ = extinction_factor(keV, Z, hkl, microthick, th0)
    T=T/cos(th0)
    muT=mu(Z, keV) * T
    sigmaT=sigma(Z, keV, th0, eta, hkl) * (T-Tamorph)* extinction_factor_
    return sinh(sigmaT) * exp(-muT -sigmaT)

def reflectivity_layer_no_abs(keV, Z, hkl, T=False):
    """Perfect crystal reflectivity according to Erola. T is the laye thickness!"""
    tB=BraggAngle(keV, d_hkl(Z, hkl))
    Q=mat_fac(Z, hkl)*angular_factor(tB) # [1/cm]
    if T is False:
        TA = extinction_lenght(keV, Z, hkl)
        T=TA*1.e-8 # from Angstrom to cm
        A=1.
    else:
        TA=T*1.e8 # from cm to Angstrom
        A=TA/extinction_lenght(keV, Z, hkl)
    P_kin=Q*T*exp(-mu(Z, keV)*T)
    print P_kin, Q, exp(-mu(Z, keV)*T), T
    return P_kin*tanh(A)/A

def wavevector(keV, kver=(0.,0.,1.)):
    """Given the energy and the versor returns the wavevector"""
    return 2.*pi*(keV/hc)*array(kver)

def best_thickness(keV, Z, hkl, eta):
    """Given the energy and the versor returns the wavevector, returns the best
    thickness."""
    th0 = BraggAngle(keV, d_hkl(Z, hkl))
    sigma_ = sigma(Z, keV, th0, eta, hkl)
    mu_    = mu(Z,keV)
    T_best = log(1.+2.*sigma_/mu_)/2./sigma_
    return T_best

def free_path(coeff_int_lin, d_max='+inf'):
    '''returns path to the next interaction
    [Zaidi]'''
    if d_max=='+inf':
        return -log(uniform(0,1))/coeff_int_lin
    else:
        # Forces free_path to be lower equal to d_max
        return -log(uniform(0,1) * (1.-exp(-coeff_int_lin*d_max)) )/coeff_int_lin

def DWCu300(hkl):
    """Debye Waller factor. Assumed copper Debye T, 320 K"""
    sys.stderr.write("Implemented for Copper and room temperature\n")
    debyefunction300=0.750013710877
    t=343./300.
    M29=63.546*amu # [keV]
    T_part = debyefunction300/t+0.25
    constant_part=6.*(hc**2)/(M29 * (kB/1000.) * 343.) # [A**2], k/1000 [eV -> keV]
    M=constant_part*T_part/(2.*d_hkl(29, hkl))**2
    return math.exp(-2.*M)

# Debye Waller function
# Little bit tricky because of pygsl module...
try:
    from pygsl.sf import debye_1
    def DWfactor(Z, T, hkl):
        """Debye Waller factor for the Z element and hkl Miller indices at temperature T"""
        t=TDebye[Z]/T
        Debye_function=debye_1(t)[0]
        T_part = Debye_function/t+0.25
        constant_part=6.*(hc**2)/(atomic_mass[Z]* amu * (kB/1000.) * TDebye[Z]) # [A**2], k/1000 [eV -> keV]
        M=constant_part*T_part/(2.*d_hkl(Z, hkl))**2
        return math.exp(-2.*M)
except:
    def DWfactor(Z, T, hkl):
        """Debye Waller factor for the Z element and hkl Miller indices at temperature T"""
        sys.stderr.write("Debye Waller factor set to 1.\n")
        return 1.

if __name__=='__main__':
    pass
    print reflectivity_layer_no_abs(200., 14, (1,1,1))
