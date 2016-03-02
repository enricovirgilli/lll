# -*- coding: iso-8859-15 -*-
import os, Physics
math=Physics
Booklet=Physics
array=Physics.array
sys=Physics.sys
pi=math.pi
from numpy import arange

class Material(object):
    """Class defining material properties useful for crystallography."""
    def __init__(self,
                 Z=29, # atomic number
                 fwhm=1., # mosaic spread [arcmin]
                 microthick=0., # microblock thickness [�]
                 ):
        self.Z, self.fwhm, self.microthick = Z, fwhm, microthick

    def __str__(self):
        """Class representation"""
        return "%s, fwhm=%s arcmin, microblock thickens=%s" % (Booklet.name[self.Z], self.fwhm, self.microthick)

    # fwhm property
    def __set_fwhm(self, fwhm):
        self.__fwhm = fwhm # arcmin
        self.eta = Booklet.arcmin2rad(Physics.fwhm2eta(fwhm)) # radians

    def __get_fwhm(self):
        return self.__fwhm

    def __del_fwhm(self):
        pass
    fwhm=property(__get_fwhm, __set_fwhm, __del_fwhm, "FWHM [arcmin]")

    # Z property
    def __set_Z(self, Z):
        self.__Z = Z
        self.density = Booklet.density[Z]
        self.lattice = Booklet.lattice[Z]
        self.volume = Booklet.volume[Z]
        self.TDebye = Booklet.TDebye[Z]
        self.atomic_mass = Booklet.atomic_mass[Z]

    def __get_Z(self):
        return self.__Z

    def __del_Z(self):
        pass

    Z=property(__get_Z, __set_Z, __del_Z, "Atomic number")

    def mu(self, keV):
        """Wrapper function"""
        return Booklet.mu(self.Z, keV)

class Xtal(Material):
    """Class containing info about a crystal."""
    def __init__(self,
                 r = (0.,0.,0.),
                 g = (0.,0.,0.),
                 dim=(0.,0.,0.),
                 Z=29,
                 hkl=(0,0,0),
                 fwhm=1., # [arcmin]
                 microthick=0.): # microblock thickness [�]
        # Init of parent class
        Material.__init__(self, Z, fwhm, microthick)
        # class specific attributes
        self.r = r
        self.g = g
        self.dim = dim
        self.cross_section=dim[0]*dim[1] # [cm2]
        # diffraction attributes and lists
        # Miller indices
        self.hkl = hkl
        self.hkls = [tuple([x*order for x in hkl]) for order in range(0,20)]
        self.hkls = tuple(self.hkls)
        # Plane distances
        self.d_hkl = Booklet.d_hkl(self.Z, self.hkl) # [Angstrom]
        self.d_hkls=[0.]
        [self.d_hkls.append(self.d_hkl/order) for order in range(1,20)] # [Angstrom]
        self.d_hkls = tuple(self.d_hkls)
        # structure factors
        self.sf = Booklet.sf(Z, hkl) # []
        self.sfs=[0.]
        [self.sfs.append(Booklet.sf(Z, miller)) for miller in self.hkls[1:]]
        self.sfs = tuple(self.sfs)
        # A_0const
        if self.Z=="GaAs":
            tmp = Booklet.rem_angstrom * Booklet.hc * microthick / 32.
        else:
            tmp = Booklet.rem_angstrom * Booklet.hc * microthick / self.Z
        self.A_0const = tmp * self.sf # Parameter related to secondary extinction []
        self.A_0consts = [tmp * sf for sf in self.sfs]
        self.A_0consts = tuple(self.A_0consts)
        # material factors
        self.mat_fac = Booklet.mat_fac(Z, hkl) # Basic value of the Material Factor
        self.mat_facs=[0.]
        [self.mat_facs.append(Booklet.mat_fac(Z, miller)) for miller in self.hkls[1:]]
        self.mat_facs = tuple(self.mat_facs)
        #

    def __str__(self):
        """Class representation"""
        return "xtal at %s, phi %s" % (self.r, self.phi)

    # r property
    def __get_r(self): return self.__r

    def __set_r(self, r):
        self.__r = array(r)
        self.__rho = math.sqrt(self.r[0]**2 + self.r[1]**2)
        self.__phi = math.atan2(self.r[1], self.r[0])
        self.__z = self.r[2]

    def __del_r(self): pass

    r=property(__get_r, __set_r, __del_r, "crystal position")

    # rho property
    def __get_rho(self):
        return self.__rho

    def __set_rho(self, rho):
        self.__rho=rho
        self.__r=(rho*math.cos(self.phi),
                rho*math.sin(self.phi),
                self.r[2])

    def __del_rho(self):
        pass

    rho = property(__get_rho, __set_rho, __del_rho, "crystal distance")

    # phi property
    def __get_phi(self):
        return self.__phi

    def __set_phi(self, phi):
        self.__phi=phi
        self.__r=(self.rho*math.cos(phi),
                self.rho*math.sin(phi),
                self.r[2])

    def __del_phi(self):
        pass

    phi = property(__get_phi, __set_phi, __del_phi, "azimuth")

    # g property
    def __get_g(self):
        return self.__g

    def __set_g(self, g):
        self.__g = array(g)
        self.gnorm = Physics.norm(g)
        self.gnormalized = Physics.normalize(g)
        # angular coordinates (azimuth and colatitude) of the g versor
        self.gtheta = math.acos(self.gnormalized[2])
        self.gphi = Physics.atan2(g[1], g[0])

    def __del_g(self):
        pass

    g = property(__get_g, __set_g, __del_g, "reciprocal lattice vector")

    # z property
    def __get_z(self):
        return self.__z

    def __set_z(self, z):
        self.__z=z
        self.__r=(self.r[0], self.r[1], z)

    def __del_z(self):
        pass
    z=property(__get_z, __set_z, __del_z, "crystal height")

    # Functions
    def keV2Bragg(self, keV, order=1):
        '''Calculates the energy related to the the Bragg angle
        for the material the lens is composed of.'''
        return math.asin(Booklet.hc / (2. * self.d_hkls[order]*keV))

    def Bragg2keV(self, theta, order=1):
        '''Calculates the Bragg angle for the input energy for the material the
        lens is composed of.'''
        return Booklet.hc / (2. * self.d_hkls[order]*math.sin(theta))

    def best_thickness(self, keV):
        """Wrapper for the Physics module function"""
        return Physics.best_thickness(keV, self.Z, self.hkl, self.eta)

    def xtal_info(self, filename=False, degree=False):
        """outputs on stdout or on the file defined by filename informations
        about the crystal.
        The angle can be optionally in degrees.
        ACHTUNG!
        The file is by default overwritten.

        TODO:
        output also:
        - dimensions
        - energy
        - size
        maybe pickle is better
        """
        if filename:
            output = file(filename, 'a').write
        else :
            output=sys.stdout.write
        # energy
        tB = math.pi/2 - self.gtheta
        output("%s\t" % self.Bragg2keV(tB))
        # height
        output("%10g\t" % self.z)
        # rho
        output("%10g\t" % self.rho)
        # phi
        output("%10g\t" % self.gphi)
        # theta
        output("%10g\t" % self.gtheta)
        # thickness
        output("%s\n" % self.dim[2])

    def mc(self):
        pass

    def gettB(self, OFFSET=0.):
        """Calculates Bragg angle for an X-ray source with a certain offset.
        NOTE: Assumes source azimuthal angle equal to zero!!!"""
        if OFFSET==0.:
            return math.pi/2-self.gtheta
        else:
            stheta=math.pi*(1.-OFFSET/60./180.)
            return Physics.anglebetween((math.sin(stheta),0.,math.cos(stheta)), self.gnormalized)-math.pi/2.

    def Erange(self, order=1, OFFSET=0., NUM_OF_SIGMA=4):
        """energy range diffracted by the crystal at a given order"""
        Dthetamax = NUM_OF_SIGMA*self.eta
        tB = self.gettB(OFFSET)
        Emin = self.Bragg2keV(tB + Dthetamax, order)
        Emax = self.Bragg2keV(tB - Dthetamax, order)
        return Emin, Emax

    def keVgen(self, order=1, NUM_OF_SIGMA=4):
        """Generates an energy list"""
        Emin, Emax = self.Erange(order, NUM_OF_SIGMA)
        return arange(int(Emin), int(Emax)+1)

    def extinction_factor(self, keV, theta_0, order=1):
        """Wrapper for Physics.extinction_factor function"""
        if self.microthick==0.:
            return 1.
        else:
            # This part can be speed up by using Xtal facilities!!!
            return Physics.extinction_factor(keV, self.Z, self.hkls[order], self.microthick, theta_0)

    def getDelta(self, tB, OFFSET=0.):
        """Function to get the angle between the most probable planes and the
        planes that scatter photons with a certain Bragg angle.
        Contains the obvious answer for in axis photons and the more complex
        calculatin when photons are off axis"""
        if OFFSET==0.:
            return pi/2-self.gtheta-tB
        else:
            ktheta=pi - Physics.arcmin2rad(OFFSET)
            stk, ctk = math.sin(ktheta), math.cos(ktheta)
            stB = math.sin(tB)
            # All starts from a 2nd order equ
            # A=1.
            Bhalf=(math.cos(self.gphi)*stk)/(stB-ctk)
            C=(stB+ctk)/(stB-ctk)
            Delta = Bhalf**2-C
            # In case there are non solutions (energy too low, return None)
            # if Delta<0.: return None
            sqrtDelta = math.sqrt(Delta)
            Ts= (-Bhalf - sqrtDelta, -Bhalf + sqrtDelta)
            Deltas = [abs(2.*math.atan(T) - self.gtheta) for T in Ts]
            return min(Deltas)

    def sigma(self, keV, OFFSET=0., order=1, distribution="gaussian"):
        """Replacement for Physics.sigma function."""
        tB = self.keV2Bragg(keV, order)
        Q = self.mat_facs[order]*Physics.angular_factor(tB)
        Dtheta = self.getDelta(tB, OFFSET)
        if distribution=="hat":
            weight=Physics.hat(Dtheta, Physics.eta2fwhm(self.eta))
        elif distribution=="gaussian":
            weight=Physics.gaussian(Dtheta, self.eta)
        else: return 0.
        return Q*weight/math.cos(tB)

    def reflectivity(self, keV, OFFSET=0., order=1):
        """Very fast reflectivity calculation based on Xtal attributes."""
        extinction_factor = self.extinction_factor(keV, self.keV2Bragg(keV), order)
        T = self.dim[2]/math.cos(pi/2-self.gtheta)
        muT=self.mu(keV) * T
        sigmaT= self.sigma(keV, OFFSET, order) * T * extinction_factor
        return math.sinh(sigmaT) * math.exp(-muT -sigmaT)

    def reflectivity_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4, MAX_ORDER=5):
        """Calculates the reflectivity of a Xtal for the energies given by self.keVgen"""
        tB=self.gettB(OFFSET)
        reflectivity=array([(e,self.reflectivity(e, OFFSET, order)) for order in range(1, MAX_ORDER+1) for e in self.keVgen(order)])
        reflectivity.transpose()
        return reflectivity

    def dump_reflectivity(self, OFFSET=0., NUM_OF_SIGMA=4, MAX_ORDER=5):
        print "\n".join(["\t".join(map(str, x)) for x in self.reflectivity_at_keVs(OFFSET, NUM_OF_SIGMA, MAX_ORDER)])

if __name__=='__main__':
    gaas=Xtal(Z="GaAs", dim=[1,1,0.], hkl=(1,1,1), fwhm=1.)
    germ=Xtal(Z=32, dim=[1,1,0.], hkl=(1,1,1), fwhm=2./3)
    e = 100.
    print e, gaas.best_thickness(e), germ.best_thickness(e)