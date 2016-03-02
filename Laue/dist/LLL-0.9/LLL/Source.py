"""Source definition"""
import Booklet, Physics, Xtal
import Physics
Booklet=Physics
math=Physics
pi=math.pi
from random import uniform
from numpy import zeros, array

EPSILON=1.e-10

class Photon(object):
    """Class containing info about a Photon.

    You can get information about the photon and much more.
    At init time you have to declare the Photon with three optional
    parameters:
    - i the number of photons (intensity);
    - r it's position;
    - k it's wavevector.

    After k is defined also relateed quantities derived from it are
    calculated, namely:
    - lam:  the wavelength
    - e:    the energy
    - knormalized: the unit vector of k
    - knorm: the norm of k

    These parameters can be changed also later.
    """
    def __init__(self,
                 i = 1.,
                 r = zeros(3),
                 k = zeros(3),
                 e = 100.):
        self.r = array(r) # position of the photon
        self.k = array(k) # wavevector
        self.i=i # intensity
        # set up some things that vary together with k
        self.k=k
        self.e=e

    # k property
    def __get_k(self):
        return self.__k

    def __set_k(self, k):
        self.__k = array(k)
        self.knorm=Physics.norm(k)
        self.knormalized = Physics.normalize(k)
        # Calculates wavelength and energy of the photon
        if self.knorm==0.:
            self.__e, self.__lam=0., '+inf'
        else:
            self.__lam=2.*pi/self.knorm
            self.__e=Booklet.hc/self.__lam

    def __del_k(self):
        pass

    k=property(__get_k, __set_k, __del_k, "wavevector")

    # e property
    def __get_e(self):
        return self.__e

    def __set_e(self, e):
        '''Changes the energy of the photon and related quantities'''
        self.__e=e
        self.__lam=Booklet.hc/e
        self.knorm=2.*pi/self.__lam
        self.__k=self.knorm*self.knormalized

    def __del_e(self):
        pass

    e = property(__get_e, __set_e, __del_e, "energy")

    # lam property
    def __get_lam(self):
        return self.__lam

    def __set_lam(self, lam):
        '''Changes the energy of the photon and related quantities'''
        self.__lam=lam
        if lam==0.: self.__e="+inf"
        else: self.__e=Booklet.hc/lam
        self.knorm=2.*pi/self.__lam
        self.__k=self.knorm*self.knormalized

    def __del_lam(self):
        pass

    lam = property(__get_lam, __set_lam, __del_lam, "wavelength")

    # Functions to provide rotation and translation to the photon.
    def travelz(self, z=0):
        """Move the photon until it reaches the z coordinate reaches the desired value (default z=0).
        Note that this transport function can make the photon go back!"""
        deltaz=z - self.r[2]
        self.r=self.r + deltaz/self.k[2]*self.k

    def travel(self, path):
        '''Makes the photon go forth or back of a to desired path'''
        self.r= self.r + path*self.knormalized

    def deviate(self, theta, phi):
        '''Changes the photon wavevector after a scattering [Penelope]'''
        ctheta, stheta = math.cos(theta), math.sin(theta)
        cphi, sphi = math.cos(phi), math.sin(phi)
        if abs(self.knormalized[2])==1.:
            self.k=(self.knormalized[2]*stheta*cphi,
                    self.knormalized[2]*stheta*sphi,
                    self.knormalized[2]*ctheta)
        else:
            rad=math.sqrt(1.-self.knormalized[2]**2)
            self.k=(self.knormalized[0]*ctheta+stheta/rad*(self.knormalized[0]*self.knormalized[2]*cphi-self.knormalized[1]*sphi),
                    self.knormalized[1]*ctheta+stheta/rad*(self.knormalized[1]*self.knormalized[2]*cphi-self.knormalized[0]*sphi),
                    self.knormalized[2]*ctheta - rad*stheta*cphi)

class PhoXtal(object):
    def __init__(self, photon, xtal):
        """Init function for the class defining the interactions between a photon and a crystal. It is really simple at the moment..."""
        self.photon=photon
        self.xtal=xtal

    def g_needed(self, order=1):
        """Calculates the reciprocal lattice vector that a photon needs to be
        scattered by a known crystal with planes oriented in a particular
        direction that is already known"""
        #When the energy is too low the returned value is None.
        if self.photon.knorm < self.xtal.gnorm*order: return None
        # For the calculation see the notes of 09/07/2004 modified the 12/07/04
        A = (self.photon.k[0]*math.cos(self.xtal.gphi)+self.photon.k[1]*math.sin(self.xtal.gphi))
        B = self.photon.k[2]
        C = self.xtal.gnorm/2.
        if abs(A)<=EPSILON:
            # B can't be 0 together with A unless gnorm=0, so this function should be almost safe
            # The result can be otained with the calculation in the else block
            # but this is faster and common (is the case of paraxial photons).
            if abs(C/B)>1.: return None
            else: theta=math.acos(-C/B)
        else:
            A2B2=A**2+B**2
            BC=B*C
            # Two solutions, but only one satisfies the Laue equation
            Mean=-BC/A2B2
            try:
                Delta=A*math.sqrt(A2B2-C**2)/A2B2
            except ValueError:
                return None
            cosp, cosm = Mean+Delta, Mean-Delta
            thp, thm = math.acos(cosp), math.acos(cosm)
            sinp, sinm = math.sin(thp), math.sin(thm)
            checkp = abs(A*sinp + B*cosp + C)
            checkm = abs(A*sinm + B*cosm + C)
            # Da ricontrollare
            if  checkp<checkm:
                theta=thp
            else:
                theta=thm
        return Physics.spherical2cartesian(self.xtal.gnorm, self.xtal.gphi, theta)

    def sigma(self, G=False, order=1, distribution="gaussian"):
        """Interface to the sigma function from Physics."""
        tB = self.xtal.keV2Bragg(self.photon.e)
        ang_fac = Physics.angular_factor(tB)
        Q = self.xtal.mat_facs[order]*ang_fac
        if G is None: G=self.g_needed(order)
        Dtheta = Physics.anglebetween(G, self.xtal.g)
        if distribution=="hat":
            weight=Physics.hat(Dtheta, Physics.eta2fwhm(self.xtal.eta))
        elif distribution=="gaussian":
            weight=Physics.gaussian(Dtheta, self.xtal.eta)
        else: return 0.
        cosine = abs(self.photon.k[2]/self.photon.knorm)
        return Q*weight/cosine

    def sigmas(self):
        '''Returns the probabilities per unit of length that a photon is diffracted.'''
        MAX_ORDER=5 # Maximum diffraction order
        sigmas_=[]
        for i in range(1, MAX_ORDER+1):
            # Forward diffraction
            G = self.g_needed(order=i)
            if G is None: sigma_= 0.
            else: sigma_=self.sigma(G, order=i)
            sigmas_.append(sigma_)

            # Backward diffraction: reversing for a while xtal.g!!!!!
            self.xtal.g=-self.xtal.g
            G = self.g_needed(order=i)
            if G is None: sigmaback=0.
            else: sigmaback=self.sigma(G, order=i)
            self.xtal.g=-self.xtal.g # xtal.g returns normal              !!!!!
            sigmas_.append(sigmaback)
        return sigmas_

    def CSs(self):
        '''Gives all the CSs of a photon in a Xtal'''
        CSs_=[Physics.CS_Photo(self.xtal.Z, self.photon.e),
              Physics.CS_Compt(self.xtal.Z, self.photon.e),
              Physics.CS_Rayl(self.xtal.Z, self.photon.e)]
        # The first three elements depend only because the Photon interacts with a material
        # The other are the scattering cross section coming out from diffractin theory
        CSs_.extend([x/xtal.rho for x in self.sigmas(xtal)])
        return CSs_

    def isinside(self):
        '''Boolean and self explaining... Works if the crystal has a
        parallelepipedon shape.'''
        if not 0.-EPSILON <= -self.photon.r[2] <= self.xtal.dim[2]+EPSILON: return False
        if not abs(self.photon.r[0])<=self.xtal.dim[0]/2+EPSILON: return False
        if not abs(self.photon.r[1])<=self.xtal.dim[1]/2+EPSILON: return False
        else: return True

    def Compton(self):
        '''Penelope MC algorithm for Compton scattering [Penelope]'''
        phi=uniform(0.,2.*pi)
        k=self.photon.e/Booklet.MEC2
        k2=2.*k
        a_1 = log(1.+k2)
        a_2 = k2*(1.+k) / (1. + k2)**2
        tau_min=1./(1.+k2)
        while True:
            if uniform(0., a_1+a_2)<a_1: tau=tau_min**uniform(0.,1.)
            else: tau=math.sqrt(tau_min**2+uniform(0.,1.)*(1.-tau_min**2))
            costheta=1.-(1.-tau)/(k*tau)
            Tnum=(1.-tau)*((k2+1.)*tau-1.)
            Tden=k**2*tau*(1.+tau**2)
            Tcostheta=(1.- Tnum/Tden)
            if uniform(0.,1.)<Tcostheta:
                theta=acos(costheta)
                self.photon.deviate(theta, phi)
                self.photon.e=self.photon.e/(1.+k*(1.-costheta))
                break

    def Rayleigh(self):
        '''Scattering Rayleigh calculated through rejection method'''
        while True:
            theta=uniform(0, pi)
            ctheta, stheta  = math.cos(theta), math.sin(theta)
            FF=Physics.FF_Rayl(self.xtal.Z, stheta/self.photon.lam)
            ptheta = (1. + ctheta**2) * FF**2
            if ptheta <=uniform(0., self.xtal.Z**2):
                self.photon.deviate(theta, uniform(0., 2.*pi))
                break

    def which_interaction(self):
        '''Decide the interaction of the photon that interacts after a path t.'''
        CSs= self.CSs()
        # Calculate interaction point and move the photon there
        mu = sum(CS*self.xtal.rho for CS in CSs)
        t = Physics.free_path(mu)
        self.photon.travel(t)
        # If the photon exits return -1 (No interaction) else evaluate interaction type
        if not self.isinside(): return -1
        else:
            partial_sums, total_sum = [0.], sum(CSs)
            # Random variable to discriminate the interaction
            rnd=uniform(0., total_sum)
            for i, CS in enumerate(CSs):
                partial_sums.append(partial_sums[-1]+CS)
                if rnd < partial_sums[-1]: return i

    def interact(self, interaction_type):
        '''Once known the interaction proceed to the Photon treatment.'''
        if interaction_type==0: pass # Dead photon do nothing
        elif interaction_type==1: self.Compton() # Compton
        elif interaction_type==2: self.Rayleigh() # Rayleigh
        elif interaction_type==-1: pass # Exits without interaction
        # if i>2 we have diffraction (frontward if i is odd).
        # i/2 or i/2-1 is the order of the diffraction in the two cases
        elif interaction_type%2:
            G=xtal.g*(interaction_type/2)
            self.photon.k=self.photon.k+G # Diffraction frontward
        else:
            G=-self.xtal.g*(interaction_type/2-1)
            self.photon.k=self.photon.k+G # Diffraction Back
        self.int_history.append(interaction_type) # Update history

    def interaction(self):
        """Interaction between photon and xtal."""

        # Lambda function for coordinate transformation
        theta, phi = pi/2-self.xtal.gtheta, self.xtal.gphi
        drag_in_xtal = lambda x: x-self.xtal.r
        drag_out_xtal = lambda x: x+self.xtal.r
        rot_in_xtal  = lambda x: Physics.rot(x, (('z', phi), ('y', theta)))
        rot_out_xtal = lambda x: Physics.rot(x, (('y', -theta), ('z', -phi)))
        in_xtal = lambda x: rot_in_xtal(drag_in_xtal(x))
        out_xtal= lambda x: drag_out_xtal(rot_out_xtal(x))

        # Changing reference frame for photon and xtal
        self.photon.k=rot_in_xtal(self.photon.k)
        self.photon.r = in_xtal(self.photon.r)
        self.xtal.g=rot_in_xtal(self.xtal.g)

        self.int_history=[] # History init

        # Photon adventure begins...
        while self.isinside():
            interaction_type=self.which_interaction()
            self.interact(interaction_type)
            if interaction_type==0: break # Photoabsorption: the photon dies in to xtal

        # Photon adventure is over...
        if interaction_type==0: pass # Dead photon (photo-absorption)
        elif self.photon.knormalized[2]<=0.: # photon points towards focal plane... (no backscattering)
            # Back to the original reference frame
            self.photon.k=rot_out_xtal(self.photon.k)
            self.photon.r=out_xtal(self.photon.r)
            # transport the photon to the focal plane
            self.photon.travelz(0.)
        else: pass # Photon is backscattered and thus lost

        # anyway the xtal goes backs to the previous reference frame
        self.xtal.g=rot_out_xtal(self.xtal.g)

    def erandom(self, Erange=(1.,1000.), Exp=False):
        """Random energy from powerlaw spectrum."""
        if Exp:
            R = uniform(0., 1.)
            Exp_ = 1.+Exp
            Emin, Emax = (E**Exp_ for E in Erange)
            E=(Emin + R*(Emax-Emin))**(1./Exp_)
        else:
            E = uniform(Erange[0], Erange[1])
        self.photon.e=E

    def generator(self, Emin=1., Emax=1000., Exp=False,  sphi=0., stheta=pi):
        '''Generates randomly a photon that will go on to xtal.
        Exp is the powerlaw spectral index (-2.1 for the Crab Nebula).
        '''
        theta, phi = pi/2-self.xtal.gtheta, self.xtal.gphi
        # lambda functions #
        drag_out_xtal = lambda x: array(x) + array(self.xtal.r)
        rot_out_xtal = lambda x: Physics.rot(x, (('y', -theta), ('z', -phi)))

        # Puts a photon an EPSILON under the xtal surface in a random 2D point
        random_pos_on_surface=[uniform(-x/2., x/2.) for x in self.xtal.dim[:2]]
        random_pos_on_surface.append(-EPSILON)
        self.photon.r=random_pos_on_surface

        # The photon position in observer reference frame (rotation + translation)
        self.photon.r=rot_out_xtal(self.photon.r)
        self.photon.r=drag_out_xtal(self.photon.r)
        if self.photon.i !=1: self.photon.i=1.
        # wavevector definition from source angular parameters
        self.photon.k=Physics.makenormalvec(sphi, stheta)
        # Energy calculation (powerlaw if Exp is False)
        self.erandom((Emin, Emax), Exp=Exp)

    def diffracted_eranges(self, Emin=1., Emax=1000., sphi=0., stheta=pi):
        """Energy ranges that are diffracted by the crystals."""
        order, reached_max_order = 0,False
        Deltatheta=2.5*(self.xtal.fwhm/60/180*pi)
        thetamin, thetamax = self.xtal.gtheta-Deltatheta, self.xtal.gtheta+Deltatheta

        constant=-Booklet.hc/(2.*self.xtal.d_hkl)
        eranges=[]
        while not reached_max_order:
            order+=1
            cosdeltaphi=math.cos(sphi-self.xtal.gphi)
            sinstheta=math.sin(stheta)
            cosstheta=math.cos(stheta)
            sinthetamin=math.sin(thetamin)
            costhetamin=math.cos(thetamin)
            sinthetamax=math.sin(thetamax)
            costhetamax=math.cos(thetamax)
            Emin_= constant*order/(cosdeltaphi*sinstheta*sinthetamin+cosstheta*costhetamin)
            Emax_= constant*order/(cosdeltaphi*sinstheta*sinthetamax+cosstheta*costhetamax)
            if Emax_ >= Emax:
                if Emin_>=Emax: return eranges
                else:
                    reached_max_order=True
                    Emax_=Emax
            if Emin_ <= Emin:
                if Emax_ <= Emin: pass
                else:
                    Emin_=Emin
                    eranges.append((Emin_, Emax_))
            else: eranges.append((Emin_, Emax_))
        return eranges

    def generator_fast(self, Emin=1., Emax=1000., Exp=0.,  sphi=0., stheta=pi):
        """Generates a photon that is diffracted for sure."""
        self.photon.i=0.
        eranges=self.diffracted_eranges(Emin=Emin, Emax=Emax, sphi=sphi, stheta=stheta)
        weigths, cumweights = [], []
        partial_sums=0
        Exp_=Exp+1.
        # normalization costant proportional to the incident photons
        normalization=1./(Emin**(Exp_)-Emax**(Exp_))
        for i, erange in enumerate(eranges):
            # calculating weights, i.e. fraction of photons considered vs all
            # each weight is calculated for a diffraction order
            # cumweights are the sums of the weights
            weight=normalization*(erange[0]**(Exp_)-erange[1]**(Exp_))
            weigths.append(weight)
            partial_sums+=weight
            cumweights.append(partial_sums)
        # random number to decide which are the orders to consider
        rnd=uniform(0., partial_sums)
        for i, cumweight in enumerate(cumweights):
            if rnd<=cumweight:
                self.generator(Emin=eranges[i][0], Emax=eranges[i][1],
                                Exp=Exp, sphi=sphi, stheta=stheta)
                self.photon.i=weigths[i]
                break

    def fast_diffraction(self):
        """Coerce a photon to diffract, changing it's intensity."""
        # Calculate most efficient diffraction order
        order, sigma_ = 0,0.
        i=0.
        for sigma in self.sigmas():
            i+=.5
            if sigma > sigma_: order, sigma_ = int(i)+1, sigma
        self.photon.k=self.photon.k+self.xtal.g*order
        T=self.xtal.dim[2]
        PoS=0.5*math.exp(-Booklet.mu(self.xtal.Z, self.photon.e)*T)*(1.-math.exp(-sigma_*T))
        self.photon.i = self.photon.i*PoS

class Source(PhoXtal):
    def __init__(self, phi=0., offset=0., powerindex=0., crab=1.):
        """Class defining X-ray source with some characteristics"""
        pass

if __name__=='__main__':
    ph=Photon()
    print ph.lam, ph.e, ph.k, ph.r, ph.knormalized, ph.knorm
    pass