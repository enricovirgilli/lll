# -*- coding: iso-8859-15 -*-
import os, time, thread, numpy
from random import uniform, gauss
import string
#from pylab import plot, show, load, imshow
import pyfits
import numpy as np
import Xtal, Source, Detector, Physics
from mycounter import Progressbar
from numpy import convolve, save, savetxt
# Playing with namespaces
Physics=Xtal.Physics
sys=Physics.sys
Booklet=Physics
math=Physics
pi = math.pi
save=Detector.pylab.save
zeros=Detector.pylab.zeros
array=Detector.pylab.array
arange=Detector.pylab.arange
import time
import sys
from progressbar import  ProgressBar, Percentage, Bar
from Tkinter import *
import threading


class progressBar():
    def __init__(self):
        self.pp=pBar()
        self.pp.percent=0
        self.pp.start()
    def setPercent(self,perc):
        self.pp.setPercent(perc)
    def getPercent():
        return self.pp.percent

    
#http://www.pierotofy.it/pages/guide_tutorials/Python/ProgressBar_con_GUI__Python/
class pBar(threading.Thread):
    def run(self):
        self.main=Tk()
        self.frame=Frame(self.main)
        self.frame.pack()
        self.label = Label(self.frame,text="Photons")
        self.label.pack()
        self.title="Photons generation"
        self.main.title(self.title)
        self.imm=Canvas(self.frame,height=50,width=300)
        self.imm.pack()
        self.back=self.imm.create_rectangle(1,3,300,20)
        self.barr=self.imm.create_rectangle(1,1,1,20,fill="darkgreen")
        self.refresh()
        self.main.mainloop()
    def refresh(self):
        self.imm.coords(self.barr,1,1,self.percent*3,20)
        self.label.configure(text=str(self.percent)+"%")
        self.main.title(self.title+" - "+str(self.percent)+"%")
        if self.percent>=100:
            self.main.quit()
        self.main.after(100, self.refresh)
    def setPercent(self,perc):
        self.percent=perc




class Generic(Xtal.Xtal):
    """Definition of the class ottica that should contain, even with redundance
    ALL the informations needed to calculate EVERYTHING

    TODO: --
    """
    def __init__(self,
                 name = "lens",
                 geo = "spiral",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor = 1.,
                 datadir="data",
                 Focal = 2000.,                 # [cm]
                 source_lens_distance = 2600.,  # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04,        # [cm] diameter of the X ray tube source focus (valid for finite source)
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 Emin=100.,                    # [keV]
                 Emax=600.,                    # [keV]
                 dim=[1., 1., 0.4],            # [cm,cm,cm]
                 rad_framewidth = 0.20,        # [cm] radial and tangential framewitdh,
                 tan_framewidth = 0.20,        # [cm] they can be different
                 curvature =  1/4000.,         # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 external_curvature = 40.,
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,                      # [arcmin] mosaicity 
                 hkl=[1,1,1],                  # CURVCORRECTION ## I erased the flag external_curvature=40.
                 Z=29,                         # it was made for curved crystals, for flat crystal does not affect.
                 microthick=0.,                # Actually it depends on Focal, then is a self.variable=Focal*2
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        """Init of the class. Every parameter has already a default.
        The user is invited to call this constructor naming all parameters, for the
        sake of clarity.

        INPUT PARAMETERS
        - name: lens name mainly used for the output
        - geo: geometry (rings, spiral, ...)
        - datadir: output directory
        - Focal: lens focal length [cm]
        - Emin,Emax: nominal lens working bandwidth [keV]
        - dim: xtal tile size [cmxcmxcm]
        - rad_framewidth: extraspace around each tile [cm] in the radial direction
        - tan_framewidth: extraspace around each tile [cm] in the tangential direction
        - fwhm: mosaic spread [arcmin]
        - hkl: Miller indices
        - Z: atomic number
        - microthick: microblocks thickness [A]
        - pixels: detector pixels
        - pitches: detector pixel size [cm]
        - error: FWHM of the misplacement of crystal [arcmin]

        TODO:
        - Extracting right values of tBmin, tBmax, rmin and rmax from the tile list
        """
        Xtal.Xtal.__init__(self, r=(0.,0.,0.), g=(0.,0.,0.),
                            dim=dim, Z=Z,
                            hkl=hkl, fwhm=fwhm, microthick=microthick)
        # rules to build the lens
        self.name = name
        self.geo  = geo
        self.profile = profile
        self.structure=structure
        self.thickness_profile = thickness_profile
        self.thickness_threshold = thickness_threshold
        self.thickness_factor = thickness_factor
        self.error=Booklet.arcmin2rad(error)
        self.Focal = Focal
        self.source_lens_distance = source_lens_distance
        self.source_dimension = source_dimension
        self.detector_distance = detector_distance 
        self.Emin = Emin
        self.Emax = Emax
        self.rad_framewidth = rad_framewidth
        self.tan_framewidth = tan_framewidth
        self.curvature = curvature
########################################################## CURVCORRECTION
        self.external_curvature = external_curvature #37. #(2./100.)*self.Focal   ## external_curvature is 2 * Focal, but everywhere
########################################################## I used meters instead of cm, then 100 is to convert 
        self.datadir=datadir                     ######### Focal (cm) to (m)
        #
        # Here I begin to calculate some useful
        # values on the basis of the input parameters
        #
        self.tBmin = self.keV2Bragg(self.Emax) # Minimum and maximum Bragg angles (approximation) [radians]
        self.tBmax = self.keV2Bragg(self.Emin)
        self.rmin = 2.*self.Focal*math.tan(self.tBmin) # Minimum and maximum radii of the lens (approximation) [cm]
        self.rmax = 2.*self.Focal*math.tan(self.tBmax)
        thicknesses=[self.dim[2] for x in self.xtals()] # A little bit long calculation
        self.noofxtals=len(thicknesses)
        self.averagethickness = sum(thicknesses)/self.noofxtals
        self.volume = self.cross_section*self.averagethickness*self.noofxtals
        self.weight = self.volume*Booklet.density[self.Z]*1.e-3
        self.filling_factor = (self.noofxtals*self.cross_section)/(pi*((self.rmax+dim[0]/2+rad_framewidth/2)**2-(self.rmin-dim[0]/2-rad_framewidth/2)**2))
        # keV dependent parameters
        self.keVlist = []
        self.arealist = []
        self.senslist = []
        # Detector part
        cartesian=Detector.PSD_cartesian(pixels=pixels, pitches=pitches)
        radius=math.sqrt(cartesian.sides[0]**2+cartesian.sides[1]**2)
        polar=Detector.PSD_polar(int(radius/min(pitches)), min(pitches), 24)
        self.detector = Detector.PSD(cartesian, polar)

    def __str__(self, filename=False):
        """Displays ottica values"""
        dim=(self.dim[0], self.dim[1], self.averagethickness)
        string=("Geometry:                 %s" % self.geo,
                "Profile:                  %s" % self.profile,
                "Thickness profile:        %s" % self.thickness_profile,
                "Thickness threshold:      %s" % repr(self.thickness_threshold),
                "Thickness factor:         %s" % self.thickness_factor,
                "Focal Length:             %s m" % (self.Focal/100),
                "Source Lens distance:     %s" % self.source_lens_distance,
                "Source Dimension:     %s" % self.source_dimension,
                "Lens Detector distance:   %s" % self.detector_distance,
                "Inner/Outer radius:       %s/%s m" % (self.rmin/100, self.rmax/100),
                'Z (Atomic number):        %s' % self.Z,
                "Energy range:             [%g, %g] keV" % (self.Emin, self.Emax),
                'hkl:                      (%d,%d,%d)' % tuple(self.hkl),
                "d_hkl:                    %s A" % self.d_hkl,
                'Cell volume:              %s A^3' % Booklet.volume[self.Z],
                'Structure factor:         %s' % self.sf,
                'Microcrystal size:        %s micron' % (self.microthick/10000),
                'A_0 constant:             %s keV' % self.A_0const,
                'Material factor constant: %s cm^-1' % self.mat_fac,
                'Crystal tile size:        %gx%gx%g cm^3' % tuple(dim),
                'Number of crystals:       %s' % self.noofxtals,
                'Filling factor:           %s' % self.filling_factor,
                'Crystals volume:          %s cm^3' % self.volume,
                'Crystals weight:          %s kg' % self.weight,
                'mosaicity (FWHM):         %s arcmin' % self.fwhm,
                'Radial Frame width:       %s cm' % self.rad_framewidth,
                'Tangential Frame width:   %s cm' % self.tan_framewidth,
                'Curvature radius of each tile: %s cm' % self.curvature,
                'External curvature:   %s m' % self.external_curvature,
                'Misorientation (FWHM):    %s arcmin' % Booklet.rad2arcmin(self.error*2.35))
        return "\n".join(string)

    def set_position(self, keV, phi=0.0):
        """Places a crystal with the optimum orientation to concentrate photons
        of energy keV on the focal spot for a Focal distace equal to F"""
        self.phi=phi #+ 0.52
        tB=self.keV2Bragg(keV)

        if self.profile=="spherical":
            self.z=(2.*math.cos(tB)-1.)*self.Focal
            self.rho=2.*self.Focal*math.sin(tB)
 
        elif self.profile=="parabola":
            tg2tB=(math.tan(2.*tB))**2
            self.z=2.*self.Focal*(math.sqrt(1.+tg2tB)-1)/tg2tB
            self.rho=self.z*math.tan(2.*tB)

        elif self.profile=="flat":

            self.z=self.Focal
            self.rho=self.Focal*math.tan(2.*tB)

        else:
            sys.stderr.write("Wrong profile!\nExiting...\n")
            sys.exit()
        gnorm=2.*pi/self.d_hkl
        #print "GNORM", gnorm
        if self.error==0.:
            self.g=Physics.spherical2cartesian(gnorm, self.phi+pi, pi/2.-tB)

        else:
            # gauss is taken from import random. Is a function getting a random number
            # with gaussian probability around the 0.0 value, having width = self.error
            self.g=Physics.spherical2cartesian(gnorm, self.phi+pi, pi/2.-tB+gauss(0.0, self.error))


    def set_thickness(self, keV):
        """Changes crystal thickness according to the lens rules"""
        if self.thickness_profile=="optimized":
            tbest = self.best_thickness(keV)
            t=tbest*self.thickness_factor
            if t <= self.thickness_threshold[0]:
                self.dim[2] = self.thickness_threshold[0]
            elif t >= self.thickness_threshold[1]:
                self.dim[2] = self.thickness_threshold[1]
            else:
                self.dim[2] = t

    def gen_spiral_fixed(self, COUNT=False):
        #
        #
        # function that is also used for the SingleXtal class
        #
        #
        '''Generator for returning crystals that are placed on a spiral and
        with fixed thickness'''
        # kickstart variables

        crystalloffsettato = 0.

        phi = crystalloffsettato
        keV = self.Emax

        # Not so good calculation of the spacing between the tiles
        # TODO: Do it better
        shift_rad = self.dim[0]/2. + self.rad_framewidth
        shift_tan = self.dim[1]/2. + self.tan_framewidth
        spiral_step = self.dim[0] + 2.*self.rad_framewidth

        # START placing the tiles
        while(keV >= self.Emin):
            # This line has to been done because self.rho is used for calculations
            # It should be placed inside the if COUNT... condition to speed up calculations
            self.set_position(keV, phi)
            if COUNT is not True:
                # modifying crystal
                self.set_thickness(keV)
            # Preparing for new step
            phistep= 2.* math.atan((shift_tan) / (self.rho - shift_rad))
            phi += phistep
            next_rho = self.rmin + spiral_step * phi/(2.*pi)
            keV=self.radius2keV(next_rho)
            yield None

    def gen_one_ring(self, keV=False, xtalnumber=False, COUNT=False):
        #######################################################################
        ####### this is the correct gen_one ring which is used for
        ####### colculate the number of crystals foreach ring
        ########################################################################
        '''Generator for returning crystals that are placed on a single ring
        optimized for a defined energy (default Emin+Emax/2)'''
        if not keV: keV = 0.5*(self.Emin+self.Emax)
        if xtalnumber is False:
            tB=self.keV2Bragg(keV)
            #r=2.*self.Focal*math.sin(tB)
            r=self.Focal*math.tan(2*tB)
            num=self.tan_framewidth+self.dim[1]/2
            den=r-self.dim[0]/2 - self.rad_framewidth   ########## correction for number of crystals rad_frame
            xtalnumber=int(pi/math.atan(num/den))
        self.set_position(keV, 0.)
        if COUNT is True:
            for i in range(xtalnumber):
                yield None
        else:
            Delta_phi= 2.*pi/xtalnumber
            for i in range(xtalnumber):
                self.set_position(keV, Delta_phi*i)
                self.set_thickness(keV)
                yield None




    def gen_rings(self, COUNT=False):
        '''Generator for returning crystals that are placed on concentric rings
        to diffract energies in the nominal range [Emin, Emax]
        '''
        keV=float(self.Emax)
        while(keV >= self.Emin):
            # yield xtals of ONE ring
            for xtal in self.gen_one_ring(keV=keV, COUNT=COUNT): yield None
            # Next ring: radius and energy calculation (PhD thesis)
            
            part1 = (self.rho+self.dim[0]/2.+self.rad_framewidth)**2
            part2 = (self.dim[1]/2.+self.tan_framewidth)**2

            sqrt=math.sqrt(part1+part2)
            #####################################################################################
            # NEXTRHO CALCOLATO DA PISA SECONDO ME TROPPO COMPLICATO E SBAGLIATO, CORREGGO SOTTO
            next_rho = self.dim[0]/2. + self.rad_framewidth + sqrt
            #CORREZIONE AL CALCOLO DEL RAGGIO RHO (ENRICO)
            next_rho_enrico = self.rho + self.dim[0] + 2 * self.rad_framewidth
            ####################################################################################
            # next energy
            keV = self.radius2keV(next_rho)


    def radius2keV(self, radius, order=1):
        if self.profile=="spherical":
            angle=math.asin(radius/2./self.Focal)
        elif self.profile=="flat":
            angle=0.5*math.atan(radius/self.Focal)
        return self.Bragg2keV(angle, order)

##############################################################
#    def keV2radius(self, keV, focal=20, order=1):
#        # being hc and d_hkl expressed in Angstrom, focal
        # is espressed in cm then r will be in cm 
#        d_hkl= Booklet.d_hkl(self.Z, self.hkl)
#        result = ( Booklet.hc * self.Focal ) / (d_hkl * keV)
#        return result
###############################################################

    def xtals(self, xtalnumber=False, COUNT=False):
        '''Generator that yields the crystals.
        Fast and not memory demanding'''
        if self.geo=="spiral":
            xtals_generator=self.gen_spiral_fixed
        elif self.geo=="one ring":
            xtals_generator=self.gen_one_ring
        elif self.geo=="rings":
            xtals_generator=self.gen_rings
        elif self.geo=="lauegame":
            xtals_generator=self.gen_laue_game
        else:
            sys.stderr.write('Wrong geometry!\nExiting...\n')
            sys.exit()
        for xtal in xtals_generator(COUNT=COUNT):
            yield None

    def info(self, filename=False):

        """Displays data using external to the class function display_value,
        TODO:
        - pretty showing all in one
        - auto resize of the table in function of the output
        - add the chance to exclude values
        - output on file
        """
        if filename: sys.stdout=file(filename, 'w')
        sys.stdout.write(self.__str__())
        sys.stdout=file('/dev/stdout', 'w')
################################################################################
    def smooth(self,lista):
        degree=2
        window=degree*2-1  
        weight=numpy.array([1.0]*window)  
        weightGauss=[]  
        for i in range(window):  
            i=i-degree+1  
            frac=i/float(window)  
            gauss=1/(numpy.exp((4*(frac))**2))  
            weightGauss.append(gauss)  
        weight=numpy.array(weightGauss)*weight  
        smoothed=[0.0]*(5-window)  
        for i in range(len(smoothed)):  
            smoothed[i]=sum(numpy.array(lista[i:i+window])*weight)/sum(weight)  
        return lista
#################################################################################

    def EA_ring_at_keVs(self, xtal_in_ring=1., NUM_OF_SIGMA=4):
        """Calculates effective area for a ring known the number of crystals
        composing it.
        NOTE: works only on axis!!!"""
        # Finding incidence angle with respect to crystal planes
        tB = math.pi/2-self.gtheta # always on axis
        # determining maximum diffraction order
        max_order = int(self.keVlist[-1]/self.Bragg2keV(tB)+1)
        # geometrical area of the ring
        ring_cross_section=self.cross_section * xtal_in_ring
        #for  i,keV in enumerate(self.keVlist):
            #Physics.best_thickness(keV, self.Z, self.hkl, self.eta)            
        for order in range(1,max_order):
            # execute something only if structure factor is not zero
            if not self.sfs[order]==0.:
                Emin, Emax = self.Erange(order, OFFSET=0., NUM_OF_SIGMA=NUM_OF_SIGMA)
                keVs=[(i,keV) for i,keV in enumerate(self.keVlist) if Emin<=keV<=Emax]
                
                if self.structure=="perfect" and  self.curvature !=0.:
                    poss=[(i, self.peakreflecurved(keV, 0., order)) for i, keV in keVs]
                else:              
                    poss=[(i, self.reflectivity(keV, 0., order)) for i, keV in keVs]
                for i, pos in poss:
                    if pos!=0.: self.arealist[i] = self.arealist[i] + pos * ring_cross_section



    def smoothListGaussian(list,degree=10): 
        window=degree*2-1
        weight=numpy.array([1.0]*window)  
        weightGauss=[] 
        for i in range(window):  
            i=i-degree+1  
            frac=i/float(window)  
            gauss=1/(numpy.exp((4*(frac))**2))  
            weightGauss.append(gauss)  
        weight=numpy.array(weightGauss)*weight  
        smoothed=[0.0]*(len(list)-window)  
        for i in range(len(smoothed)):  
            smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
        return smoothed

    
    def bentperfectcrystals(self,x,eta,fwhm,NUM_OF_SIGMA=4.):
        g = []
        h = []
        width = fwhm
        #if width < 1: width = 1.0 
        lower = int(x - NUM_OF_SIGMA * fwhm)  
        upper = int(x + NUM_OF_SIGMA * fwhm)
        for i in range(0,lower):
            g.append(0.0)
            h.append(Physics.hat1(i, fwhm, x0=x))       
        for i in range(lower,upper):
            g.append(math.exp(-4*(i-x)**2))
            #g.append(self.peakreflecurved(i, 0., 1))
            h.append(Physics.hat1(i, fwhm, x0=x)*self.peakreflecurved(i, 0., 1))
        hatg = convolve(h,g)     
        sl = slice(x, 2000+x, 1)      
        outfile = open("hatg.dat","w")
        for item in hatg[sl]:
            outfile.write("%s\n" % item)
        outfile = open("g.dat","w")
        for item1 in g:
            outfile.write("%s\n" % item1)
        outfile = open("h.dat","w")
        for item2 in h:
            outfile.write("%s\n" % item2)
        return hatg[sl]
               
###############################################################################################
#  To build effective area for a ring made of bend crystal Quasimosaic
###############################################################################################
    def EA_bendring_at_keVs(self, xtal_in_ring=1., NUM_OF_SIGMA=4):
       
        tB = math.pi/2-self.gtheta      # on axis
        
        d_hkl = Booklet.d_hkl(self.Z, self.hkl)

        EB=Physics.BraggEnergy(tB, d_hkl)
    
        # alpha_nominal is needed to renormalize the reflectivity if the total energy band
        # is covered by a single crystal or more than 1. alpha_nominal is the nominal angle
        # of the whole lens you want to build see below the normalization factor 
        # alpha_nominal/alpha (that of 1 single crystal)
        
        #alpha_nominal= 2 * math.atan((self.rmax-self.rmin)/(2*self.external_curvature*100) )
   
        rout= self.rmax + self.dim[0]/2
        rin = self.rmin - self.dim[0]/2

        alpha_nominal = (rout-rin) * self.curvature 

        eta = (self.dim[2] * self.curvature) / (2.59) 
       
        fwhm = self.dim[0] * self.curvature  

        eta_energy=Physics.BraggEnergy(tB-eta/2, d_hkl)-Physics.BraggEnergy(tB+eta/2, d_hkl)

        fwhm_energy=Physics.BraggEnergy(tB-fwhm/2, d_hkl)-Physics.BraggEnergy(tB+fwhm/2, d_hkl)

        #Dthetamax = fwhm#NUM_OF_SIGMA*fwhm
        
        # determining maximum diffraction order                                  # corrected alpha_nominal using rin and rout 
        max_order = 2 #int(self.keVlist[-1]/self.Bragg2keV(tB)+1)                # because in the case of single ring rmin 
        # geometric  area of the ring                                            # rmax (lens ring) are the same i.e. the 
        ring_cross_section=self.cross_section * xtal_in_ring                     # central radius at which the cristal is posed
        Emin = self.Emin - fwhm_energy/2
        Emax = self.Emax + fwhm_energy/2
       
        #for  i,keV in enumerate(self.keVlist):
        for order in range(1,max_order):
            # execute something only if structure factor is not zero
            if not self.sfs[order]==0.:
                alpha=self.dim[0]/(self.external_curvature * 100.) # factor 100 for conversion to cm  
                renorm_factor = (alpha_nominal/alpha)
                ###Emin = self.Bragg2keV(tB + alpha/2)   # Emin & Emax diffract by this curved crystal
                ###Emax = self.Bragg2keV(tB - alpha/2)   # due to dim[0] and to external curvature radius
                keVs=[(i,keV) for i,keV in enumerate(self.keVlist) if Emin<=keV<=Emax]
                a=[]
                #for keV in range (Emin,Emax):
                #    a.append(self.peakreflecurved(keV, 0., order))
                #b=self.smoothListGaussian(int(a))
                poss=[(i, self.peakreflecurved(keV, 0., order)) for i, keV in keVs]
                #poss=[(i, self.smoothListGaussian(a)) for i, keV in keVs]
                #poss=[(i, self.bentperfectcrystals(EB,eta_energy,fwhm_energy,NUM_OF_SIGMA=4.)[keV]) for i, keV in keVs]               
                for i, pos in poss:
                    if pos!=0.: self.arealist[i] = self.arealist[i] + pos * self.dim[1] * self.dim[2] * xtal_in_ring * (1/2.59)  
                    #if pos!=0.: self.arealist[i] = self.arealist[i] + pos * self.dim[1] * xtal_in_ring * (self.external_curvature * 100)

###############################################################################################


    def conf(self,x, y): 
        P, Q, N = len(x), len(y), len(x)+len(y)-1 
        z = [] 
        for k in range(N): 
            t, lower, upper = 0, max(0, k-(Q-1)), min(P-1, k) 
            for i in range(lower, upper+1): 
                t = t + x[i] * y[k-i] 
            z.append(t) 
        return z 

###############################################################################################
    def largegaussss(self,x,eta,fwhm,NUM_OF_SIGMA=4.):
        g = []
        h = []
        lower = (x - NUM_OF_SIGMA * fwhm/2)  
        upper = (x + NUM_OF_SIGMA * fwhm/2)
        l10 = int(lower*10)
        u10 = int(upper*10)
        for k in range(0,l10):
            i = k/10.
            g.append(0.0)
            h.append(0.0)
        for k in range(l10,u10):
            i = k/10.
            muT=self.mu(i) * self.dim[2]
            g.append(self.reflectivitywithoutabsorption(i, 0., 1))
            h.append(Physics.hat(i, fwhm, x0=x)* math.exp(-muT))
        #hatg = self.conf(h,g)
        hatg =convolve(h,g)
        a = []
        b = []
        c = []
        norm = max(g)/max(hatg)
        for d in range(len(hatg)):
            a.append(hatg[d]*norm)
        sl = slice(int(x*10), 10000,1)
        for g in range(0,len(a[sl])):
            b.append(a[sl][g])
        for f in range(0,len(b),10):
            c.append(b[f])   
        return c#b#b[sl]
######################################################

    def largegaussssfake(self,x,eta,fwhm,NUM_OF_SIGMA=4.):
        g = []
        h = []
        width = fwhm
        if width < 1: width = 1.0 
        decimali = x - int(x)
        lower = (x - NUM_OF_SIGMA * fwhm)  
        upper = (x + NUM_OF_SIGMA * fwhm)
        l10 = int(lower * 10)
        u10 = int(upper * 10)
        for k in range(10,l10):
            i = k/10. -1
            g.append(0.0)
            h.append(Physics.hat(i, fwhm, x0=x))       
        for k in range(l10,u10):
            i=k/10. -1
            muT=self.mu(i) * self.dim[2]
            g.append(self.reflectivitywithoutabsorption(i, 0., 1))
            h.append(Physics.hat(i, fwhm, x0=x) * math.exp(-muT))
        hatg = convolve(h,g)
        sl = slice(x*10, 2000+(x*10), 1)
        return hatg[sl]#*fwhm


    
###############################################################################################
#  To build effective area for a ring made of bent mosaic crystals 
###############################################################################################
    def EA_bentmosaic_at_keVs(self, xtal_in_ring=1., NUM_OF_SIGMA=4):
        # Finding incidence angle with respect to crystal planes
        tB = math.pi/2-self.gtheta # always on axis
        EB=Physics.BraggEnergy(tB, Booklet.d_hkl(self.Z, self.hkl))
        th0=0
        #eta = self.fwhm #0.0002  ###per la gaussiana  ##################�ATTENZIONE
        sigma = (self.fwhm/60.)*(math.pi/180.)#################*(1/2.35)
        hatwidth = self.dim[0]/(self.external_curvature * 100.)
        #fwhmm = self.dim[0]/(self.external_curvature * 100.) #0.003 ### per la cassetta
        Dthetamax = NUM_OF_SIGMA*hatwidth
        deltaEgauss=Physics.BraggEnergy(tB-sigma/2, Booklet.d_hkl(self.Z, self.hkl))-Physics.BraggEnergy(tB+sigma/2, Booklet.d_hkl(self.Z, self.hkl))
        deltaEhat=Physics.BraggEnergy(tB-hatwidth/2, Booklet.d_hkl(self.Z, self.hkl))-Physics.BraggEnergy(tB+hatwidth/2, Booklet.d_hkl(self.Z, self.hkl))
        #rout= self.rmax + self.dim[0]/2
        #rin = self.rmin - self.dim[0]/2
        max_order = 2 #int(self.keVlist[-1]/self.Bragg2keV(tB)+1)  
        ring_cross_section=self.cross_section * xtal_in_ring               
        for order in range(1,max_order):
            # execute something only if structure factor is not zero
            if not self.sfs[order]==0.:
                ####external_curvature = 40.   # meters
                #alpha=self.dim[0]/(self.external_curvature * 100.) # factor 100 for conversion to cm 
                #internal_spread = self.fwhm
                Emin=Physics.BraggEnergy(tB+Dthetamax/2, Booklet.d_hkl(self.Z, self.hkl))
                Emax=Physics.BraggEnergy(tB-Dthetamax/2, Booklet.d_hkl(self.Z, self.hkl))
                #Emin = self.Bragg2keV(tB + alpha/2)   # Emin & Emax diffract by this curved crystal
                #Emax = self.Bragg2keV(tB - alpha/2)   # due to dim[0] and to external curvature radius
                keVs=[(i,keV) for i,keV in enumerate(self.keVlist) if Emin<=keV<=Emax]
                poss=[(i, self.reflectivity(keV, 0., order)) for i, keV in keVs]
                #teposs=[(i,Physics.largegauss(EB,deltaEgauss,deltaEhat,NUM_OF_SIGMA=3.)[keV]) for i, keV in keVs]
                teposs=[(i,self.largegaussss(EB,deltaEgauss,deltaEhat,NUM_OF_SIGMA=4.)[int(keV)]) for i, keV in keVs]
                for i, pos in teposs:
                    if pos!=0.: self.arealist[i] = self.arealist[i] + pos * (self.external_curvature * 100) * sigma * self.dim[1] * xtal_in_ring   
###############################################################################################


    def EA_xtal_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4.):
        """Calculates effective area for a single crystal."""
        # Finding incidence angle with respect to crystal planes
        tB = self.gettB(OFFSET)
        # determining maximum diffraction order
        max_order = int(self.keVlist[-1]/self.Bragg2keV(tB)+1)
        for order in range(1,max_order+1):
            # check if this order has a non zero structure factor
            if not self.sfs[order]==0.:
                Emin, Emax = self.Erange(order, OFFSET=OFFSET, NUM_OF_SIGMA=NUM_OF_SIGMA)
                keVs=((i,keV) for i,keV in enumerate(self.keVlist) if Emin<=keV<=Emax)
                poss=((i, self.reflectivity(keV, OFFSET, order)) for i, keV in keVs)
                for i, pos in poss:
                    if pos!=0: self.arealist[i] = self.arealist[i] + pos * self.cross_section
        

    def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4):
        """Calculates effective area for the energy values defined by self.keVlist"""
        self.arealist = zeros(len(self.keVlist), 'f')
        # monitoring tools
        countercondition= (self.noofxtals > 1000)
        if countercondition:
            counter=Progressbar(update_time=1.,
                                    steps=self.noofxtals,
                                    title=self.name,
                                    name='Effective\ area',
                                    )
        # effective calculation
        for xtal in self.xtals():
            self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)
            if countercondition: counter.index+=1
        if countercondition: counter.close()

    def EA_at_keV(self, keV):
        """Effective area [cm^2]."""
        reflex, max_order = 0., 10
        for xtal in self.xtals():
            tB = pi/2-self.gtheta
            for order in range(1,max_order+1): reflex += self.reflectivity(keV, tB, order)
        return reflex * self.cross_section

    

    def sensitivity_at_keV(self, keV, areaeff=None, fraction=0.5, radius=None, Tobs=1.e6, nsigma=3., effdet=1., BG="default"):
        """Estimates the sensitivity of the lens at energy keV and returns the number
        photons that are expected to be counted in an observation time Tobs [ph/(cm^2 s keV)].
        """
        eqS=self.Sfactor_at_keV(keV, areaeff=areaeff, radius=radius, fraction=fraction)
        DeltaE=0.5*keV
        if eqS==0: return 1.
        if BG=="default":
            BG = Booklet.bg(keV)
        elif BG=="ezio":
            BG = Booklet.spibg.read(keV)/Booklet.mu("CdTe", keV)
        return nsigma/effdet*math.sqrt(2.*BG / (Tobs*DeltaE*eqS))

    def Gfactor_at_keV(self, keV, areaeff=None, radius=None, fraction=0.5):
        """Focusing factor at a given energy (dimensionless)."""
        if areaeff is None: areaeff=self.EA_at_keV(keV)
        if radius is None: radius=self.detector.rforfraction(fraction)
        Ad=pi*radius**2
        return fraction*areaeff/Ad

    def Sfactor_at_keV(self, keV, areaeff=None, radius=None, fraction=0.5):
        """Equivalent surface: G^2*Ad. [cm^2]"""
        if areaeff is None: areaeff=self.EA_at_keV(keV)
        if radius is None: radius=self.detector.rforfraction(fraction)
        Ad=pi*radius**2
        return (fraction*areaeff)**2/Ad



    
    def mc_fast(self, stheta=pi, sphi=0., Erange=False, Exp=0., photons=1., oufile='/dev/stdout', watch=False):

        '''stheta and sphi are the source angular parameters, Erange is the Energy range of the source flat 
        spectrum (defaults to self.Emin/2, 2*self.Emax). vainf. the total Photons simulated and randomly 
        distributed all over the crystals'''
        #stheta = pi - phoxtal.divergence_angle()
        # Determining spectrum energy range
        if Erange is False: Emin, Emax = 0.5*self.Emin, 2.*self.Emax
        else: Emin, Emax = Erange
        # Flow controls       
        #if watch:
        #    thread.start_new_thread(self.mc_control_thread, ())
        #counter=Progressbar(update_time=1., steps=int(photons), title=self.name, name='Monte Carlo')
        
        # Starting simulation
        photon=Source.Photon()
        
        ppx=float(photons)/self.noofxtals   # Photons per crystals
        ppx_frac, ppx_int = math.modf(ppx)  # give integer and fractional part of number of photons
        xtalnow=0
        x = 1
        ph  = 0
        #app=progressBar()
        curr=0.1
        for xtal in self.xtals():
            
            phoxtal=Source.PhoXtal(photon, self)
            
            xtalnow+=1
            # Determing number of photons for this crystal
            if uniform(0,1) < ppx_frac: 
                ppx_now=int(ppx_int+1)
            else: 
                ppx_now=int(ppx_int)
            # Diffracting these photons
            
            for i in xrange(ppx_now):
                #counter.index+=1
                deltasource = uniform ( - self.source_dimension/2 , + self.source_dimension/2 )
                phoxtal.generator_fast(Emin=Emin, Emax=Emax, Exp=Exp, stheta=stheta, sphi=sphi, deltasource = deltasource)
                
                if phoxtal.photon.i!=0.:
                    #PROGRESSBAR COUNTER 100/(number of total photons)
                    curr=curr+100/photons
                    # ORIGINAL DIFFRACTION FUNCTION ONLY FOR FLAT CRYSTALS. IF THE fast_diff_bent WORKS 
                    # PROPERLY THIS FUNCTION WILL BE DELETED
                     
           
                    # This function include both mosaic or perfect crystals 
                    # with the possibility to select between flat or bent

                    phoxtal.fast_diffraction()       
                    #phoxtal.fast_diff_bent(Focal=self.Focal)
                    #phoxtal.fast_diff_bent(Focal=self.Focal, Emin=Emin, Emax=Emax, Exp=Exp, stheta=stheta, sphi=sphi, deltasource = deltasource)
  
                    # Function to move the photons as diffracted and reaching the detector plane (z=0) 
                    difference = self.Focal - self.detector_distance
                    phoxtal.photon.travelz(z = difference)  # era z = 0 
                    
                    nominal_angle = (math.asin(self.rho/(self.Focal)))
                    xshift = difference * math.tan(nominal_angle)
                    yshift = difference * math.tan(nominal_angle)

                    # use this function if you want to shift the detector for having the diffratced
                    #at its center, or the second function if want to have the detector
                    # fixed at the nominal posiiton of 20 m
                    self.detector.measure(photon, xshift, yshift)
                    #self.detector.measure(photon, 0, 0)
                    
                    #app.setPercent(int(curr))
                    ph = ph+1
            x=x+1
        # counter death
        #counter.close()

    def PSF(self, stheta=pi, sphi=0., Erange=False, Exp=0., Nevents=1., oufile='/dev/stdout', watch=False):
        #"""Numerical calculation of the PSF! No MC.
        #TODO: everything!"""
        pass


    def xshift(self):
        difference = self.Focal - self.detector_distance
        nominal_angle = (math.asin(self.rho/(self.Focal)))
        xshift = difference * math.tan(nominal_angle)
        return xshift


    def yshift(self):
        difference = self.Focal - self.detector_distance
        nominal_angle = (math.asin(self.rho/(self.Focal)))
        yshift = difference * math.tan(nominal_angle)
        return yshift

    
        

    # Data writers
    def dump_mcdataevery(self, t=10):
        while True:
            time.sleep(t)
            self.dump_mcdata()
            self.plots()

    def dump_mcdata(self, append=""):
        """Dumps all the informations of the matrix"""
        basename=os.path.join(self.datadir, self.name+'_mc')
        x=array((self.detector.xtics(self.xshift()), self.detector.cartesian.data.sum(axis=0)))
        y=array((self.detector.ytics(self.yshift()), self.detector.cartesian.data.sum(axis=1)))
        
        #############
        center = len(self.detector.cartesian.data[:,0])/2
        xcenter=array((self.detector.xtics(self.xshift()), self.detector.cartesian.data[:,center]))
        ycenter=array((self.detector.ytics(self.yshift()), self.detector.cartesian.data[center,:]))
        #############

        r=array((self.detector.rtics(), self.detector.polar.data.sum(axis=0)))
        encircled=array((self.detector.rtics(), self.detector.polar.data.sum(axis=0).cumsum()))
 
       
        savetxt('%s%s.dat' % (basename, append), self.detector.cartesian.data     , fmt='%g')
        savetxt('%s_x%s.dat' % (basename, append), x.transpose()                  , fmt='%g')
        savetxt('%s_y%s.dat' % (basename, append), y.transpose()                  , fmt='%g')
        savetxt('%s_r%s.dat' % (basename, append), r.transpose()                  , fmt='%g')
        savetxt('%s_xcenter%s.dat' % (basename, append), xcenter.transpose()      , fmt='%g')
        savetxt('%s_ycenter%s.dat' % (basename, append), ycenter.transpose()      , fmt='%g')
        savetxt('%s_encircled%s.dat' % (basename, append), encircled.transpose()  , fmt='%g')




        file('%s%s_R50.dat' % (basename, append), 'w').write('%s' % (self.detector.rforfraction(0.5)))
        
        hdu = pyfits.PrimaryHDU(data=self.detector.cartesian.data)
        hdu.writeto('%s%s.fits' % (basename, append))

    def dump_arealist(self,bands=False, append=""):
        '''Dumps effective area list'''
        if bands:
            filename=os.path.join(self.datadir, self.name+'_area_'+bands+'%s.dat' % append)
        else:
            filename=os.path.join(self.datadir, self.name+'_area%s.dat' % append)
        lines=['%s\t%s' % (x, y) for x,y in zip(self.keVlist, self.arealist)]
        file(filename, 'w').write('\n'.join(lines))

    def dump_senslist(self, append=""):
        '''Dumps sensitivity area list'''
        filename=os.path.join(self.datadir, self.name+'_sens%s.dat' % append)
        lines=['%s\t%s' % x for x in self.senslist]
        file(filename, 'w').write('\n'.join(lines))

    def dump_Slist(self, append=""):
        '''Dumps equivalent surface list'''
        filename=os.path.join(self.datadir, self.name+'_S%s.dat' % append)
        lines=['%s\t%s' % x for x in self.Slist]
        file(filename, 'w').write('\n'.join(lines))

    def dump_Glist(self, append=""):
        '''Dumps focusing factor list'''
        filename=os.path.join(self.datadir, self.name+'_G%s.dat' % append)
        lines=['%s\t%s' % x for x in self.Glist]
        file(filename, 'w').write('\n'.join(lines))

    def macro(self, INFO=False, XTALINFO=False, ERANGE=False, AREA=False, MC=False, SENS=False, S=False, G=False, PHOTONS=0, OFFSET=0.0, FRACTION=0.5, Tobs=1.e6):
        """Executes a series of useful operations"""
        # name.log is the file where all the information are logged.
        LOGFILENAME=os.path.join(self.datadir, self.name+'.log')
        logfile=file(LOGFILENAME, 'a').write

        MACRONUMBER=len([l for l in file(LOGFILENAME) if l[:5]=="MACRO"])
        APPEND="_%s" % MACRONUMBER

        logfile("MACRO:     %s\n" % MACRONUMBER)
        logfile("TIME:      %s\n" % time.asctime())
        logfile("INFO:      %s\n" % INFO)
        logfile("XTALINFO:  %s\n" % XTALINFO)
        logfile("ERANGE:    %s\n" % repr(ERANGE))
        logfile("AREA:      %s\n" % AREA)
        logfile("MC:        %s\n" % MC)
        logfile("SENS:      %s\n" % SENS)
        logfile("S:         %s\n" % S)
        logfile("PHOTONS:   %s\n" % PHOTONS)
        logfile("OFFSET:    %s\n" % OFFSET)
        logfile("FRACTION:  %s\n" % FRACTION)

        # Lens Informations:
        # name_info.dat stores all the info about the optics
        if INFO:
            INFOFILE=os.path.join(self.datadir, self.name+'_info%s.dat' % APPEND)
            sys.stdout=file(INFOFILE, 'w')
            self.info()
            sys.stdout=file('/dev/stdout', 'w')

        # Xtal Informations:
        # name_xtalinfo.dat stores all the info about the crystal positions
        XTALINFOFILE=os.path.join(self.datadir, self.name+'_xtalinfo%s.dat' % APPEND)
        if XTALINFO:
            logfile('%s: Writing crystal table\n' % time.asctime())
            for xtal in self.xtals(): self.xtal_info(XTALINFOFILE)
            self.plot_sketch(append=APPEND)

        # Changes the energy range for effective calculations
        # keVlist is modified anyway
        if not ERANGE:
            Emin, Emax, Estep = self.Emin, self.Emax, 1.
        elif len(ERANGE)==2:
            Emin, Emax = ERANGE
            Estep=1.                       # I changed the step size from 1 to 0.1 #######################################
        elif len(ERANGE)==3:
            Emin, Emax, Estep = ERANGE
        self.keVlist=arange(Emin, Emax+Estep/2, Estep)

        # Calculation of the effective area:
        # Area is stored in name_area.dat
        if AREA:
            logfile('%s: Calculating effective area\n' % time.asctime())
            self.EA_at_keVs(OFFSET)
            self.dump_arealist(append=APPEND)
            self.plot_area(append=APPEND)

        # Monte Carlo
        # MCFILE=os.path.join(self.datadir, self.name+'_mc%s.dat' % APPEND)
        if MC:
            logfile('%s: Calculating MC\n' % time.asctime())
            #thread.start_new_thread(self.dumpmcdataevery, (5,))
            self.mc_fast(photons=PHOTONS,
                        Erange=(min(self.keVlist), max(self.keVlist)),
                        stheta=pi*(1.-OFFSET/60/180), Exp=0., watch=False)
            self.dump_mcdata(append=APPEND)
            self.plot_mcs(append=APPEND)

        # Calculation of the equivalent Surface
        if (OFFSET==0.0 and (S or SENS or G)):
            logfile('%s: Calculating equivalent surface\n' % time.asctime())
            radius=self.detector.rforfraction(fraction=FRACTION)
            if radius==0.:
                sys.stderr.write("Unable to determine radius!\n")
                radius = self.detector.radius
            self.Slist=[(x, self.Sfactor_at_keV(x, areaeff=y, radius=radius, fraction=FRACTION))
                        for x, y in zip(self.keVlist, self.arealist)]
            if S:
                self.dump_Slist(append=APPEND)
                self.plot_S(append=APPEND)
            if G:
                self.Glist=[(x, self.Gfactor_at_keV(x, areaeff=y, radius=radius, fraction=FRACTION))
                        for x, y in zip(self.keVlist, self.arealist)]
                self.dump_Glist(append=APPEND)
                self.plot_G(append=APPEND)
            if SENS:
                logfile('%s: Calculating sensitivity\n' % time.asctime())
                self.senslist=[(x, self.sensitivity_at_keV(x, areaeff=y, radius=radius, fraction=FRACTION, Tobs=Tobs)) for x, y in zip(self.keVlist, self.arealist)]
                self.dump_senslist(append=APPEND)
                self.plot_sens(append=APPEND)
        # Close logfile
        logfile('%s: End\n\n' % time.asctime())

    # Data plotters
    def plots(self, append=""):
        """Writes all the gnuplot file for plotting matrix and area information"""
        if os.path.isfile(os.path.join(self.datadir, "%s_area%s.dat" % (self.name, append))):
            self.plot_area(append=append)
        if os.path.isfile(os.path.join(self.datadir, "%s_sens%s.dat" % (self.name, append))):
            self.plot_sens(append=append)
        if os.path.isfile(os.path.join(self.datadir, "%s_S%s.dat" % (self.name, append))):
            self.Splot(append=append)
        if os.path.isfile(os.path.join(self.datadir, "%s_mc%s.dat" % (self.name, append))):
            self.plot_mc(append=append)
            self.plot_mcx(append=append)
            self.plot_mcxcenter(append=append)
            self.plot_mcy(append=append)
            self.plot_mcycenter(append=append)
            self.plot_mcr(append=append)
            self.plot_mcxy(append=append)
            self.plot_mcencircled(append=append)
        if os.path.isfile(os.path.join(self.datadir, "%s_xtalinfo%s.dat" % (self.name, append))):
            self.plot_sketch()

    def plot_area(self, append=""):
        """Plot of the effective area"""
        basename="%s_area%s" % (self.name, append)
        file(os.path.join(self.datadir, "%s.plot" % basename), 'w').write("""set sty data histeps
set bmargin 4
set xlabel 'Energy - keV' font 'Times, 24'
set ylabel 'A_{eff} - cm^2' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "%(basename)s.eps"
pl "%(basename)s.dat" u 1:2 t "" lw 3
!okular %(basename)s.eps &
""" % {"basename" : basename,})

    def plot_S(self, append=""):
        basename="%s_S%s" % (self.name, append)
        file(os.path.join(self.datadir, "%s.plot" % basename), 'w').write("""set sty data histeps
set bmargin 4
set xlabel 'Energy - keV' font 'Times, 24'
set ylabel 'Equivalent Surface - cm^2' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "%(basename)s.eps"
pl "%(basename)s.dat" u 1:2 t "" lw 3
!okular %(basename)s.eps &
""" % {"basename" : basename,})

    def plot_G(self, append=""):
        basename="%s_G%s" % (self.name, append)
        file(os.path.join(self.datadir, "%s.plot" % basename), 'w').write("""set sty data histeps
set bmargin 4
set xlabel 'Energy - keV' font 'Times, 24'
set ylabel 'Focusing factor' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "%(basename)s.eps"
pl "%(basename)s.dat" u 1:2 t "" lw 3
!okular %(basename)s.eps &
""" % {"basename" : basename,})

    def plot_sens(self, append=""):
        basename="%s_sens%s" % (self.name, append)
        file(os.path.join(self.datadir, "%s.plot" % basename), 'w').write("""set sty data histeps
set log
set bmargin 4
set xlabel 'Energy - keV' font 'Times, 24'
set ylabel 'A_{eff} - cm^2' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "%(basename)s.eps"
pl "%(basename)s.dat" u 1:2 t "" lw 3
!okular %(basename)s.eps &
""" % {"basename" : basename,})

    def mc_control_thread(self):
        '''Plots the photons measured on the detector'''
        while True:
            time.sleep(10)


    def plot_mc(self, append=""):
        """Plots the data matrix"""
        basename=self.name+'_mc%s' % append
        file(os.path.join(self.datadir, basename+'.plot'), 'w').write('''set pm3 at s
set size square
set origin 0,-0.13
set size 0.98,1.4

set format z ""
set format cb ""

unset ztics
unset cbtics
unset cbmtics

set xlabel "y - cm"
set ylabel "x - cm"

set ticslevel 0
unset surface

set yrange [*:*]# reverse
set border 0

set view 180,90
set nocolorbox
sp [][][:] '%(name)s.dat' u %(using)s matr w l t ""
pause -1
set font 'Times new Roman, 14'
set ter png enh font 'Times new Roman, 14'
set ou "%(name)s.png"
rep
''' % {'name':basename, 'using':self.detector.using_string()})


    def plot_mcx(self, append=""):
        """Plots the xprojection of the matrix"""
        basename=self.name+'_mc_x%s' % append
        file(os.path.join(self.datadir, basename+'.plot'), 'w').write('''set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 20'
set ter pos enh colo solid 'Times' 24
set ou "%(name)s.eps"
pl [:] \
"%(name)s.dat" u 1:($2/%(normalization)s) t "" lw 3
!okular %(name)s.eps &'''
% {'name':basename, 'normalization':self.detector.cartesian.data.sum()})



    def plot_mcxcenter(self, append=""):
        """Plots the xprojection of the matrix"""
        basename=self.name+'_mc_xcenter%s' % append
        file(os.path.join(self.datadir, basename+'.plot'), 'w').write('''set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 20'
set ter pos enh colo solid 'Times' 24
set ou "%(name)s.eps"
pl [:] \
"%(name)s.dat" u 1:($2/%(normalization)s) t "" lw 3
!okular %(name)s.eps &'''
% {'name':basename, 'normalization':self.detector.cartesian.data.sum()})




    def plot_mcy(self, append=""):
        """Plots the yprojection of the matrix"""
        basename=self.name+'_mc_y%s' % append
        file(os.path.join(self.datadir, basename+'.plot'), 'w').write('''set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "%(name)s.eps"
pl [:] \
"%(name)s.dat" u 1:($2/%(normalization)s) t "" lw 3
!okular %(name)s.eps &''' % {'name':basename, 'normalization':self.detector.cartesian.data.sum()})




    def plot_mcycenter(self, append=""):
        """Plots the yprojection of the matrix"""
        basename=self.name+'_mc_ycenter%s' % append
        file(os.path.join(self.datadir, basename+'.plot'), 'w').write('''set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "%(name)s.eps"
pl [:] \
"%(name)s.dat" u 1:($2/%(normalization)s) t "" lw 3
!okular %(name)s.eps &''' % {'name':basename, 'normalization':self.detector.cartesian.data.sum()})


    


    def plot_mcxy(self, append=""):
        """Plots the x and y projections of the matrix"""
        xname = self.name+'_mc_x%s.dat' % append
        yname = self.name+'_mc_y%s.dat' % append
        basename=self.name+'_mc_xy%s' % append
        file(os.path.join(self.datadir, basename+'.plot'), 'w').write('''set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 24
set ou "%(name)s.eps"
pl [:] \
"%(xname)s" u 1:($2/%(normalization)s) t "x" lw 3,\
"%(yname)s" u 1:($2/%(normalization)s) t "y" lw 3 lt 3
!okular %(name)s.eps &''' % {'xname':xname, 'yname':yname, 'name':basename, 'normalization':self.detector.cartesian.data.sum()})

    def plot_mcr(self, append=""):
        """Plots the x and y projections of the matrix"""
        basename=self.name+'_mc_r%s' % append
        file(os.path.join(self.datadir, basename+'.plot'), 'w').write('''set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Normalized counts' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "%(name)s.eps"
pl [:] \
"%(name)s.dat" u 1:($2/%(normalization)s) t "" lw 3
!okular %(name)s.eps &''' % {'name':basename, 'normalization':self.detector.polar.data.sum()})

    def plot_mcencircled(self, append=""):
        """Plots the x and y projections of the matrix"""
        basename=self.name+'_mc_encircled%s' % append
        file(os.path.join(self.datadir, basename+'.plot'), 'w').write('''set sty data histeps
set xlabel 'cm' font 'Times, 24'
set ylabel 'Encircled photon fraction' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set ou "%(name)s.eps"
pl [:] \
"%(name)s.dat" u 1:($2/%(normalization)s) t "" lw 3
!okular %(name)s.eps &''' % {'name':basename, 'normalization':self.detector.polar.data.sum()})

    def plot_mcs(self, append=""):
        """Prepares for plotting all the MC output"""
        self.plot_mc(append=append)
        self.plot_mcx(append=append)
        self.plot_mcxcenter(append=append)
        self.plot_mcy(append=append)
        self.plot_mcycenter(append=append)
        self.plot_mcr(append=append)
        self.plot_mcxy(append=append)
        self.plot_mcencircled(append=append)

    def plot_sketch(self, append=""):
        """Plots a sketch of the lens using xtalinfo file"""
        names={"plot":self.name+'_sketch%s.plot' % append, "xtalinfo":self.name+'_xtalinfo%s' %append}
        file(os.path.join(self.datadir, names["plot"]), 'w').write('''set data sty dots
set bmargin 4
set size square
set xlabel 'cm' font 'Times, 24'
set ylabel 'cm' font 'Times, 24'
set title '' font 'Times, 24'
set ter pos enh colo solid 'Times' 20
set polar
set ou "%(xtalinfo)s.eps"
pl  \
"%(xtalinfo)s.dat" u 4:3 t "" lw 3
!okular %(xtalinfo)s.eps &''' % names)

# Derived classes
class Ring(Generic):
    """Derived class, redefines only effective_at_keVs."""
    def __init__(self,
                 name = "ring",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 source_lens_distance = 2600.,  # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04, 
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 energy=None,
                 radius=None,
                 dim=[1., 1., 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20,
                 curvature =  1/4000.,   # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 structure = "mosaic",   # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        if radius is not None:
            tB = .5*Physics.atan(radius/Focal)
            d_hkl = Booklet.d_hkl(Z, hkl)
            energy = Physics.BraggEnergy(tB, d_hkl)
        elif energy is None:
            energy=200.


        #if radius is None:
        #    radius = (Booklet.hc * Focal)/(Booklet.d_hkl(Z, hkl)*energy) 


        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="one ring", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                         datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance, source_dimension=source_dimension, detector_distance=detector_distance, Emin=energy, Emax=energy,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature, structure=structure, fwhm=fwhm,
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)

    def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4):
        """Calculates the effective area for a single xtal and then multiplies
        by the number of xtals."""
        self.arealist = zeros(len(self.keVlist), 'f')
        

        if self.curvature == 0.: ###############  for a flat crystal
            if OFFSET==0.:
                self.EA_ring_at_keVs(self.noofxtals, NUM_OF_SIGMA)
            else:  
                for xtal in self.xtals():
                    self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)

        else: #################  for a bent crystal
            if OFFSET==0.:
                self.EA_bendring_at_keVs(self.noofxtals, NUM_OF_SIGMA)
            else:
                for xtal in self.xtals():
                    self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)

##########################################################################################################
##########################################################################################################
##########################################################################################################
class Bendedring(Generic):
    """Derived class, redefines only effective_at_keVs."""
    def __init__(self,
                 name = "ring",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 source_lens_distance = 2600.,  # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04,        # [cm] 
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 energy=None,
                 radius=None,
                 dim=[1., 1., 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20,     # Erased the flag external_curvature=40. It depends upon 
                 curvature =  1/4000.,      # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 external_curvature = 40.,
                 structure = "mosaic",      # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,                   # the Focal, then it will be self.Focal*2 ##### CURVCORRECTION
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        if radius is not None:
            tB = .5*Physics.atan(radius/Focal)
            d_hkl = Booklet.d_hkl(Z, hkl)
            energy = Physics.BraggEnergy(tB, d_hkl)
        elif energy is None:
            energy=200.
        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="one ring", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance,source_dimension=source_dimension, detector_distance=detector_distance, Emin=energy, Emax=energy,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature,
                        external_curvature=external_curvature, structure=structure, fwhm=fwhm,         ######## CURVCORRECTION #####
                        hkl=hkl, Z=Z, microthick=microthick,                      # before fwhm=fwhm 
                        pixels=pixels, pitches=pitches, error=error)              # there was added 
                                                                                  # external_curvature=external_curvature
    def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4): 
        """Calculates the effective area for a single xtal and then multiplies
        by the number of xtals."""
        self.arealist = zeros(len(self.keVlist), 'f')
        # effective calculation
        if OFFSET==0.:
            self.EA_bendring_at_keVs(self.noofxtals, NUM_OF_SIGMA)
        else:
            for xtal in self.xtals():
                self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)                
####################################################################################################################
####################################################################################################################
####################################################################################################################







##########################################################################################################
##########################################################################################################
##########################################################################################################
class Bentmosaicring(Generic):
    """Derived class, redefines only effective_at_keVs."""
    def __init__(self,
                 name = "ring",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 source_lens_distance = 2600., # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04,        # [cm] 
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 energy=None,
                 radius=None,
                 dim=[1., 1., 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20,    # Erased the flag external_curvature=40. It depends upon 
                 curvature =  1/4000.,     # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 external_curvature = 40.,
                 structure = "mosaic",     # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,                  # the Focal, then it will be self.Focal*2 ##### CURVCORRECTION
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        if radius is not None:
            tB = .5*Physics.atan(radius/Focal)
            d_hkl = Booklet.d_hkl(Z, hkl)
            energy = Physics.BraggEnergy(tB, d_hkl)
        elif energy is None:
            energy=200.
        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="one ring", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance, source_dimension = source_dimension, detector_distance=detector_distance, Emin=energy, Emax=energy,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature,
                        external_curvature=external_curvature, structure=structure, fwhm=fwhm,          ######## CURVCORRECTION #####
                        hkl=hkl, Z=Z, microthick=microthick,                      # before fwhm=fwhm 
                        pixels=pixels, pitches=pitches, error=error)              # there was added 
                                                                                  # external_curvature=external_curvature
    def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4): 
        """Calculates the effective area for a single xtal and then multiplies
        by the number of xtals."""
        self.arealist = zeros(len(self.keVlist), 'f')
        # effective calculation
        if OFFSET==0.:
            self.EA_bentmosaic_at_keVs(self.noofxtals, NUM_OF_SIGMA)
        else:
            for xtal in self.xtals():
                self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)                
####################################################################################################################
####################################################################################################################
####################################################################################################################





class SingleXtal(Generic):
    """Derived class, redefines only effective_at_keVs."""
    def __init__(self,
                 name = "one",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 source_lens_distance = 2600., # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04,        # [cm]
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 energy=None,
                 radius=None,
                 dim=[1., 1., 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20,
                 curvature =  1/4000.,   # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        if radius is not None:
            tB = .5*Physics.atan(radius/Focal)
            d_hkl = Booklet.d_hkl(Z, hkl)
            energy = Physics.BraggEnergy(tB, d_hkl)
        elif energy is None:
            energy=100.
        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="spiral", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance, source_dimension = source_dimension, detector_distance=detector_distance, Emin=energy, Emax=energy,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature,
                        structure=structure, fwhm=fwhm, hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)

class Rings(Generic):
    """This will be a class  to deal with systems of rings"""
    def __init__(self,
                 name = "rings",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 source_lens_distance = 2600.,  # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04,         # [cm]
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 Emin=100.,
                 Emax=600.,
                 dim=[1.5, 1.5, 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20,
                 curvature =  1/4000.,   # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="rings", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance, source_dimension = source_dimension, detector_distance=detector_distance, Emin=Emin, Emax=Emax,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature, structure=structure, fwhm=fwhm,
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)
        self.noofrings=1+max([x[0] for x in enumerate(self.rings(COUNT=True))])

    def __str__(self, filename=False):
        """Displays ottica values"""
        dim=(self.dim[0], self.dim[1], self.averagethickness)
        string=("Geometry:                 %s" % self.geo,
                "Profile:                  %s" % self.profile,
                "Thickness profile:        %s" % self.thickness_profile,
                "Thickness threshold:      %s" % repr(self.thickness_threshold),
                "Thickness factor:         %s" % self.thickness_factor,
                "Focal Length:             %s m" % (self.Focal/100),
                "Source lens distance:     %s" % self.source_lens_distance,
                "Source Dimension:         %s" % self.source_dimension,
                "Detector Lens distance:   %s" % self.detector_distance,
                "Inner/Outer radius:       %s/%s m" % (self.rmin/100, self.rmax/100),
                'Z (Atomic number):        %s' % self.Z,
                "Energy range:             [%g, %g] keV" % (self.Emin, self.Emax),
                'hkl:                      (%d,%d,%d)' % tuple(self.hkl),
                "d_hkl:                    %s A" % Booklet.d_hkl(self.Z, self.hkl),
                'Cell volume:              %s A^3' % Booklet.volume[self.Z],
                'Structure factor:         %s' % self.sf,
                'Microcrystal size:        %s micron' % (self.microthick/10000),
                'A_0 constant:             %s keV' % self.A_0const,
                'Material factor constant: %s cm^-1' % self.mat_fac,
                'Crystal tile size:        %gx%gx%g cm^3' % tuple(dim),
                'Number of crystals:       %s' % self.noofxtals,
                'Filling factor:           %s' % self.filling_factor,
                'Number of rings:          %s' % self.noofrings,
                'Crystals volume:          %s cm^3' % self.volume,
                'Crystals weight:          %s kg' % self.weight,
                'mosaicity (FWHM):         %s arcmin' % self.fwhm,
                'Radial Frame width:       %s cm' % self.rad_framewidth,
                'Tangential Frame width:   %s cm' % self.tan_framewidth,
                'Curvature radius of each tile:  %s cm' % self.curvature,)
        return "\n".join(string)

    def rings(self, COUNT=False):
        '''Generator for returning crystals that are placed on concentric rings
        to diffract energies in the nominal range [Emin, Emax]'''
        keV=float(self.Emax)
        while(keV >= self.Emin):
            # yield xtals of ONE ring
            xtal_in_ring = max(x[0] for x in enumerate(self.gen_one_ring(keV=keV, COUNT=COUNT)))+1
            # Next ring: radius and energy calculation (PhD thesis)
            part1 = (self.rho+self.dim[0]/2.+self.rad_framewidth)**2
            part2 = (self.dim[1]/2.+self.tan_framewidth)**2
            sqrt = math.sqrt(part1+part2)
            #####################################################################################
            # NEXTRHO CALCOLATO DA PISA SECONDO ME TROPPO COMPLICATO E SBAGLIATO, CORREGGO SOTTO
            next_rho = self.dim[0]/2. + self.rad_framewidth + sqrt
            #CORREZIONE AL CALCOLO DEL RAGGIO RHO (ENRICO)
            next_rho_enrico = self.rho + self.dim[0] + 2 * self.rad_framewidth
            ####################################################################################
            # next energy
            
            keV = self.radius2keV(next_rho)
        
            yield xtal_in_ring

    def rings_info(self):
        """returns info about rings structure"""
        return [(ring, self.rho, self.radius2keV(self.rho), self.dim[2]) for ring in self.rings()]
        

    def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4):
        """Calculates the effective area for a single xtal and then multiplies
        by the number of xtals for each ring."""
        self.arealist = zeros(len(self.keVlist), 'f')
        # monitoring tools
        countercondition=(self.noofxtals>1000)
        if countercondition:
            counter=Progressbar(update_time=1.,
                               steps=self.noofxtals,
                               title=self.name,
                               name='Effective area',
                               )
        # effective calculation
        if OFFSET==0. and self.error==0.:
            for xtal_in_ring in self.rings():
                self.EA_ring_at_keVs(xtal_in_ring, NUM_OF_SIGMA)
                if countercondition: counter.index+=xtal_in_ring
        else:
            for xtal in self.xtals():
                self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)
                if countercondition: counter.index+=1
        if countercondition: counter.close()
        
######################################################################################
class Bend(Generic):
    """Lens made of bend crystals to be dealed like systems of rings"""
    def __init__(self,
                 name = "rings",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 source_lens_distance = 2600.,  # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04,          # [cm]
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 Emin=100.,
                 Emax=600.,
                 dim=[1.5, 1.5, 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20, #### CURVECORRECTION ### I erased external_curvature because actually
                 curvature =  1/4000.,   # [1/cm] external curvature of the crystals (1/curvature is the curvature radius)
                 external_curvature = 40.,
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,  # mosaic spread [arcmin]         ## depends upon the focal (is self.Focal * 2)
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="rings", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance, source_dimension = source_dimension, detector_distance=detector_distance, Emin=Emin, Emax=Emax,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature,
                        external_curvature=external_curvature, structure=structure, fwhm=fwhm,       ## CURVCORRECTION  there was also external_curvature=external_curvature
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)
        self.noofrings=1+max([x[0] for x in enumerate(self.rings(COUNT=True))])

    def __str__(self, filename=False):
        """Displays ottica values"""
        dim=(self.dim[0], self.dim[1], self.averagethickness)
        fwhm = (self.dim[2] * 180 * 60 ) / (self.external_curvature * 2.6 * 100 * math.pi)  # mosaic spread [arcmin]
        string=("Geometry:                 %s" % self.geo,
                "Profile:                  %s" % self.profile,
                "Thickness profile:        %s" % self.thickness_profile,
                "Source Dimension:         %s" % self.source_dimension,
                "Thickness threshold:      %s" % repr(self.thickness_threshold),
                "Thickness factor:         %s" % self.thickness_factor,
                "Focal Length:             %s m" % (self.Focal/100),
                "Inner/Outer radius:       %s/%s m" % (self.rmin/100, self.rmax/100),
                'Z (Atomic number):        %s' % self.Z,
                "Energy range:             [%g, %g] keV" % (self.Emin, self.Emax),
                'hkl:                      (%d,%d,%d)' % tuple(self.hkl),
                "d_hkl:                    %s A" % Booklet.d_hkl(self.Z, self.hkl),
                'Cell volume:              %s A^3' % Booklet.volume[self.Z],
                'Structure factor:         %s' % self.sf,
                'Microcrystal size:        %s micron' % (self.microthick/10000),
                'A_0 constant:             %s keV' % self.A_0const,
                'Material factor constant: %s cm^-1' % self.mat_fac,
                'Crystal tile size:        %gx%gx%g cm^3' % tuple(dim),
                'Number of crystals:       %s' % self.noofxtals,
                'Filling factor:           %s' % self.filling_factor,
                'Number of rings:          %s' % self.noofrings,
                'Crystals volume:          %s cm^3' % self.volume,
                'Crystals weight:          %s kg' % self.weight,
                'mosaicity (FWHM):         %s arcmin' % fwhm,
                'external curvature:         %s m' % self.external_curvature,
                'Radial Frame width:       %s cm' % self.rad_framewidth,
                'Tangential Frame width:   %s cm' % self.tan_framewidth,
                'Curvature radius of each tile:  %s cm' % self.curvature,)
        return "\n".join(string)

    def rings(self, COUNT=False):
        '''Generator for returning crystals that are placed on concentric rings
        to diffract energies in the nominal range [Emin, Emax]'''
        keV=float(self.Emax)
        while(keV >= self.Emin):
            # yield xtals of ONE ring
            xtal_in_ring = max(x[0] for x in enumerate(self.gen_one_ring(keV=keV, COUNT=COUNT)))+1
            # Next ring: radius and energy calculation (PhD thesis)
            part1 = (self.rho+self.dim[0]/2.+self.rad_framewidth)**2
            part2 = (self.dim[1]/2.+self.tan_framewidth)**2
            sqrt = math.sqrt(part1+part2)
            #####################################################################################
            # NEXTRHO CALCOLATO DA PISA SECONDO ME TROPPO COMPLICATO E SBAGLIATO, CORREGGO SOTTO
            next_rho = self.dim[0]/2. + self.rad_framewidth + sqrt
            #CORREZIONE AL CALCOLO DEL RAGGIO RHO (ENRICO)
            next_rho_enrico = self.rho + self.dim[0] + 2 * self.rad_framewidth
            ####################################################################################
            # next energy
            keV = self.radius2keV(next_rho)
            yield xtal_in_ring

    def rings_info(self):
        """returns info about rings structure"""
        return [(ring, self.rho, self.radius2keV(self.rho), self.dim[2]) for ring in self.rings()]
        

    def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4):
        """Calculates the effective area for a single xtal and then multiplies
        by the number of xtals for each ring."""
        self.arealist = zeros(len(self.keVlist), 'f')
        # monitoring tools
        countercondition=(self.noofxtals>1000)
        if countercondition:
            counter=Progressbar(update_time=1.,
                               steps=self.noofxtals,
                               title=self.name,
                               name='Effective area',
                               )
        # effective calculation
        if OFFSET==0. and self.error==0.:
            for xtal_in_ring in self.rings():
                self.EA_bendring_at_keVs(xtal_in_ring, NUM_OF_SIGMA)
                if countercondition: counter.index+=xtal_in_ring
        else:
            for xtal in self.xtals():
                self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)
                if countercondition: counter.index+=1
        if countercondition: counter.close()
######################################################################################


############################################################################
### BENT MOSAIC 
#############################################################################
class Bentmosaic(Generic):
    """Lens made of bend crystals to be dealed like systems of rings"""
    def __init__(self,
                 name = "rings",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 2000.,
                 Emin=100.,
                 Emax=600.,
                 dim=[1.5, 1.5, 0.4],
                 source_dimension = 0.04,          # [cm]
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20, #### CURVECORRECTION ### I erased external_curvature because actually
                 curvature =  1/4000.,   # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 external_curvature = 40.,
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,  # mosaic spread [arcmin]         ## depends upon the focal (is self.Focal * 2)
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="rings", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance,  source_dimension = source_dimension, detector_distance=detector_distance, Emin=Emin, Emax=Emax,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature,
                        external_curvature=external_curvature, structure=structure, fwhm=fwhm,       ## CURVCORRECTION  there was also external_curvature=external_curvature
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)
        self.noofrings=1+max([x[0] for x in enumerate(self.rings(COUNT=True))])

    def __str__(self, filename=False):
        """Displays ottica values"""
        dim=(self.dim[0], self.dim[1], self.averagethickness)
        fwhm = (self.dim[2] * 180 * 60 ) / (self.external_curvature * 2.6 * 100 * math.pi)  # mosaic spread [arcmin]
        string=("Geometry:                 %s" % self.geo,
                "Profile:                  %s" % self.profile,
                "Thickness profile:        %s" % self.thickness_profile,
                "Thickness threshold:      %s" % repr(self.thickness_threshold),
                "Thickness factor:         %s" % self.thickness_factor,
                "Source Dimension:         %s" % self.source_dimension,
                "Focal Length:             %s m" % (self.Focal/100),
                "Inner/Outer radius:       %s/%s m" % (self.rmin/100, self.rmax/100),
                'Z (Atomic number):        %s' % self.Z,
                "Energy range:             [%g, %g] keV" % (self.Emin, self.Emax),
                'hkl:                      (%d,%d,%d)' % tuple(self.hkl),
                "d_hkl:                    %s A" % Booklet.d_hkl(self.Z, self.hkl),
                'Cell volume:              %s A^3' % Booklet.volume[self.Z],
                'Structure factor:         %s' % self.sf,
                'Microcrystal size:        %s micron' % (self.microthick/10000),
                'A_0 constant:             %s keV' % self.A_0const,
                'Material factor constant: %s cm^-1' % self.mat_fac,
                'Crystal tile size:        %gx%gx%g cm^3' % tuple(dim),
                'Number of crystals:       %s' % self.noofxtals,
                'Filling factor:           %s' % self.filling_factor,
                'Number of rings:          %s' % self.noofrings,
                'Crystals volume:          %s cm^3' % self.volume,
                'Crystals weight:          %s kg' % self.weight,
                'mosaicity (FWHM):         %s arcmin' % fwhm,
                'external curvature:         %s m' % self.external_curvature,
                'Radial Frame width:       %s cm' % self.rad_framewidth,
                'Tangential Frame width:   %s cm' % self.tan_framewidth,
                'Curvature radius of each tile:  %s cm' % self.curvature,)
        return "\n".join(string)

    def rings(self, COUNT=False):
        '''Generator for returning crystals that are placed on concentric rings
        to diffract energies in the nominal range [Emin, Emax]'''
        keV=float(self.Emax)
        while(keV >= self.Emin):
            # yield xtals of ONE ring
            xtal_in_ring = max(x[0] for x in enumerate(self.gen_one_ring(keV=keV, COUNT=COUNT)))+1
            # Next ring: radius and energy calculation (PhD thesis)
            part1 = (self.rho+self.dim[0]/2.+self.rad_framewidth)**2
            part2 = (self.dim[1]/2.+self.tan_framewidth)**2
            sqrt = math.sqrt(part1+part2)
            next_rho = self.dim[0]/2. + self.rad_framewidth + sqrt
            # next energy
            keV = self.radius2keV(next_rho)
            yield xtal_in_ring

    def rings_info(self):
        """returns info about rings structure"""
        return [(ring, self.rho, self.radius2keV(self.rho), self.dim[2]) for ring in self.rings()]
        

    def EA_at_keVs(self, OFFSET=0., NUM_OF_SIGMA=4):
        """Calculates the effective area for a single xtal and then multiplies
        by the number of xtals for each ring."""
        self.arealist = zeros(len(self.keVlist), 'f')
        # monitoring tools
        countercondition=(self.noofxtals>1000)
        if countercondition:
            counter=Progressbar(update_time=1.,
                               steps=self.noofxtals,
                               title=self.name,
                               name='Effective area',
                               )
        # effective calculation
        if OFFSET==0. and self.error==0.:
            for xtal_in_ring in self.rings():
                self.EA_bentmosaic_at_keVs(xtal_in_ring, NUM_OF_SIGMA)  #######################!!!!!!!!!!!!!!!!!!!!!
                if countercondition: counter.index+=xtal_in_ring
        else:
            for xtal in self.xtals():
                self.EA_xtal_at_keVs(OFFSET, NUM_OF_SIGMA)
                if countercondition: counter.index+=1
        if countercondition: counter.close()
######################################################################################

###############��

class Sector(Rings):
    """Class to deal with a circular sector of a ring lens"""
    def __init__(self,
                 name = "sector",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 10000.,
                 Emin=60.,
                 Emax=200.,
                 minphi=0.,
                 maxphi=pi/2.,
                 dim=[1., 1., 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20,
                 curvature =  1/4000.,   # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        """Init of the class. The same of the parent Lens."""
        self.minphi, self.maxphi = minphi, maxphi
        Rings.__init__(self, name=name, profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance, detector_distance=detector_distance, Emin=Emin, Emax=Emax,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature, structure=structure, fwhm=fwhm,
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)

    def __str__(self, filename=False):
        """Displays ottica values"""
        dim=(self.dim[0], self.dim[1], self.averagethickness)
        string=("Geometry:                 %s" % self.geo,
                "Profile:                  %s" % self.profile,
                "Thickness profile:        %s" % self.thickness_profile,
                "Thickness threshold:      %s" % repr(self.thickness_threshold),
                "Thickness factor:         %s" % self.thickness_factor,
                "Focal Length:             %s m" % (self.Focal/100),
                "Inner/Outer radius:       %s/%s m" % (self.rmin/100, self.rmax/100),
                "Angular sector:           [%s,%s]" % (self.minphi/pi, self.maxphi/pi),
                'Z (Atomic number):        %s' % self.Z,
                "Energy range:             [%g, %g] keV" % (self.Emin, self.Emax),
                'hkl:                      (%d,%d,%d)' % tuple(self.hkl),
                "d_hkl:                    %s A" % Booklet.d_hkl(self.Z, self.hkl),
                'Cell volume:              %s A^3' % Booklet.volume[self.Z],
                'Structure factor:         %s' % self.sf,
                'Microcrystal size:        %s micron' % (self.microthick/10000),
                'A_0 constant:             %s keV' % self.A_0const,
                'Material factor constant: %s cm^-1' % self.mat_fac,
                'Crystal tile size:        %gx%gx%g cm^3' % tuple(dim),
                'Number of crystals:       %s' % self.noofxtals,
                'Filling factor:           %s' % self.filling_factor,
                'Number of rings:          %s' % self.noofrings,
                'Crystals volume:          %s cm^3' % self.volume,
                'Crystals weight:          %s kg' % self.weight,
                'mosaicity (FWHM):         %s arcmin' % self.fwhm,
                'Radial Frame width:       %s cm' % self.rad_framewidth,
                'Tangential Frame width:   %s cm' % self.tan_framewidth,
                'Curvature radius of each tile:  %s cm' % self.curvature,)
        return "\n".join(string)

    def gen_one_ring(self, keV=False, xtalnumber=False, COUNT=False):
        '''Generator for returning crystals that are placed on a single ring
        optimized for a defined energy (default Emin+Emax/2)

        Differ from the Base class one because the angular range is different
        '''
        angular_range=(self.maxphi-self.minphi)
        if not keV: keV = 0.5*(self.Emin+self.Emax)
        if not xtalnumber:
            tB=self.keV2Bragg(keV)
            r=2.*self.Focal*math.sin(tB)    ##### CORRECTION FOR FLAT OR SPHERICAL LENS
            #r=self.Focal*math.tan(2*tB)
            num=self.tan_framewidth+self.dim[1]/2
            den=r-self.dim[0]/2
            xtalnumber=int(angular_range/2./math.atan(num/den))
        Delta_phi=angular_range/xtalnumber
        for i in range(xtalnumber):
            self.set_position(keV, self.minphi+Delta_phi*i)
            self.set_thickness(keV)
            yield None

    def rings(self, COUNT=False):
        '''Generator for returning crystals that are placed on a sector of a concentric rings lens
        to diffract energies in the nominal range [Emin, Emax]
        '''
        keV=float(self.Emax)
        while(keV >= self.Emin):
            # yield xtals of ONE ring
            xtal_in_ring = max(x[0] for x in enumerate(self.gen_one_ring(keV=keV, COUNT=True)))+1

            # Next ring: radius and energy calculation
            part1 = (self.rho+self.dim[0]/2.+self.rad_framewidth)**2
            part2 = (self.dim[1]/2.+self.tan_framewidth)**2
            sqrt_=sqrt(part1+part2)
            next_rho = self.dim[0]/2. + self.rad_framewidth + sqrt_
            # next energy
            keV = self.radius2keV(next_rho)
            yield xtal_in_ring

class Petal(Generic):
    """This will be a class  to deal with a petal composed of rings"""
    def __init__(self,
                 name = "petals",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 200.,
                 source_lens_distance = 2600.,  # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04,          # [cm]
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 Emin=60.,
                 Emax=200.,
                 dim=[1., 1., 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20,
                 curvature =  1/4000.,  # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="spiral", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal, source_lens_distance = source_lens_distance, source_dimension = source_dimension, detector_distance=detector_distance, Emin=Emin, Emax=Emax,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature, structure=structure, fwhm=fwhm,
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)

    def gen_spiral_fixed(self, COUNT=False):
        '''Generator for returning crystals that are placed on a spiral and
        with fixed thickness'''
        # kickstart variables
        phi, keV = 0., self.Emax

        # Not so good calculation of the spacing between the tiles
        # TODO: Do it better
        shift_rad = self.dim[0]/2. + self.rad_framewidth
        shift_tan = self.dim[1]/2. + self.tan_framewidth
        spiral_step = self.dim[0] + 2.*self.rad_framewidth

        # START placing the tiles
        while(keV >= self.Emin):
            # modifying crystal
            self.set_position(keV, phi)
            self.set_thickness(keV)

            # Preparing for new step
            phistep= 2.* atan((shift_tan) / (self.rho - shift_rad))
            phi += phistep
            next_rho = self.rmin + spiral_step * phi/(2.*pi)
            keV=self.radius2keV(next_rho)
            if 0<=self.r[1]<=50 and self.r[0]>=0.: yield None

class GRI(Rings):
    """Class for GRI lens"""
    def __init__(self,
                 name = "GRI",
                 thickness_profile = "optimized",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Rmin=88.,
                 Rmax=185.,
                 dim=[1.5, 1.5, 0.4],
                 fraction=1.,
                 rad_framewidth = 0.05,
                 tan_framewidth = 0.05,
                 curvature =  1/4000.,   # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 error=0.,
                 ):
        # Preparing to init
        FOCAL=1e3 # 100 m
        tBmin, tBmax = .5*Physics.atan(Rmin/FOCAL), .5*Physics.atan(Rmax/FOCAL)
        d_hkl = Booklet.d_hkl(Z, hkl)
        Emin, Emax = Physics.BraggEnergy(tBmax, d_hkl), Physics.BraggEnergy(tBmin, d_hkl)
        self.fraction=fraction
        # init
        Rings.__init__(self, name=name, profile="flat",
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=FOCAL, Emin=Emin, Emax=Emax,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature, structure=structure, fwhm=fwhm,
                        hkl=hkl, Z=Z, microthick=microthick,
                        pixels=(256,256), pitches=(0.1,0.1), error=error)

    def gen_one_ring(self, keV=False, COUNT=False):
        '''Modified generator to consider the 8 sectors of 36.4 deg'''
        numofmodules=8
        angular_range=36.4*(pi/180)
        angular_gap=8.6*(pi/180)
        tB=self.keV2Bragg(keV)
        r=2.*self.Focal*math.sin(tB)
        num=self.tan_framewidth+self.dim[1]/2
        den=r-self.dim[0]/2
        xtalinsector=int(angular_range/(2.*math.atan(num/den))*self.fraction)
        xtalnumber=xtalinsector*numofmodules
        Delta_phi=angular_range/xtalinsector
        phi=0.
        self.set_position(keV, phi)
        if False:#COUNT is True:
            for i in range(xtalnumber): yield None
        else:
            for module in range(numofmodules):
                # add space before even sector
                if module%2==0: phi+=angular_gap
                for xtal in range(xtalinsector):
                    phi+=Delta_phi
                    self.set_position(keV, phi)
                    self.set_thickness(keV)
                    yield None
                # add space after odd sector
                if module%2==1: phi+=angular_gap

class Spiral(Generic):
    """This will be a class  to deal with systems of rings"""
    def __init__(self,
                 name = "spiral",
                 profile = "spherical",
                 thickness_profile = "fixed",
                 thickness_threshold = (0, 10),
                 thickness_factor=1.,
                 datadir="data",
                 Focal = 10000.,
                 source_lens_distance = 2600.,  # [cm] or put 'infinitive' for astrophysics sources
                 source_dimension = 0.04,          # [cm]
                 detector_distance = 2000.,     # [cm] In focus: detector_distance = Focal. If not --> lens out of focus
                 Emin=60.,
                 Emax=200.,
                 dim=[1.5, 1.5, 0.4],
                 rad_framewidth = 0.20,
                 tan_framewidth = 0.20,
                 curvature =  1/4000.,   # [1/cm] external curvature of the crystals (1/curvature is the curvature radius) 
                 structure = "mosaic",         # crystalline structure: "mosaic" or "perfect"
                 fwhm=1.,
                 hkl=[1,1,1],
                 Z=29,
                 microthick=0.,
                 pixels=(50,50),
                 pitches=(0.2,0.2),
                 error=0.,
                 ):
        """Init of the class. The same of the parent Lens."""
        Generic.__init__(self, name=name, geo="spiral", profile=profile,
                        thickness_profile=thickness_profile,
                        thickness_threshold=thickness_threshold,
                        thickness_factor=thickness_factor,
                        datadir=datadir, Focal=Focal,source_lens_distance = source_lens_distance,  source_dimension = source_dimension, detector_distance=detector_distance, Emin=Emin, Emax=Emax,
                        dim=dim, rad_framewidth=rad_framewidth, tan_framewidth=tan_framewidth, curvature=curvature,
                        structure=structure, fwhm=fwhm, hkl=hkl, Z=Z, microthick=microthick,
                        pixels=pixels, pitches=pitches, error=error)

if __name__=='__main__':
    for fwhm in (0.4,):#,0.6,0.8,1.0,1.5,):
        b=GRI(Z=32, fwhm=fwhm,error=0.1,Rmax=90.)
        b.macro(INFO=True,
                XTALINFO=True,
                AREA=True,
                MC=True,
                SENS=True,
                S=True,
                PHOTONS=5e3,
                ERANGE=(140., 800.01, 1.),
                OFFSET=0.,
                )
