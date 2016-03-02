from Xtal import Physics
Booklet=Physics
math=Physics
pi=math.pi
from random import uniform, randint
from numpy import zeros, array
from Xtal import Source
from Xtal import Detector
from Xtal import Xtal
Photon=Source.Photon
PhoXtal=Source.PhoXtal

ANGLES=[0., pi/4, pi/2, pi]

#def circular_collimator(p, r=1., x=0., y=0.):
    #"""Blocks the photon (p) outside an hole of radius r and center (x,y)"""
    #if (p.r[0]-x)**2+(p.r[1]-y)**2>r**2: p.i=0.

#def horizontal_slit(p, w=1., y=0.):
    #"""Simulate an horizontal slit of width w at an height y"""
    #if abs(p.r[1]-y)>w/2.: p.i=0.

#def vertical_slit(p, w=1., x=0.):
    #"""Simulate an horizontal slit of width w aranslated by x"""
    #if abs(p.r[0]-x)>w/2.: p.i=0.

#def wire(p, w=.3, m=1., q=0., x=0., y=0.):
    #"""Simulate an absorbing wire described by the equation y=mx+q. Vertical wires have m equal to '+inf'"""
    #if m=="inf":
        #if abs(p.r[0]-x)<w/2.:p.i=0
        #return
    #d=abs(m*(p.r[0]-x)+q-(p.r[1]-y))/math.sqrt(1+m**2)
    #if d<w/2: p.i=0.

def XRT(xtal, RSPOT=0.0150):
    """Yields photons produced by the a tube with """
    RSPOTSQUARE=RSPOT**2
    while True:
        # Random coordinates of the photon production point
        r = array((uniform(-RSPOT,RSPOT), uniform(-RSPOT,RSPOT), 0.))
        if r[0]**2+r[1]**2<RSPOTSQUARE:
            # Random coordinates of the photon destination point
            R = array((xtal.r[0]+uniform(-xtal.dim[0]/2, xtal.dim[0]/2), xtal.r[1]+uniform(-xtal.dim[1]/2, xtal.dim[1]/2), 600.))
            e= uniform(70., 110.)
            return Photon(r=R, k=(R-r), e=e)

def XTALS():
    # Set crytals positions
    LENS_RADIUS=17.8
    xtals=[Xtal.Xtal(r = (LENS_RADIUS*math.cos(ANGLES[i]), LENS_RADIUS*math.sin(ANGLES[i]), 600.),
                    dim=(1.5,1.5,0.2),
                    Z=29,
                    hkl=(1,1,1),
                    fwhm=4., # [arcmin]
                    microthick=0.) for i in range(4)]
    # Define orientation error
    error_arcmin=10.#10./60.
    error=error_arcmin/60./180.*pi
    # Set crystals orientation
    gnorm=2.*pi/xtals[0].d_hkl
    for i in range(len(xtals)):
        # Computing angular coordinates
        theta_error=uniform(-error, error)
        phi_error=uniform(-error, error)
        theta=pi/2.+error#theta_error
        phi=ANGLES[i]+phi_error-pi
        g=Physics.spherical2cartesian(gnorm, phi, theta)
        xtals[i].g=g
    return xtals

def macro(d):
    xtals=XTALS()
    pxs=[PhoXtal(None,x) for x in xtals]
    try:
        while True:
            for px in pxs:
                for j in range(1):
                    px.photon=XRT(px.xtal)
                    px.fast_diffraction()
                    px.photon.travelz(1100.)
                    d.measure(px.photon)
    except KeyboardInterrupt:
        pass
    d.plot()

if __name__=='__main__':
    d=Detector.PSD_cartesian(pixels=(128,128), pitches=(.1,.1))
    macro(d)
