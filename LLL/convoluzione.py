
def largegauss(x,eta,fwhm,NUM_OF_SIGMA=3.):
    Dthetamax = NUM_OF_SIGMA*fwhm
    g = []
    h = []
    #r=Xtal.reflectivity(100., 0., 1)
    for i in range(-400,400):
        i=i*1.0
        #g.append(gaussian(i, eta, x0=0))
        g.append(self.reflectivity(i, 0., order))
        h.append(Phisics.hat(i, fwhm, x0=0))
    hatg = convolve(h,g)
    sl = slice(800-x, 800+x, 1)
    return hatg[sl]#*fwhm
