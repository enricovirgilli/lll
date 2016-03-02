
def EA_ring_at_keVs(self, xtal_in_ring=1., NUM_OF_SIGMA=4):
    order = 1
    tB = math.pi/2-self.gtheta
    ring_cross_section=self.cross_section * xtal_in_ring        
    Emin, Emax = self.Erange(order, OFFSET=0., NUM_OF_SIGMA=NUM_OF_SIGMA)    
    keVs=[(i,keV) for i,keV in enumerate(self.keVlist) if Emin<=keV<=Emax]    
    poss=[(i, self.reflectivity(keV, 0., order)) for i, keV in keVs]    
    for i, pos in poss:
        if pos!=0.: self.arealist[i] = self.arealist[i] + pos * ring_cross_section

def EA_bendring_at_keVs(self, xtal_in_ring=1., NUM_OF_SIGMA=4):
    order = 1
    tB = math.pi/2-self.gtheta
    ring_cross_section=self.cross_section * xtal_in_ring
    alpha=self.dim[0]/(self.external_curvature * 100.)
    internal_spread = self.dim[2] / (2.59 * self.external_curvature * 100. )
    Emin = self.Bragg2keV(tB + alpha/2)
    Emax = self.Bragg2keV(tB - alpha/2)
    keVs=[(i,keV) for i,keV in enumerate(self.keVlist) if Emin<=keV<=Emax]
    poss=[(i, self.peakreflecurved(keV, 0., order)) for i, keV in keVs]
    for i, pos in poss:
        if pos!=0.: self.arealist[i] = self.arealist[i] + pos * internal_spread * (self.external_curvature * 100) * self.dim[1] * xtal_in_ring 



def EA_bentmosaic_at_keVs(self, xtal_in_ring=1., NUM_OF_SIGMA=4):
        tB = math.pi/2-self.gtheta
        ring_cross_section=self.cross_section * xtal_in_ring    
        alpha=self.dim[0]/(self.external_curvature * 100.) 
        internal_spread = self.fwhm
        Emin = self.Bragg2keV(tB + alpha/2)
        Emax = self.Bragg2keV(tB - alpha/2)
        keVs=[(i,keV) for i,keV in enumerate(self.keVlist) if Emin<=keV<=Emax]
        poss=[(i, self.reflectivity(keV, 0., order)) for i, keV in keVs]
        th0=0
        for i, pos in poss:
            if pos!=0.: self.arealist[i] = self.arealist[i] + pos #* internal_spread * (self.external_curvature * 100) * self.dim[1] * xtal_in_ring

