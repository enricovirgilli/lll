def reflectivity(self, keV, OFFSET=0., order=1):
    extinction_factor = self.extinction_factor(keV, self.keV2Bragg(keV), order)
    T = self.dim[2]/math.cos(pi/2-self.gtheta)
    muT=self.mu(keV) * T
    sigmaT = abs(self.sigma(keV, OFFSET, order) * T * extinction_factor)
    r = 0.5 * (1 - math.exp(-2 * sigmaT)) * math.exp(-muT)
    return r

def reflectivity_of_bent_perfect(self, keV, OFFSET=0., order=1):
        internal_curvature = 2.56 * self.external_curvature
        curvature_s = 1/internal_curvature
        real_fwhm = (self.dim[2] * 180 * 60 ) / (internal_curvature * 100 * math.pi)
        Z=self.Z
        tB = self.keV2Bragg(keV)
        polarize = (1 + abs(math.cos(2 * tB )))
        up_for_lambda = math.pi * volumes[Z] * math.cos(tB)
        down_for_lambda =  (Booklet.hc / keV) * Booklet.rem_angstrom * abs(self.sf) #*polarize
        lambda0 = up_for_lambda / down_for_lambda
        T = self.dim[2] * 100000000  ### OK
        cp = curvature_s/10000000000.0
        omega = cp*T
        diffraction = 1 - math.exp(-(math.pi**2 * self.d_hkls[order] ) / (cp * lambda0**2) )  ##### OK
        alpha=self.dim[0]/(self.external_curvature * 100.)
        DeltaTheta = tB - th0
        weight=Physics.hat(DeltaTheta, alpha)        
        return diffraction * self.absorption_in_curved(keV, order) * weight

def reflectivity_of_bent_mosaic(self, keV, OFFSET=0., order=1):
        internal_curvature = 2.56 * self.external_curvature
        curvature_s = 1/internal_curvature
        real_fwhm = (self.dim[2] * 180 * 60 ) / (internal_curvature * 100 * math.pi)
        Z=self.Z
        tB = self.keV2Bragg(keV)
        polarize = (1 + abs(math.cos(2 * tB )))
        up_for_lambda = math.pi * volumes[Z] * math.cos(tB)
        down_for_lambda =  (Booklet.hc / keV) * Booklet.rem_angstrom * abs(self.sf) #*polarize
        lambda0 = up_for_lambda / down_for_lambda
        T = self.dim[2] * 100000000  ### OK
        cp = curvature_s/10000000000.0
        omega = cp*T
        diffraction = 1 - math.exp(-(math.pi**2 * self.d_hkls[order] ) / (cp * lambda0**2) )  ##### OK
        alpha=self.dim[0]/(self.external_curvature * 100.)
        DeltaTheta = tB - th0
        weight=Physics.hat(DeltaTheta, alpha)        
        return diffraction * self.absorption_in_curved(keV, order) * weight
