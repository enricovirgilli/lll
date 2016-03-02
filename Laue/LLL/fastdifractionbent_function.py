     def fast_diff_bent(self, Focal): 
        print "now diffraction:"
        T = self.xtal.dim[2]
        #order, sigma_ = 0,0.
        # Estimate which is the most important diffrction order 
        #R = math.sqrt(self.photon.r[0]**2 + self.photon.r[1]**2)
        factor = 2.5     
        
        #raggio = math.sqrt (self.photon.r[0]**2 + self.photon.r[1]**2)
        #delta_angle =  (self.xtal.rho - raggio ) * self.xtal.curvature
        #g_theta = self.xtal.gtheta - delta_angle        
        #new_g = Physics.spherical2cartesian(norm(self.xtal.g),self.xtal.gphi, g_theta)
        #print new_g 
        new_g = self.change_the_g(self.xtal.g)
 
        order, sigma_ = 0,0.
        i=0.
        
        if self.xtal.structure == "mosaic":

            for sigma in self.sigmas("gaussian"):
                i+=.5
                if sigma > sigma_: order, sigma_ = int(i)+1, sigma
            internal_spread =  self.xtal.fwhm   
            self.photon.k =  self.photon.k + new_g * order  

        else:

            for sigma in self.sigmas("hat"):
                i+=.5
                if sigma > sigma_: order, sigma_ = int(i)+1, sigma
            internal_spread =  (T * self.xtal.curvature) / factor
            self.photon.k =  self.photon.k + new_g * order 

            #e_window=[]
            #e_window.append(self.xtal.Bragg2keV(t_b_centroid + internal_spread/2, 1)) 
            #e_window.append(self.xtal.Bragg2keV(t_b_centroid - internal_spread/2, 1))
            #if self.photon.e > e_window[0] and self.photon.e < e_window[1]:
            #    self.photon.k =  self.photon.k + new_g * order 
            #else:
            #    self.photon.k = self.photon.k

        tB = self.xtal.keV2Bragg(self.photon.e,1)
        mu = Booklet.mu(self.xtal.Z, self.photon.e)
        up_lambda = math.pi * self.xtal.volume * math.cos(tB)
        down_lambda =  (Booklet.hc /  self.photon.e) * Booklet.rem_angstrom * abs(self.xtal.sf) * ( 1 + math.cos(tB))
        lambda0 = up_lambda / down_lambda
        absorption = math.exp(-(mu * T) / (math.cos(tB) ))
        diffraction = 1 - math.exp(-(math.pi**2 * self.xtal.d_hkl ) / (1e-8 * self.xtal.curvature * lambda0**2) )
        if self.xtal.curvature != 0 and self.xtal.structure == "perfect":
            PoS = diffraction * absorption * sigma_
            self.photon.i = self.photon.i * PoS
        else:
           
            PoS = 0.5 * (1.-math.exp(-sigma_*T)) * math.exp(-Booklet.mu(self.xtal.Z, self.photon.e)*T)            
            self.photon.i = self.photon.i * PoS
        
