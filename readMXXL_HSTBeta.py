#######################
# Reads in MXXL files for processing
# This version modified betas to match 
#######################

import numpy as np
import binaryutils
import nfwutils
import readMXXL
import scipy.interpolate as interpolate

#######################



class MXXLHSTSimReader(object):

    def __init__(self, config):

        self.config = config

    #########

    def effectiveBeta(self):

        betatable = np.array([[0.4, 0.6],
                              [0.65,  0.43],  
                              [0.70,  0.40 ],
                              [0.75,  0.38 ],  
                              [0.80,  0.36 ],      
                              [0.85,  0.34 ],  
                              [0.90,  0.32 ],  
                              [0.95,  0.30 ],  
                              [1.00,  0.28 ],
                              [1.05,  0.26 ],
                              [1.10,  0.24 ],
                              [1.15,  0.23 ],
                              [1.20,  0.21 ],
                              [1.25,  0.195],
                              [1.30,  0.18 ]])

        model = interpolate.interp1d(betatable[:,0], betatable[:,1])
        effbeta = model(self.config.targetz)

        return effbeta

        


    #########

    def getCosmology(self):

        return nfwutils.Cosmology(omega_m = 0.25, omega_l = 0.75, h=0.73)

    #########

    def load(self, filebase):


        return MXXLHSTSim(filebase, self.effectiveBeta())

    
###############################


class MXXLHSTSim(object):


    #########

    def __init__(self, filebase, effectiveBeta):

        # returns set of angle and mpc distances from the true center (flattened)
        # for each position, provide reduced shear g1 and g2, as well as z and beta for each source

        kappafile = '{0}.convergence_map'.format(filebase)
        gamma1file = '{0}.shear_1_map'.format(filebase)
        gamma2file = '{0}.shear_2_map'.format(filebase)
        answerfile = '{0}.answer'.format(filebase)

        kappa = readMXXL.MXXLBinary(kappafile)
        gamma1 = readMXXL.MXXLBinary(gamma1file)
        gamma2 = readMXXL.MXXLBinary(gamma2file)
        answerfile = asciireader.read(answerfile)

        self.zcluster = kappa.redshift

        delta_mpc, delta_arcmin = kappa.grid()
        delta_mpc = [x.flatten() for x in delta_mpc]
        delta_arcmin = [x.flatten() for x in delta_arcmin]

        r_arcmin = np.sqrt(delta_arcmin[0]**2 + delta_arcmin[1]**2)
        r_mpc = np.sqrt(delta_mpc[0]**2 + delta_mpc[1]**2)


        cosphi = delta_mpc[0] / r_mpc
        sinphi = delta_mpc[1] / r_mpc
    
        sin2phi = 2.0*sinphi*cosphi
        cos2phi = 2.0*cosphi*cosphi-1.0

        self.x_arcmin = delta_arcmin[0]
        self.y_arcmin = delta_arcmin[1]
        self.r_mpc = r_mpc
        self.r_arcmin = r_arcmin
        self.cos2phi = cos2phi
        self.sin2phi = sin2phi
        self.m500 = answerfile['m500c'][0]
        self.m200 = answerfile['m200c'][0]
        self.c200 = answerfile['c200c'][0]


        self.redshifts = 2*np.ones_like(kappa.data).flatten()
        self.beta_s = (effectiveBeta / nfwutils.global_cosmology.beta([1e6], self.zcluster))*np.ones_like(self.redshifts)
        
        betas = effectiveBeta*np.ones_like(self.redshifts)

        self.g1 = betas*gamma1.data.flatten() / (1-betas*kappa.data.flatten())
        self.g2 = betas*gamma2.data.flatten() / (1-betas*kappa.data.flatten())



        
        
    


    
    
    

    
    
    
