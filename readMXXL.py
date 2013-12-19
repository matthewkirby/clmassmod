#######################
# Reads in MXXL files for processing
#######################

import numpy as np
import binaryutils
import nfwutils

#######################


class MXXLBinary(object):

    def __init__(self, filename):

        self.filename = filename
        self.parseBinary()

    ####

    def parseBinary(self):

        #The files are binary (little endian). As pseudo-C struct, the file looks like:
        with open(self.filename, 'rb') as input:
            input = open(self.filename, 'rb')

            #double lower_bound[2];                          // lower bound of area represented by plane (in this case in comoving Mpc/h)
            self.lower_bound = binaryutils.readArray(input, 'd', (2,))

            #double upper_bound[2];                          // upper bound of area represented by plane (in this case in comoving Mpc/h)
            self.upper_bound = binaryutils.readArray(input, 'd', (2,))

            #double plane_angle;                             // angle of coordinate axes of plane (should be 0.5 * M_PI)
            self.plane_angle = binaryutils.readVal(input, 'd')

            #double redshift;                                // redfhift of the plane (should be 0.9887)
            self.redshift = binaryutils.readVal(input, 'd')

            #double particle_mass;                           // mass of simulation particle
            self.particle_mass = binaryutils.readVal(input, 'f')

            #int    N_pixels[2];                             // number of pixels in each dimension (should both be 1024 in this case)
            self.npixels = binaryutils.readArray(input, 'i', (2,))

            #float  mass_density[N_pixels[0] * N_pixels[1]]; // surface mass density (in simulation units, i.e. 10^10 M_solar/h per comoving (Mpc/h)^2)
            self.data = binaryutils.readArray(input, 'f', self.npixels)

    def grid(self):

        gridDelta = (self.upper_bound - self.lower_bound)/(self.npixels*nfwutils.global_cosmology.h*(1.+self.redshift))  # grid delta in angular mpc

        X2,X1 = np.meshgrid(np.arange(self.npixels[1]), np.arange(self.npixels[0]))
        
        deltaX1_mpc = self.lower_bound[0] + (X1 + 0.5)*gridDelta[0]
        deltaX2_mpc = self.lower_bound[1] + (X2 + 0.5)*gridDelta[1]

        dL = nfwutils.global_cosmology.angulardist(self.redshift)

        deltaX1_arcmin = (deltaX1_mpc / dL)*(180.*60/np.pi)
        deltaX2_arcmin = (deltaX2_mpc / dL)*(180.*60/np.pi)

        return (deltaX1_mpc, deltaX2_mpc), (deltaX1_arcmin, deltaX2_arcmin)
        
        

        


    

#################################

class MXXLSimReader(object):

    #########

    def getCosmology(self):

        return nfwutils.Cosmology(omega_m = 0.25, omega_l = 0.75, h=0.73)

    #########

    def load(self, filebase):

        return MXXLSim(filebase)

    
###############################


class MXXLSim(object):


    #########

    def __init__(self, filebase):

        # returns set of angle and mpc distances from the true center (flattened)
        # for each position, provide reduced shear g1 and g2, as well as z and beta for each source

        kappafile = '{0}.convergence_map'.format(filebase)
        gamma1file = '{0}.shear_1_map'.format(filebase)
        gamma2file = '{0}.shear_2_map'.format(filebase)

        kappa = MXXLBinary(kappafile)
        gamma1 = MXXLBinary(gamma1file)
        gamma2 = MXXLBinary(gamma2file)

        self.zcluster = kappa.redshift

        self.delta_mpc, self.delta_arcmin = kappa.grid()
        self.delta_mpc = [x.flatten() for x in self.delta_mpc]
        self.delta_arcmin = [x.flatten() for x in self.delta_arcmin]

        self.redshifts = 2*np.ones_like(kappa.data).flatten()
        self.beta_s = nfwutils.global_cosmology.beta_s([2.], self.zcluster)*np.ones_like(self.redshifts)
        
        betas = nfwutils.global_cosmology.beta([2.], self.zcluster)*np.ones_like(self.redshifts)

        self.g1 = betas*gamma1.data.flatten() / (1-betas*kappa.data.flatten())
        self.g2 = betas*gamma2.data.flatten() / (1-betas*kappa.data.flatten())



        
        
    


    
    
    

    
    
    
