#######################
# Reads in MXXL files for processing
#######################

import numpy as np
import binaryutils
import nfwutils
import astropy.io.ascii as asciireader

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
            #NOTE: For "mass_map", units of Mpc/h
            #NOTE: For convergence_map, units of arcsec

            self.lower_bound = binaryutils.readArray(input, 'd', (2,))

            #double upper_bound[2];                          // upper bound of area represented by plane (in this case in comoving Mpc/h)

            #NOTE: For "mass_map", units of Mpc/h
            #NOTE: For convergence_map, units of arcsec
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

        gridDelta = (self.upper_bound - self.lower_bound)/(self.npixels)  # grid delta in arcsec

        X2,X1 = np.meshgrid(np.arange(self.npixels[1]), np.arange(self.npixels[0]))
        
        deltaX1_arcmin = (self.lower_bound[0] + (X1 + 0.5)*gridDelta[0])/60.
        deltaX2_arcmin = (self.lower_bound[1] + (X2 + 0.5)*gridDelta[1])/60.

        dL = nfwutils.global_cosmology.angulardist(self.redshift)

        deltaX1_mpc = (deltaX1_arcmin * dL * np.pi)/(180.*60)
        deltaX2_mpc = (deltaX2_arcmin * dL * np.pi)/(180.*60)
        return (deltaX1_mpc, deltaX2_mpc), (deltaX1_arcmin, deltaX2_arcmin)
        
        

        


    

#################################

class MXXLSimReader(object):

    def __init__(self, *args, **keywords):

        pass

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
        answerfile = '{0}.answer'.format(filebase)

        kappa = MXXLBinary(kappafile)
        gamma1 = MXXLBinary(gamma1file)
        gamma2 = MXXLBinary(gamma2file)
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
        self.beta_s = nfwutils.global_cosmology.beta_s([2.], self.zcluster)*np.ones_like(self.redshifts)
        
        betas = nfwutils.global_cosmology.beta([2.], self.zcluster)*np.ones_like(self.redshifts)

        self.g1 = betas*gamma1.data.flatten() / (1-betas*kappa.data.flatten())
        self.g2 = betas*gamma2.data.flatten() / (1-betas*kappa.data.flatten())



        
        
    


    
    
    

    
    
    
