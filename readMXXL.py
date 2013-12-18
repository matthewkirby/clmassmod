#######################
# Reads in MXXL files for processing
#######################

import numpy as np
import binaryutils

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


    

#################################
