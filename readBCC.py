######################
# Reads in BCC files so that nfwfit.py can use them
######################

import h5py, h5pyutils
import nfwutils, clusterTools
import numpy as np


class BCCSimReader(object):

    def getCosmology(self):
        return nfwutils.Cosmology(omega_m = 0.23, omega_l = 0.77, h=.72)

    def load(self, filebase):

        return BCCSim(filebase)

####

class BCCSim(object):

    def __init__(self, filename):


        rawcat = h5py.File(filename, 'r')
        shearcat = h5pyutils.getShearCat(rawcat)

        cluster_ra, cluster_dec, zlens, m200c = h5pyutils.getClusterProperties(rawcat)

        #REMOVE FOREGROUND AND CLUSTER CONTAMINATION (FOR NOW)
        shearcat = shearcat[shearcat['z'] > zlens + 0.1]


        r_arcmin = clusterTools.greatCircleDistance(shearcat['ra'], shearcat['dec'], cluster_ra, cluster_dec)*60

        dL = nfwutils.global_cosmology.angulardist(zlens)

        r_mpc = (r_arcmin/60.)*(np.pi/180.)*dL
    
        phi = clusterTools.positionAngle(shearcat['ra'], shearcat['dec'], cluster_ra, cluster_dec)

        self.cos2phi = np.cos(2*phi)
        self.sin2phi = np.sin(2*phi)

        self.g1 = shearcat['gamma1']
        self.g2 = shearcat['gamma2']

        self.zcluster = zlens
        
        self.r_mpc = r_mpc
        self.r_arcmin = r_arcmin

        self.redshifts = shearcat['z']
        self.beta_s = nfwutils.global_cosmology.beta_s(self.redshifts, self.zcluster)


        

