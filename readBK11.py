######################
# Reads in BK11 files so that nfwfit.py can use them
######################

import pyfits, ldac


class BK11SimReader(object):

    def getCosmology(self):
        return nfwutils.Cosmology(omega_m = 0.27, omega_l = 0.73, h = 0.7)

    def load(self, filebase):

        return BK11Sim(filebase)


class BK11Sim(object):

    def __init__(self, filename):


        sim = ldac.LDACCat(pyfits.open(catalogname)[1])

        boxlength = sim['BOXWIDTHCOMOVINGHINVMPC'] / nfwutils.global_cosmology.h
        gridsize = np.sqrt(float(sim['A00'].shape[1]))

        self.zcluster = float(sim['ZLENS'])

        ####becker angles

        Ng = 512

        dc_readout = nfwutils.global_cosmology.comovingdist(clusterz)+sim['BOXLENGTHCOMOVINGHINVMPC'][0]/2.0
        da_arcmin = np.arctan2(sim['BOXWIDTHCOMOVINGHINVMPC'][0],dc_readout)/np.pi*180.0*60.0/Ng
        da_mpc = (da_arcmin/60.)*(np.pi/(180.))*nfwutils.global_cosmology.angulardist(clusterz)
        

        xcen = Ng/2.0
        ycen = Ng/2.0

        radii_arcmin1 = np.zeros((Ng, Ng))
        radii_arcmin2 = np.zeros(Ng,Ng))
        radii_mpc1 = np.zeros((Ng,Ng))
        radii_mpc2 = np.zeros((Ng,Ng))

        cosphi = np.zeros((Ng,Ng))
        sinphi = np.zeros((Ng,Ng))

        for i in range(Ng):
            for j in range(Ng):
                radii_arcmin1 = (i + 0.5 - xcen)*da_arcmin
                radii_arcmin2 = (j + 0.5 - ycen)*da_arcmin
                radii_mpc1 = radii_arcmin1*(da_mpc / da_arcmin)
                radii_mpc2 = radii_arcmin2*(da_mpc / da_arcmin)



        A00 = sim['A00'].reshape(gridsize, -1)
        A11 = sim['A11'].reshape(gridsize, -1)
        A10 = sim['A10'].reshape(gridsize, -1)
        A01 = sim['A01'].reshape(gridsize, -1)

        kappa = 0.5*(2- A00 - A11)
        gamma1 = 0.5*(A11 - A00)
        gamma2 = -0.5*(A01+A10)
        
        self.g1 = (gamma1/(1-kappa)).flatten()
        self.g2 = (gamma2/(1-kappa)).flatten()

        
        self.delta_mpc = [radii_mpc1.flatten(), radii_mpc2.flatten()]
        self.delta_arcmin = [radii_arcmin1.flatten(), radii_arcmin2.flatten()]

        self.redshifts = np.ones(len(self.g1))*self.zcluster
        self.beta_s = nfwutils.global_cosmology.beta_s(self.redshifts, self.zcluster)


        

