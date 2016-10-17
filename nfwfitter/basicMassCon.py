#########
# Defines mass-con relations
#########

import numpy as np
import nfwutils
import nfwmodeltools as tools



class constant(object):

    def configure(self, config):

        self.concentration = config['concentration']

    def __call__(self, mass, z, overdensity = 200):

        return self.concentration

#######################################################

class OldDuffy(object):

    def ratioResid(self, ratio, conc, deltaRatio):
        """Helper function for duffyConcentration."""
        return ratio**3 * deltaRatio \
            - (np.log(1. + conc * ratio) \
                   - conc * ratio / (1. + conc * ratio)) \
                   / (np.log(1. + conc) - conc / (1. + conc))


    def __call__(self, m, z, overdensity = 200):

        """Compute the Duffy et al. mass concentration relation for a halo
        with
        m           - in M_sun/h
        z           - redshift
        overdensity - wrt to the critical density

        Borrowed from J. Dietrich, 2013
        """
        n = 1e5
        A = 5.71
        B = -0.084
        C = -0.47

        if overdensity != 200:
        # Get a first estimate of the concentration to convert to M200crit
            c0 = A * (m / 2e12)**B * (1. + z)**C
            delta200c = overdensity / 200.
            minval = 0.2
            maxval = 5.
            ratio = minval + np.arange(n) / n * (maxval - minval)
            res = self.ratioResid(ratio, c0, delta200c)
            rRatio = ratio[(res**2).argmin()]
            mRatio = delta200c * rRatio**3
            m200c = m / mRatio
            # Convert input to M200c using the updated concentration
            nIter = 2
            for i in range(nIter):
                c = A * (m200c / 2e12)**B * (1. + z)**C
                res = self.ratioResid(ratio, c, delta200c)
                rRatio = ratio[(res**2).argmin()]
                mRatio = delta200c * rRatio**3
                m200c = m / mRatio
        else:
            m200c = m
        return A * (m200c / 2e12)**B * (1. + z)**C


#########


class Duffy(object):

    def __call__(self, m, z, overdensity = 200):

        """Compute the Duffy et al. mass concentration relation for a halo
        with
        m           - in M_sun/h
        z           - redshift
        overdensity - wrt to the critical density

        Inspired by J. Dietrich, 2013
        Heavily optimized by D. Applegate
        """
        n = 1e5
        A = 5.71
        B = -0.084
        C = -0.47

        rho_crit = nfwutils.global_cosmology.rho_crit(z)

        c0 = A * (m / 2e12)**B * (1. + z)**C

        if overdensity == 200:
            return c0

        oldC = 100.
        while (np.abs(c0 - oldC)/oldC >= 0.05):
            oldC = c0
            rscale = tools.rscaleConstM(m, c0, rho_crit, overdensity)
            r200 = c0*rscale
            m200c = 200*rho_crit*(4*np.pi/3)*r200**3
            c0 = A * (m200c / 2e12)**B * (1. + z)**C

        return c0


#########


class Bhattacharya(object):

    def __call__(self, m, z):

        """Compute the Duffy et al. mass concentration relation for a halo
        with
        m           - in M_sun/h
        z           - redshift
        overdensity - wrt to the critical density
        '''

        Borrowed from J. Dietrich, 2013
        """
        
        Dz = nfwutils.global_cosmology.GrowthFactor(z)
        nu = (1.12*(m/5e13)**0.3 + 0.53)/Dz

        c200 = 5.9*(Dz**0.54)*(nu**-0.35)

        return c200



