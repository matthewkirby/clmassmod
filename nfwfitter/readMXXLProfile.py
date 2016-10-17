##################
# Utility to parse Stefan's 2D and 3D mass profiles
##################


import numpy as np
import readtxtfile
import varcontainer
import re, os.path
import nfwutils
import scipy.interpolate as interp

#################

h=.73
#m_p = 7.18*8.60657e8/h  #M_sol
m_p = 8.456e9
cosmo = nfwutils.Cosmology(omega_m = 0.25, omega_l = 0.75, h = 0.73)

################

profile_start_template = re.compile('r_lo')

class MXXLProfile(object):

    def __init__(self, filename):

        self.parseFile(filename)
        self._massEnclosed = None
        self._overdensitymass = None
        self._overdensityradius = None

    ###

    def parseFile(self, filename):

        basedir = os.path.dirname(filename)
        self.redshift = float(readtxtfile.readtxtfile('%s/redshift' % basedir))

        with open(filename) as input:

            lines = input.readlines()

            interior_particles = int(lines[1].split()[1])

            for curline, line in enumerate(lines):
                if profile_start_template.search(line) is not None:
                    break
            assert(curline < len(lines))
            rawprofile = np.array([map(float,x) for x in [x.split() for x in lines[curline+1:]] \
                          if x != []])



        innermass = interior_particles*m_p
        diffMass = rawprofile[:,3]*m_p

        inner_radius = rawprofile[:,0]/h #Mpc
        median_radius = rawprofile[:,1]/h
        outer_radius = rawprofile[:,2]/h
        diff_radius = outer_radius - inner_radius

        self.inner_mass = innermass
        self.inner_radius = inner_radius
        self.median_radius = median_radius
        self.outer_radius = outer_radius
        self.diff_radius = diff_radius
        self.diff_mass = diffMass

    ###

    def massEnclosed(self, radius):
    
        if self._massEnclosed is None:
            mass_enclosed = self.inner_mass + np.cumsum(self.diff_mass)
            self._massEnclosed = interp.interp1d(self.outer_radius, mass_enclosed,
                                                 kind = 'cubic',
                                                 bounds_error = False)

        return self._massEnclosed(radius)

    ###
    
    def overdensityMass(self, delta):

        if self._overdensitymass is None:
            mass_enclosed = self.inner_mass + np.cumsum(self.diff_mass)
            volume_enclosed = (4./3.)*np.pi*self.outer_radius**3
            density = mass_enclosed / volume_enclosed
            rho_crit = cosmo.rho_crit(self.redshift)
            overdensity = density / rho_crit
            self._overdensitymass = interp.interp1d(np.log(overdensity), 
                                                    np.log(mass_enclosed),
                                                    kind = 'linear',
                                                    bounds_error = False)
        return np.exp(self._overdensitymass(np.log(delta)))


    ###

    def overdensityRadius(self, delta):

        if self._overdensityradius is None:
            mass_enclosed = self.inner_mass + np.cumsum(self.diff_mass)
            volume_enclosed = (4./3.)*np.pi*self.outer_radius**3
            density = mass_enclosed / volume_enclosed
            rho_crit = cosmo.rho_crit(self.redshift)
            overdensity = density / rho_crit
            self._overdensityradius = interp.interp1d(np.log(overdensity), 
                                                      np.log(self.outer_radius),
                                                      kind = 'cubic',
                                                      bounds_error = False)
        
        return np.exp(self._overdensityradius(np.log(delta)))
            
