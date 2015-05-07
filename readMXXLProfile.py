##################
# Utility to parse Stefan's 2D and 3D mass profiles
##################


import numpy as np
import readtxtfile
import varcontainer
import re

#################

h=.73
m_p = 8.60657e8/h  #M_sol

################

profile_start_template = re.compile('r_lo')

def parseFile(filename):

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

    profile = varcontainer.VarContainer()
    profile.inner_mass = innermass
    profile.inner_radius = inner_radius
    profile.median_radius = median_radius
    profile.outer_radius = outer_radius
    profile.diff_radius = diff_radius
    profile.diff_mass = diffMass

    return profile

#############
    

    
    
