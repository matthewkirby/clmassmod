###############
# Script to compare 2D surface density profiles from MXXL to SD profile predicted by NFW and Diemer
###############

import os
import numpy as np
import matplotlib.pyplot as plt
import cPickle
import pkg_resources



import colossus.cosmology.cosmology as cCosmo
import colossus.halo.concentration as chc
import colossus.halo.profile_dk14 as dk14prof
import colossus.defaults as cDefaults

import nfwfitter.nfwutils as nfwutils
import nfwfitter.colossusMassCon as cmc
import nfwfitter.readMXXLProfile as readMXXLProfile
import nfwfitter.readAnalytic as readAnalytic

########
# mass tables

answers = {
    'halo_41' : cPickle.load(pkg_resources.resource_stream('nfwfitter',
                                                           'data/mxxlsnap41_answers.pkl')),
    'halo_54' : cPickle.load(pkg_resources.resource_stream('nfwfitter',
                                                           'data/mxxlsnap54_answers.pkl'))
}

########

def computeSDProfiles(halobase, mcmodel='diemer15'):

    #read 2D profile
    simprofile = readMXXLProfile.MXXLProfile('{}.radial_profile.txt'.format(halobase))
    simarea = np.pi*(simprofile.outer_radius**2 - simprofile.inner_radius**2)
    simdensity = simprofile.diff_mass / simarea  #M_sol / Mpc**2
    r_mpc = simprofile.median_radius

    #read halo mass
    sim_and_haloid = os.path.basename(halobase)
    tokens = sim_and_haloid.split('_')
    simid = '_'.join(tokens[:2])
    haloid = '_'.join(tokens[2:])

    #make sure cosmology always matches
    curcosmo = readMXXLProfile.cosmo
    nfwutils.global_cosmology.set_cosmology(curcosmo)
    cmc.matchCosmo()

    #compute Diemer SD prediction
    r_kpch = (r_mpc*1000*curcosmo.h)
    
    m200 = answers[simid][haloid]['m200'] #M_sol/h
    zcluster = answers[simid][haloid]['redshift']
    c200 = chc.concentration(m200, '200c', zcluster, model=mcmodel)

    diemer_profile = dk14prof.getDK14ProfileWithOuterTerms(M = m200, c = c200, z = zcluster, mdef = '200c')
    surfacedensity_func, deltaSigma_func = readAnalytic.calcLensingTerms(diemer_profile, np.max(r_kpch))
    convert_units = 1./(curcosmo.h*1e6) #M_sol / Mpc^2 -> diemer units
    diemer_surfacedensity = surfacedensity_func(r_kpch)/convert_units

    return dict(radius=r_mpc,
                simprofile=simdensity,
                diemerprofile=diemer_surfacedensity)
    

########


def plotSDProfile(ax, halobase):

    profiles = computeSDProfiles(halobase)
    r_mpc = profiles['radius']
    simdensity = profiles['simprofile']
    diemer_surfacedensity = profiles['diemerprofile']

    #and plot results

    ax.loglog(r_mpc, simdensity, 'ko', markersize=5)
    ax.loglog(r_mpc, diemer_surfacedensity, 'b-')
    ax.axvline(3., c='k', linewidth=3, linestyle='--')


#######

def plotAvgProfile(ax, halobases, mcmodel='diemer15'):

    
    profiles = [computeSDProfiles(halo, mcmodel=mcmodel) for halo in halobases]


    #average up the profiles
    averadius = np.mean(np.row_stack([x['radius'] for x in profiles]), axis=0)
    avesimprofile = np.mean(np.row_stack([x['simprofile'] for x in profiles]), axis=0)
    avediemerprofile = np.mean(np.row_stack([x['diemerprofile'] for x in profiles]), axis=0)

    ax.loglog(averadius, avesimprofile, 'ko', markersize=5)
    ax.loglog(averadius, avediemerprofile, 'b-')
    ax.axvline(3., c='k', linewidth=3, linestyle='--')
