#!/usr/bin/env python

import nfwutils, nfwnoise, nfwfit
import numpy as np
import scipy.interpolate
import scipy.integrate
import shelve
from multiprocessing import Pool as ThreadPool

#####

def doMassUncert(argdict):

    return nfwnoise.calcMassUncert(**argdict)


#####

if __name__ == '__main__':


    nthreads = 4
    pool = ThreadPool(nthreads)

    nmasses = 10000
    nperthread = nmasses/nthreads



    config = nfwfit.readConfiguration('mxxl_imperial/snap41/c4_r10/config.sh')
    r_mpc = np.arange(0.75, 3.0, 0.2)
    beta_s = 0.5
    zcluster = 0.5


    masses = 10**np.random.uniform(14, 16, nmasses)
    concens = np.random.uniform(1.1, 19.9, nmasses)
    testmasses = np.logspace(14., 16., 400)
    log10testmasses = np.log10(testmasses)

    shearprofiles_p, shearerr_p = nfwnoise.createClusterSet(config, masses, zcluster, r_mpc, beta_s, 0.01)


    profile_threadgroups = [dict(config = config,
                                 r_mpc = r_mpc,
                                 beta_s = beta_s, 
                                 zcluster = zcluster, 
                                 shearprofiles = shearprofiles_p[i:i+nperthread], 
                                 shearerr = shearerr_p, 
                                 testmasses = testmasses) for i in range(0, nmasses, nperthread)]

    results = pool.map(doMassUncert, profile_threadgroups)

    masserrs_p = np.vstack([x[0] for x in results])
    massprobs_p = np.vstack([x[1] for x in results])

    print masserrs_p.shape
    print massprobs_p.shape

    cum_massprobs_p = np.vstack([scipy.integrate.cumtrapz(massprobs_p[i], log10testmasses, initial=0.) for i in range(nmasses)])

    for i in range(nmasses):
        cum_massprobs_p[i,:] = cum_massprobs_p[i,:] / cum_massprobs_p[i,-1]

    cdfs = np.hstack([np.interp([np.log10(masses[i])], log10testmasses, cum_massprobs_p[i]) for i in range(nmasses)])

    output = shelve.open('cdfcheck.r400.shelve', protocol=-1)

    output['config'] = config
    output['r_mpc'] = r_mpc
    output['beta_s'] = beta_s
    output['zcluster'] = zcluster
    output['masses'] = masses
    output['concens'] = concens
    output['testmasses'] = testmasses
    output['shearprofiles'] = shearprofiles_p
    output['shearerr'] = shearerr_p
    output['massprobs'] = massprobs_p
    output['cum_massprobs'] = cum_massprobs_p
    output['cdfs'] = cdfs

    output.close()
