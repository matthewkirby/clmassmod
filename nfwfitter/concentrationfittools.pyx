################################
# Likelihood for lognormal scatter with normal noise
################################

# Compiling info: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /u/ki/dapple/include/python2.7/ -I /u/ki/dapple/lib/python2.7/site-packages/numpy/core/include/ -o deconvolvedlognormtools.so deconvolvedlognormtools.c

#aifa: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /users/dapple/anaconda/pkgs/numpy-1.7.1-py27_2/lib/python2.7/site-packages/numpy/core/include/ -I /users/dapple/anaconda/include/python2.7/ -o deconvolvedlognormtools.so deconvolvedlognormtools.c


########################

# cython: profile=False

import numpy as np
cimport numpy as np
cimport cython

import scipy.stats


cdef extern from "math.h":
    double exp(double)
    double log(double)
    double sqrt(double)

sqrt2pi = sqrt(2*np.pi)
twopi = 2*np.pi




###############################

@cython.boundscheck(False)
@cython.wraparound(False)
def altintegral_mc(np.ndarray[np.double_t, ndim=1, mode='c'] ml_ints,
                np.ndarray[np.double_t, ndim=1, mode='c'] delta_logmls,
              double logmu, 
                double sigma,
                np.ndarray[np.double_t, ndim=1, mode='c'] cl_ints,
              double logmu_c, 
                double sigma_c,
):

    cdef Py_ssize_t i, nsamples
    nsamples = ml_ints.shape[0]

    cdef double thesum, lognormpart
    thesum = 0.

    cdef double neg2sigma2, sigmasqrt2pi
    neg2sigma2 = -2*(sigma**2)
    sigmasqrt2pi = sigma*sqrt2pi

    cdef double neg2sigma2_c, sigmasqrt2pi_c
    neg2sigma2_c = -2*(sigma_c**2)
    sigmasqrt2pi_c = sigma_c*sqrt2pi


    
    for i from nsamples > i >= 0:

        lognormpart = exp((delta_logmls[i]-logmu)**2/neg2sigma2)/(sigmasqrt2pi*ml_ints[i])
        lognormpart_c = exp((np.log(cl_ints[i])-logmu_c)**2/neg2sigma2_c)/(sigmasqrt2pi_c*cl_ints[i])

        thesum += lognormpart * lognormpart_c

    thesum = thesum / nsamples

    return thesum

#########

@cython.boundscheck(False)
@cython.wraparound(False)
def mcmcloglinearlike_mc(np.ndarray[np.double_t, ndim=2, mode='c'] ml_ints,
                      np.ndarray[np.double_t, ndim=2, mode='c'] delta_logmls,
                      np.ndarray[np.int_t, ndim=1, mode='c'] ngoodsamples,
                      double logmu, 
                      double sigma,
                      np.ndarray[np.double_t, ndim=2, mode='c'] cl_ints,
                      double logmu_c, 
                      double sigma_c,
):


    cdef Py_ssize_t i, nclusters, nsamples
    nclusters = ml_ints.shape[0]

    cdef double sumlogprob = 0.
    cdef double prob = 0.

    

    for i from nclusters > i >= 0:

        nsamples = ngoodsamples[i]
        
        prob = altintegral_mc(ml_ints[i,:nsamples],
                           delta_logmls[i,:nsamples],
                           logmu,
                           sigma,
                           cl_ints[i,:nsamples],
                           logmu_c,
                           sigma_c,
        )



        sumlogprob += log(prob)

    return sumlogprob


##############
