################################
# Likelihood for lognormal scatter with normal noise
################################

# Compiling info: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /u/ki/dapple/include/python2.7/ -I /u/ki/dapple/lib/python2.7/site-packages/numpy/core/include/ -o deconvolvedlognormtools.so deconvolvedlognormtools.c


########################

# cython: profile=False

import numpy as np
cimport numpy as np
cimport cython


cdef extern from "math.h":
    double exp(double)
    double log(double)
    double sqrt(double)

sqrt2pi = sqrt(2*np.pi)
twopi = 2*np.pi




###############################


@cython.boundscheck(False)
@cython.wraparound(False)
def integral(double mlens,
              double merr, 
              double mtrue,
              double logmu, 
              double sigma):

    cdef Py_ssize_t i, nsamples
    nsamples = 50

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] randomdeviates = np.random.standard_normal(nsamples)

    cdef double thesum, logmtrue, ml_int, normpart
    logmtrue = log(mtrue)
    thesum = 0.

    for i from nsamples > i >= 0:

        ml_int = exp(logmu + logmtrue + sigma*randomdeviates[i])

        normpart = exp(-0.5*(ml_int-mlens)**2/merr**2)/(sqrt2pi*merr)


        thesum += normpart

    thesum = thesum / nsamples

    return thesum

@cython.boundscheck(False)
@cython.wraparound(False)
def altintegral(np.ndarray[np.double_t, ndim=1, mode='c'] ml_ints,
                np.ndarray[np.double_t, ndim=1, mode='c'] delta_logmls,
              double logmu, 
              double sigma):

    cdef Py_ssize_t i, nsamples
    nsamples = ml_ints.shape[0]

    cdef double thesum, lognormpart
    thesum = 0.

    cdef double sigma2, sigmasqrt2pi
    neg2sigma2 = -2*(sigma**2)
    sigmasqrt2pi = sigma*sqrt2pi

    for i from nsamples > i >= 0:

        lognormpart = exp((delta_logmls[i]-logmu)**2/neg2sigma2)/(sigmasqrt2pi*ml_ints[i])

        thesum += lognormpart

    thesum = thesum / nsamples

    return thesum

    


@cython.boundscheck(False)
@cython.wraparound(False)
def loglinearlike(np.ndarray[np.double_t, ndim=2, mode='c'] ml_ints,
                  np.ndarray[np.double_t, ndim=2, mode='c'] delta_logmls,
                  double logmu, 
                  double sigma):


    cdef Py_ssize_t i, nclusters
    nclusters = ml_ints.shape[0]

    cdef double sumlogprob = 0.
    cdef double prob = 0.
    

    for i from nclusters > i >= 0:


        
        prob = altintegral(ml_ints[i,:],
                           delta_logmls[i,:],
                           logmu,
                           sigma)



        sumlogprob += log(prob)

    return sumlogprob

        
