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
def altintegral(double mlens,
              double merr, 
              double mtrue,
              double logmu, 
              double sigma):

    cdef Py_ssize_t i, nsamples
    nsamples = 50

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] randomdeviates = np.random.standard_normal(nsamples)

    cdef double thesum, ml_int, lognormpart, logmtrue
    logmtrue = log(mtrue)
    thesum = 0.

    for i from nsamples > i >= 0:

        ml_int = mlens + merr*randomdeviates[i]

        lognormpart = exp(-0.5*(log(ml_int)-logmtrue-logmu)**2/sigma**2)/(sqrt2pi*sigma*ml_int)

        thesum += lognormpart

    thesum = thesum / nsamples

    return thesum

    


@cython.boundscheck(False)
@cython.wraparound(False)
def loglinearlike(np.ndarray[np.double_t, ndim=1, mode='c'] mlens, 
                  np.ndarray[np.double_t, ndim=1, mode='c'] merr, 
                  np.ndarray[np.double_t, ndim=1, mode='c'] mtrue, 
                  double logmu, 
                  double sigma):


    cdef Py_ssize_t i, nclusters
    nclusters = mlens.shape[0]

    cdef double sumlogprob = 0.
    cdef double prob = 0.
    

    for i from nclusters > i >= 0:


        
        prob = integral(mlens[i],
                        merr[i],
                        mtrue[i],
                        logmu,
                        sigma)

        if prob == 0.:
            prob = altintegral(mlens[i],
                               merr[i],
                               mtrue[i],
                               logmu,
                               sigma)


        sumlogprob += log(prob)

    return sumlogprob

        
