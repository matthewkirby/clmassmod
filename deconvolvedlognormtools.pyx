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
def integral(np.ndarray[np.double_t, ndim=1, mode='c'] ml_int, 
              double ml,
              double merr, 
              double mtrue,
              double logmu, 
              double sigma):

    cdef Py_ssize_t i, nmax
    nmax = ml_int.shape[0]
    cdef double normpart, lognormpart

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] result = np.zeros(nmax, dtype=np.float64)

    for i from nmax > i >= 0:
        normpart = exp(-0.5*(ml_int[i]-ml)**2/merr**2)/(sqrt2pi*merr)
        lognormpart = exp(-0.5*(log(ml_int[i])-logmu)**2/sigma**2)/(sqrt2pi*sigma*ml_int[i])
        result[i] = normpart*lognormpart

    #trapezoid scheme
    cdef double sum = 0.
    for i from nmax-1 > i >= 1:
        sum += 2*result[i]
    sum += result[0] + result[nmax-1]


    return sum*(ml_int[nmax-1] - ml_int[0])/(2*nmax)

    


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

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] ml_int_grid
    
    cdef double low, high

    low = exp(logmu - 5*sigma)
    high = exp(logmu + 5*sigma)


    for i from nclusters > i >= 0:

        ml_int_grid = np.linspace(low*mtrue[i], high*mtrue[i], 200)
        
        prob = integral(ml_int_grid,
                        mlens[i],
                        merr[i],
                        mtrue[i],
                        logmu,
                        sigma)

        sumlogprob += log(prob)

    return sumlogprob

        
