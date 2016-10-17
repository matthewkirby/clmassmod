########
# Wrapper Around Voigt Function defined in C
#########
# Compiling info: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /u/ki/dapple/include/python2.7/ -I /u/ki/dapple/lib/python2.7/site-packages/numpy/core/include/ -I./ -o voigt_tools.so voigt_tools.c voigt.c

#

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "math.h":
    double log(double)
    double sqrt(double)
    bint isnan(double)

cdef extern from "voigt.h":
    double voigt(double, double, double)

######################

__cvs_id__ = "$Id: voigt_tools.pyx,v 1.5 2011-01-19 22:39:44 dapple Exp $"

######################

DTYPE = np.double
ctypedef np.double_t DTYPE_T

#######################

def voigtProfile(np.ndarray[DTYPE_T, ndim=1, mode='c'] x,
          double sigma,
          double gamma):

    cdef Py_ssize_t nobjs = x.shape[0]
    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] results = np.zeros(nobjs)

    cdef int i = nobjs
    for i from nobjs > i >= 0:
        results[i] = voigt(x[i], sigma, gamma)

    return results

###############

def voigtSamples(double sigma, double gamma, int size, limits=(-5,5), double binsize = 0.0001):

    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] xbins = np.arange(limits[0], limits[1], binsize)

    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] probs = voigtProfile(xbins, sigma, gamma)
    cdef int nbins = len(xbins)

    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] cdf = np.cumsum(probs)
    cdf = cdf / cdf[nbins - 1]


    cdef np.ndarray[DTYPE_T, ndim=1, mode='c'] picks = np.random.uniform(size = size)

    positions = np.searchsorted(cdf, picks, side='left')

    return xbins[positions]    


###############


    

###############
    
@cython.boundscheck(False)
@cython.wraparound(False)
def likelihood(np.ndarray[DTYPE_T, ndim=1, mode='c'] g not None, 
               np.ndarray[DTYPE_T, ndim=1, mode='c'] mu not None, 
               double sigma, 
               double gamma):

    cdef Py_ssize_t nobjs = g.shape[0]
    cdef double x = 0.
    cdef double prob = 0.

    cdef double tot_logprob = 0.
    cdef int i = nobjs
    for i from nobjs > i >=0:
        
        x = g[i] - mu[i]
        prob = voigt(x, sigma, gamma)

        if isnan(prob) or prob <= 0:
            return -np.infty

        tot_logprob = tot_logprob + log(prob)


    return tot_logprob




#####################

cdef double sqrt2pi = np.sqrt(2*np.pi)

def gauss(double x, double sigma):

    return np.exp(-0.5*(x/sigma)**2)/(sigma*sqrt2pi)


####################

def doublegauss(double x,
                double sigma1,
                double sigma2scale,
                double alpha):


    cdef double curx
    cdef double sigma2 = sigma1*sigma2scale
    cdef double comp_alpha = 1-alpha
    return comp_alpha*gauss(x, sigma1) + alpha*gauss(x, sigma2)






@cython.boundscheck(False)
@cython.wraparound(False)
def doublegauss_likelihood(np.ndarray[DTYPE_T, ndim=1, mode='c'] g not None, 
                           np.ndarray[DTYPE_T, ndim=1, mode='c'] mu not None, 
                           double sigma, 
                           double sigma2scale,
                           double alpha):

    cdef Py_ssize_t nobjs = g.shape[0]
    cdef double x = 0.
    cdef double prob = 0.

    cdef double tot_logprob = 0.
    cdef int i = nobjs
    for i from nobjs > i >=0:
        
        x = g[i] - mu[i]
        prob = doublegauss(x, sigma, sigma2scale, alpha)

        if isnan(prob) or prob <= 0:
            return -np.infty

        tot_logprob = tot_logprob + log(prob)


    return tot_logprob


