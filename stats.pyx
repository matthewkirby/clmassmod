########################
# Common statistics functions that need to run quickly
#######################
# Compiling info: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /u/ki/dapple/include/python2.7/ -I /u/ki/dapple/lib/python2.7/site-packages/numpy/core/include/ -o stats.so stats.c
#AIfA: gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I /users/dapple/anaconda/include/python2.7/ -I /users/dapple/anaconda/lib/python2.7/site-packages/numpy/core/include/ -o stats.so stats.c


########################

# cython: profile=True

import numpy as np
cimport numpy as np
cimport cython


cdef extern from "math.h":
    double exp(double)
    double log(double)
    double sqrt(double)
    double pow(double, double)

cdef double sqrt2pi = sqrt(2*np.pi)
cdef double twopi = 2*np.pi

#########################

@cython.boundscheck(False)
@cython.wraparound(False)
def Gaussian(np.ndarray[np.double_t, ndim=1, mode='c'] x,
             double mu,
             double sig):

    cdef Py_ssize_t i, nmax

    nmax = x.shape[0]
    
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] result = np.zeros(nmax, dtype = np.float64)

    for i from nmax > i >= 0:
        result[i] = exp(-0.5*(x[i]-mu)**2/sig**2)/(sqrt2pi*sig)

    return result

#########################



@cython.boundscheck(False)
@cython.wraparound(False)
def LogSumGaussian(np.ndarray[np.double_t, ndim=1, mode='c'] x,
             double mu,
             double sig):

    cdef Py_ssize_t i, nmax

    nmax = x.shape[0]
    
    cdef double sum = 0.

    for i from nmax > i >= 0:
        sum += exp(-0.5*(x[i]-mu)**2/sig**2)/(sqrt2pi*sig)

    return log(sum)


#####################

@cython.boundscheck(False)
@cython.wraparound(False)
def LogSum2DGaussian(np.ndarray[np.double_t, ndim=2, mode='c'] xs,
                     np.ndarray[np.double_t, ndim=1, mode='c'] weights,
                     np.ndarray[np.double_t, ndim=1, mode='c'] mu,
                     np.ndarray[np.double_t, ndim=2, mode='c'] invcovar,
                     double sqrtdetcovar):

    cdef Py_ssize_t i, nmax, ndim

    nmax = xs.shape[0]

    cdef double pipow = twopi
    cdef double delta1, delta2
    cdef double sum = 0.
    cdef double chisq = 0.

    for i from nmax > i >= 0:
        delta1 = xs[i,0] - mu[0]
        delta2 = xs[i,1] - mu[1]
        chisq = invcovar[0,0]*delta1*delta1 + invcovar[1,1]*delta2*delta2 + (invcovar[0,1] + invcovar[1,0])*delta1*delta2
        sum += weights[i]*exp(-0.5*chisq)/(pipow*sqrtdetcovar)

    return log(sum)

    

#####################

@cython.boundscheck(False)
@cython.wraparound(False)
def LogSumLogNormal(np.ndarray[np.double_t, ndim=1, mode='c'] x,
                    np.ndarray[np.double_t, ndim=1, mode='c'] logx,
                    double mu,
                    double sig):

    cdef Py_ssize_t i, nmax

    nmax = x.shape[0]
    
    cdef double sum = 0.

    for i from nmax > i >= 0:
        sum += exp(-0.5*(logx[i]-mu)**2/sig**2)/(sqrt2pi*sig*x[i])


    return log(sum)



########################


@cython.boundscheck(False)
@cython.wraparound(False)
def kelly_like(np.ndarray[np.double_t, ndim=1, mode='c'] x,
              np.ndarray[np.double_t, ndim=1, mode='c'] xerr2,
              np.ndarray[np.double_t, ndim=1, mode='c'] y,
              np.ndarray[np.double_t, ndim=1, mode='c'] yerr2,
              np.ndarray[np.double_t, ndim=1, mode='c'] xycovar,
              double alpha,
              double beta,
              double sigint2,
              np.ndarray[np.double_t, ndim=1, mode='c'] pis,
              np.ndarray[np.double_t, ndim=1, mode='c'] mus,
              np.ndarray[np.double_t, ndim=1, mode='c'] tau2):





    cdef Py_ssize_t ngauss = pis.shape[0]
    cdef Py_ssize_t ndat = x.shape[0]

    cdef double logp = 0.


    cdef double V00, V01, V11, delta0, delta1, detV, chisq, curp
    cdef Py_ssize_t curgauss, curdat

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] predictions = np.zeros(ngauss)
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] betatau2 = np.zeros(ngauss)    
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] beta2tau2pSigint2 = np.zeros(ngauss)
    for curgauss from ngauss > curgauss >= 0:
        predictions[curgauss] = alpha + beta*mus[curgauss]
        betatau2[curgauss] = beta*tau2[curgauss]
        beta2tau2pSigint2[curgauss] = beta*beta*tau2[curgauss] + sigint2


    for curdat from ndat > curdat >= 0:

        curp = 0.

        for curgauss from ngauss > curgauss >= 0:

            delta0 = y[curdat] - predictions[curgauss]
            delta1 = x[curdat] - mus[curgauss]

            V00 = beta2tau2pSigint2[curgauss] + yerr2[curdat]
            V01 = betatau2[curgauss] + xycovar[curdat]
            V11 = tau2[curgauss] + xerr2[curdat]

            detV = (V00*V11) - (V01*V01)

            chisq = (V11*delta0*delta0 - 2*V01*delta0*delta1 + V00*delta1*delta1)/detV

            curp = curp + pis[curgauss]*exp(-0.5*chisq)/(twopi*sqrt(detV))

        logp = logp + log(curp)


    return logp

    
