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

#########

@cython.boundscheck(False)
@cython.wraparound(False)
def pdfintegral(np.ndarray[np.double_t, ndim=1, mode='c'] ml_ints,
                np.ndarray[np.double_t, ndim=1, mode='c'] deltamasses,
                np.ndarray[np.double_t, ndim=1, mode='c'] delta_logmls,
                np.ndarray[np.double_t, ndim=1, mode='c'] pdf,
                double logmu, 
                double sigma):
    ### assumes that ml_ints = 0 is not included
    ### but that there is an implicit ml_ints = 0 calculation.
    ### deltamasses should be one less in length than ml_ints

    cdef Py_ssize_t i, nsamples
    nsamples = ml_ints.shape[0]

    cdef double thesum
    thesum = 0.

    cdef double sigma2, sigmasqrt2pi
    neg2sigma2 = -2*(sigma**2)
    sigmasqrt2pi = sigma*sqrt2pi

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] lognormpart = np.zeros(nsamples)

    lognormpart[0] = exp((delta_logmls[0]-logmu)**2/neg2sigma2)/(sigmasqrt2pi*ml_ints[0])
    thesum += 0.5*ml_ints[0]*lognormpart[0]*pdf[0]  #ml_ints[0] is deltamasses[-1]; integrand @ ml=0 is 0

    for i from 1 <= i < nsamples:

        lognormpart[i] = exp((delta_logmls[i]-logmu)**2/neg2sigma2)/(sigmasqrt2pi*ml_ints[i])

        thesum += 0.5*deltamasses[i-1]*(lognormpart[i]*pdf[i] + lognormpart[i-1]*pdf[i-1])


    return thesum



#########


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



@cython.boundscheck(False)
@cython.wraparound(False)
def mcmcloglinearlike(np.ndarray[np.double_t, ndim=2, mode='c'] ml_ints,
                      np.ndarray[np.double_t, ndim=2, mode='c'] delta_logmls,
                      np.ndarray[np.int_t, ndim=1, mode='c'] ngoodsamples,
                      double logmu, 
                      double sigma):


    cdef Py_ssize_t i, nclusters, nsamples
    nclusters = ml_ints.shape[0]

    cdef double sumlogprob = 0.
    cdef double prob = 0.

    

    for i from nclusters > i >= 0:

        nsamples = ngoodsamples[i]
        
        prob = altintegral(ml_ints[i,:nsamples],
                           delta_logmls[i,:nsamples],
                           logmu,
                           sigma)



        sumlogprob += log(prob)

    return sumlogprob


##############


@cython.boundscheck(False)
@cython.wraparound(False)
def pdfloglinearlike(np.ndarray[np.double_t, ndim=1, mode='c'] ml_ints,
                      np.ndarray[np.double_t, ndim=1, mode='c'] deltamasses,
                      np.ndarray[np.double_t, ndim=2, mode='c'] delta_logmls,
                      np.ndarray[np.double_t, ndim=2, mode='c'] pdfs,
                      double logmu, 
                      double sigma):


    cdef Py_ssize_t i, nclusters
    nclusters = pdfs.shape[0]


    cdef double sumlogprob = 0.
    cdef double prob = 0.

    

    for i from nclusters > i >= 0:
        
        prob = pdfintegral(ml_ints,
                           deltamasses,
                           delta_logmls[i,:],
                           pdfs[i,:],
                           logmu,
                           sigma)



        sumlogprob += log(prob)

    return sumlogprob

    


@cython.boundscheck(False)
@cython.wraparound(False)
def outlierloglinearlike(np.ndarray[np.double_t, ndim=2, mode='c'] ml_ints,
                         np.ndarray[np.double_t, ndim=2, mode='c'] delta_logmls,
                         np.ndarray[np.double_t, ndim=2, mode='c'] outlier_ml_ints,
                         np.ndarray[np.double_t, ndim=2, mode='c'] outlier_delta_logmls,
                         double logmu, 
                         double sigma,
                         double fracoutliers):


    cdef Py_ssize_t i, nclusters
    nclusters = ml_ints.shape[0]

    cdef double sumlogprob = 0.
    cdef double prob = 0.
    cdef outlierprob = 0.
    

    for i from nclusters > i >= 0:


        
        prob = altintegral(ml_ints[i,:],
                           delta_logmls[i,:],
                           logmu,
                           sigma)

        outlierprob = altintegral(outlier_ml_ints[i,:],
                                  delta_logmls[i,:],
                                  logmu,
                                  sigma)



        sumlogprob += log((1-fracoutliers)*prob + fracoutliers*outlierprob)

    return sumlogprob

        
###########


@cython.boundscheck(False)
@cython.wraparound(False)
def pdfGaussMix1D(np.ndarray[np.double_t, ndim=2, mode='c'] delta_mls,
                  np.ndarray[np.double_t, ndim=2, mode='c'] delta_masses,
                  np.ndarray[np.double_t, ndim=2, mode='c'] pdfs,
                  np.ndarray[np.double_t, ndim=1, mode='c'] pis,
                  np.ndarray[np.double_t, ndim=1, mode='c'] mus,
                  np.ndarray[np.double_t, ndim=1, mode='c'] tau2):
                  


    cdef Py_ssize_t i, j, nclusters, ngauss, nmasses
    nclusters = pdfs.shape[0]
    ngauss = pis.shape[0]
    nmasses = delta_mls.shape[1]

    ### build intrinsic pdf

    #normalization
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] norms = np.zeros(ngauss)
    for j from ngauss > j >= 0:
        norms[j] = 1./(sqrt2pi*sqrt(tau2[j]))
    
#    #exponentials
#    cdef np.ndarray[np.double_t, ndim=1, mode='c'] intrinsicpdf = np.zeros(nmasses)

#    for i from nmasses > i >= 0:
#        for j from ngauss > j >= 0:
#            deltamass_mu2 = (deltamasses[i] - mus[j])**2
#            gausseval = pis[j]*exp(-0.5*deltamass_mu2/tau2[j])*norms[j]
#            intrinsicpdf[i] = intrinsicpdf[i] + gausseval
#            
            
    ### convolve with pdfs

    cdef double sumlogprob = 0.
    cdef double deltamass_mu2, gausseval   
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] integrand = np.zeros(nmasses)

    for i from nclusters > i >= 0:

        for j from nmasses > j >= 0:
            integrand[j] = 0.
        
        for j from nmasses > j >= 0:

            for k from ngauss > k >= 0:

                deltamass_mu2 = (delta_mls[i,j] - mus[k])**2
                gausseval = pis[k]*exp(-0.5*deltamass_mu2/tau2[k])*norms[k]
                integrand[j] += gausseval

            integrand[j] *= pdfs[i,j]
            
        #trapezoid rule
        prob = 0.
        for j from nmasses-1 > j >= 0:
            prob += delta_masses[i,j]*(integrand[j+1] + integrand[j])
        prob *= 0.5

        sumlogprob += log(prob)

    return sumlogprob

