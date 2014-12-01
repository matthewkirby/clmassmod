import stats
import multiprocessing as mp
import numpy as np

def LogSum2DGaussianWrapper(index):

    xs = datashare.getLikeSamples(index)
    mu = datashare.getMean(index)
    invcovar = datashare.getInvCovar(index)
    sqrtdetcovar = datashare.getSqrtDetCovar(index)


    return stats.LogSum2DGaussian(xs = xs,
                                  mu = mu,
                                  invcovar = invcovar,
                                  sqrtdetcovar = sqrtdetcovar)


###################

class DataShare(object):

    def __init__(self):
        self.nclusters = None
        self.nsamples = None
        self.like_samples = None
        self.means = None
        self.invcovars = None
        self.sqrtdetcovars = None

    def setup(self, nclusters, nsamples):
        self.nclusters = nclusters
        self.nsamples = nsamples
        self.binassignment = mp.Array('i', nclusters, lock=False)
        self.like_samples = mp.Array('d', 2*nclusters*nsamples, lock=False)
        self.means = mp.Array('d', 2*nclusters, lock=False)
        self.invcovars = mp.Array('d', nclusters*4, lock=False)
        self.sqrtdetcovars = mp.Array('d', nclusters, lock=False)

        

    def getLikeSamples(self, index):
        return np.array(self.like_samples[2*self.nsamples*index:2*self.nsamples*(index+1)]).reshape((self.nsamples,2))

    def setLikeSamples(self, index, like_samples):
        self.like_samples[2*self.nsamples*index:2*self.nsamples*(index+1)] = like_samples.flatten()

    def getMean(self, index):
        return np.array(self.means[2*index:2*(index+1)])

    def setMean(self, index, means):
        self.means[2*index:2*(index+1)] = means

    def getInvCovar(self, index):
        return np.array(self.invcovars[4*index:4*(index+1)]).reshape((2,2))

    def setInvCovar(self, index, invcovar):
        self.invcovars[4*index:4*(index+1)] = invcovar.flatten()

    def getSqrtDetCovar(self, index):
        return self.sqrtdetcovars[index]

    def setSqrtDetCovar(self, index, sqrtdetcovar):
        self.sqrtdetcovars[index] = sqrtdetcovar

###################

datashare = DataShare()
