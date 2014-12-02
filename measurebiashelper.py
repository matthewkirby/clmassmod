import stats
import numpy as np

def LogSum2DGaussianWrapper(cluster):

    return stats.LogSum2DGaussian(xs = datastore.getLikeSamples(cluster),
                                  weights = datastore.getWeights(cluster),
                                      mu = datastore.getMean(cluster),
                                      invcovar = datastore.getInvCovar(cluster),
                                      sqrtdetcovar = datastore.getSqrtDetCovar(cluster))

###################



class DataStore(object):

    def __init__(self):
        self.log_mtrues = None
        self.like_samples = None
        self.weights = None
        self.bin_assignments = None
        self.bin_ratios = None
        self.bin_logc200s = None
        self.invcovars = None
        self.sqrtdetcovars = None

    def getWeights(self, cluster):
        return self.weights[cluster]

    def getLikeSamples(self, cluster):
        return self.like_samples[cluster]

    def getMean(self, cluster):

        curBin = self.bin_assignments[cluster]

        curRatio = self.bin_ratios[curBin]
        logC200 = self.bin_logc200s[curBin]

        return np.array([curRatio + self.log_mtrues[cluster],
                         logC200])

    def getInvCovar(self, cluster):

        return self.invcovars[self.bin_assignments[cluster]]

    def getSqrtDetCovar(self, cluster):
        
        return self.sqrtdetcovars[self.bin_assignments[cluster]]

    


datastore = DataStore()
