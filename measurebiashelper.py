import stats
import numpy as np

def LogSum2DGaussianWrapper(cluster):

    mu0,mu1 = datastore.getMean(cluster)

    return stats.LogSum2DGaussian(samples0 = datastore.getLogMassSamples(cluster),
                                  samples1 = datastore.getLogC200Samples(cluster),
                                  weights = datastore.getWeights(cluster),
                                  mu0 = mu0,
                                  mu1 = mu1,
                                  invcovar = datastore.getInvCovar(cluster),
                                  invsqrtdetcovar = datastore.getInvSqrtDetCovar(cluster))

###################



class DataStore(object):

    def __init__(self):
        self.log_mtrues = None
        self.logmass_samples = None
        self.logc200_samples = None
        self.weights = None
        self.bin_assignments = None
        self.bin_ratios = None
        self.bin_logc200s = None
        self.invcovars = None
        self.invsqrtdetcovars = None

    def getLogMassSamples(self, cluster):
        return self.logmass_samples[cluster]

    def getLogC200Samples(self, cluster):
        return self.logc200_samples[cluster]

    def getWeights(self, cluster):
        return self.weights[cluster]

    def getMean(self, cluster):

        curBin = self.bin_assignments[cluster]

        curRatio = self.bin_ratios[curBin]
        logC200 = self.bin_logc200s[curBin]

        return curRatio + self.log_mtrues[cluster], logC200

    def getInvCovar(self, cluster):

        return self.invcovars[self.bin_assignments[cluster]]

    def getInvSqrtDetCovar(self, cluster):
        
        return self.invsqrtdetcovars[self.bin_assignments[cluster]]

    


datastore = DataStore()
