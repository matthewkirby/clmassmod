import stats
import numpy as np

def LogSum2DGaussianWrapper(args):

    bin_number = args['number']
    bin_mass_scatter = args['mass_scatter']
    bin_mc_covar = args['mc_covar']
    bin_c200_scatter = args['c200_scatter']
    
    bin_covar = np.array([[bin_mass_scatter**2, bin_mc_covar*bin_mass_scatter*bin_c200_scatter],
                                [bin_mc_covar*bin_mass_scatter*bin_c200_scatter, bin_c200_scatter**2]])

    bin_invcovar = np.linalg.inv(bin_covar)
    bin_invsqrtdetcovar = 1./np.sqrt(np.linalg.det(bin_covar))

    bin_ratio = args['logmassratio']
    bin_logc200 = np.log(args['c200'])


    clustersinbin = datastore.getClustersInBin(bin_number)
    nInBin = len(clustersinbin)

    cluster_logprobs = np.zeros(nInBin)

    for i, cluster in enumerate(clustersinbin):

        mu0 = bin_ratio + datastore.log_mtrues[cluster]
        mu1 = bin_logc200

        cluster_logprobs[i] =  stats.LogSum2DGaussian(samples0 = datastore.getLogMassSamples(cluster),
                                                      samples1 = datastore.getLogC200Samples(cluster),
                                                      weights = datastore.getWeights(cluster),
                                                      mu0 = mu0,
                                                      mu1 = mu1,
                                                      invcovar = bin_invcovar,
                                                      invsqrtdetcovar = bin_invsqrtdetcovar)

    return np.sum(cluster_logprobs)

###################



class DataStore(object):

    def __init__(self):
        self.log_mtrues = None
        self.logmass_samples = None
        self.logc200_samples = None
        self.weights = None
        self.bin_assignments = None
        self.clusters_inbin = None

    def getLogMTrues(self, cluster):
        return self.log_mtrues[cluster]

    def getLogMassSamples(self, cluster):
        return self.logmass_samples[cluster]

    def getLogC200Samples(self, cluster):
        return self.logc200_samples[cluster]

    def getWeights(self, cluster):
        return self.weights[cluster]

    def getClustersInBin(self, bin):
        return self.clusters_inbin[bin]


    


datastore = DataStore()
