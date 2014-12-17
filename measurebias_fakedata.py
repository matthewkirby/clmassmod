##########################
# Create fake data that tests measurebias.py
#########################


import numpy as np
import scipy.stats
import cPickle
import nfwmodeltools as tools, nfwutils

nclusters = 4000

true_logmasses = np.log(10**np.random.uniform(14., 15.5, nclusters))
true_logconcens = np.log(4.*np.ones(nclusters))

true_massbias = 0.8
true_massscatter = 0.02
true_concen_scatter = 0.02
true_mc_cov = 0.65

means = np.column_stack([np.log(true_massbias) + true_logmasses, true_logconcens])
covar = np.array([[true_massscatter**2, true_mc_cov*true_massscatter*true_concen_scatter],
                  [true_mc_cov*true_massscatter*true_concen_scatter, true_concen_scatter**2]])

deltas = scipy.stats.multivariate_normal.rvs(cov = covar, size=nclusters)

ml_estimate = means + deltas

answers = {}

for i in range(nclusters):

    clusterinfo = dict(m200 = np.exp(true_logmasses[i]),
                                  c200 = np.exp(true_logconcens[i]),
                                  lm200 = np.exp(ml_estimate[i,0]),
                                  lc200 = np.exp(ml_estimate[i,1]))

    answers['halo_%d' % i] = clusterinfo

    #need to make shear profiles that I then fit.

    ###TODO AFTER LUNCH: WRITE FILEREADER FOR NFWFIT THAT GENERATES FAKE DATA FOR THIS, AND NEW CONFIG FILES TO RUN WITH NFWFIT

    with open('../measurebias_fakedata/highsn/halo_%d.info' % i, 'wb') as output:
        cPickle.dump(clusterinfo, output)
    


with open('measurebias_fakedata__highsn_answers.pkl', 'wb') as output:
    cPickle.dump(answers, output)



