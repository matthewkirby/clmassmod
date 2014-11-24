##########################
# Create fake data that tests measurebias.py
#########################


import numpy as np
import scipy.stats

nclusters = 2000

true_logmasses = np.log(10**np.random.uniform(14., 16., nclusters))
true_logconcens = np.log(4.*np.ones(nclusters))

true_massbias = 0.8
true_massscatter = 0.1
true_concen_scatter = 0.5
true_mc_cov = 0.65

means = np.column_stack([np.log(true_massbias) + true_logmasses, true_logconcens])
covar = np.array([[true_massscatter**2, true_mc_cov*true_massscatter*true_concen_scatter],
                  [true_mc_cov*true_massscatter*true_concen_scatter, true_concen_scatter**2]])

deltas = scipy.stats.rvs(cov = covar, size=nclusters)

ml_estimate = means + deltas

for i in range(nclusters):

    mass_samples = ml_estimate[i,0] + 0.2*np.random.standard_normal(size=2000)
    concen_samples = ml_estimate[i,1] + 0.2*np.random.standard_normal(size=2000)

    data = dict(logM200 = mass_samples,
                c200 = np.exp(concen_samples),
                m200 = np.exp(mass_samples))

    with open('measurebias_fakedata/halo_%d.out' % i, 'wb') as output:
        cPickle.dump(data, output)



