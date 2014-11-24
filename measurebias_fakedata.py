##########################
# Create fake data that tests measurebias.py
#########################


import numpy as np
import scipy.stats
import cPickle

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

deltas = scipy.stats.multivariate_normal.rvs(cov = covar, size=nclusters)

ml_estimate = means + deltas

answers = {}

for i in range(nclusters):

    mass_samples = ml_estimate[i,0] + 0.2*np.random.standard_normal(size=250)
    concen_samples = ml_estimate[i,1] + 0.2*np.random.standard_normal(size=250)

    data = dict(logM200 = mass_samples,
                c200 = np.exp(concen_samples),
                m200 = np.exp(mass_samples))

    answers['halo_%d' % i] = dict(m200 = np.exp(true_logmasses[i]))

    with open('../measurebias_fakedata/halo_%d.out' % i, 'wb') as output:
        cPickle.dump(data, output)


with open('measurebias_fakedata_answers.pkl', 'wb') as output:
    cPickle.dump(answers, output)



