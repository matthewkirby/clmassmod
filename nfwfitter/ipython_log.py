# IPython log file

get_ipython().magic(u'logstart')
get_ipython().magic(u'pylab')
import concentrationfittools as cfittools
import concentrationfit as cfit
import pymc
mock_data = cfit.generate_mock_data(0.95, 0.15, 4.5, 0.06, .02, .02)
model = cfit.buildMCMCModel_mc(mock_data)
pymc_maximum_aposteriori = pymc.MAP(model)
pymc_maximum_aposteriori.fit()
np.exp(pymc_maximum_aposteriori.logmu_c.value)
pymc_model = pymc.Model(model)
pymc_MCMC = pymc.MCMC(model)
pymc_MCMC.sample(1000)
pyplot.hist2d(pymc_MCMC.trace('logmu_c')[:]),pymc_MCMC.trace('logmu')[:], bins=50)
pyplot.hist2d(pymc_MCMC.trace('logmu_c')[:],pymc_MCMC.trace('logmu')[:], bins=50)
pyplot.hist2d(np.exp(pymc_MCMC.trace('logmu_c')[:]),np.exp(pymc_MCMC.trace('logmu')[:]), bins=50)
pyplot.hist2d(np.exp(pymc_MCMC.trace('logmu_c')[:]),np.exp(pymc_MCMC.trace('logmu')[:]), bins=50)
pyplot.hist2d(np.exp(pymc_MCMC.trace('logmu_c')[500:]),np.exp(pymc_MCMC.trace('logmu')[500:]), bins=50)
plt.hist(np.exp(pymc_MCMC.trace('logmu_c')[:]),bins=25)
exit()
