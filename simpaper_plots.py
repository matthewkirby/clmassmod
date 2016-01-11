##############
# Plots for the sims paper
#############

import pylab
import numpy as np

import nfwnoise
import nfwutils
import varcontainer
import publication_plots as pp
import basicMassCon


############


def likelihoodPlots():
    '''Show the lensing likelihoods, to make the point that yes, they are almost symetric in linear space, but Gaussian approximations are poor.'''

    linfig = pylab.figure()
    linax = linfig.add_subplot(1,1,1)

    logfig = pylab.figure()
    logax = logfig.add_subplot(1,1,1)

    figs = [linfig, logfig]

    true_masses = [2e14, 1e15]
    linestyles = ['-', '--']

    for true_mass, linestyle in zip(true_masses, linestyles):

        m200s = np.arange(-5.05e14, 3e15, 1e13)
        r_mpc_edges = np.linspace(0.75, 3., 13)
        r_mpc = (r_mpc_edges[1:] + r_mpc_edges[:-1])/2.
        zcluster = 0.415 #median of megacam

        duffy = basicMassCon.Duffy(None)
        c_duff = duffy(true_mass * nfwutils.global_cosmology.h, zcluster, 200.)



        beta_s = 0.49 #median of megacam

        g_truth = nfwnoise.createPerfectProfile(true_mass, c_duff, zcluster, r_mpc, beta_s)

        ngalsperarcmin = 9.7 #median of megacam set
        r_arcmin_edges = (r_mpc_edges/nfwutils.global_cosmology.angulardist(zcluster))*(180./np.pi)*60
        bin_areas = np.pi*(r_arcmin_edges[1:]**2 - r_arcmin_edges[:-1]**2)
        ngals_per_bin = ngalsperarcmin * bin_areas
        g_err = 0.25 / np.sqrt(ngals_per_bin)





        #calcualte pdf for duffy

        config = varcontainer.VarContainer()
        config.massconmodule = 'basicMassCon'

        config.massconrelation = 'Duffy'
        likelihood = nfwnoise.Likelihood(config)
        likelihood.bindData(r_mpc, g_truth, g_err, beta_s, zcluster)



        duffy_logpdf = np.array([likelihood(m200) for m200 in m200s])
        duffy_pdf = np.exp(duffy_logpdf - np.max(duffy_logpdf))
        duffy_pdf = duffy_pdf / np.trapz(duffy_pdf, m200s)


        #plot


        linax.plot(m200s, duffy_pdf, 
                   marker='None', linestyle=linestyle, 
                   color = 'k',
                   linewidth = 2., 
                   label='$M_{200} = $%1.1f x $10^{15} M_{\odot}$' % (true_mass / 1e14))

        linax.axvline(true_mass, c=(0.6, 0.6, 0.6), linewidth=2, linestyle=linestyle)


        posmass = m200s > 0
        log10m200s = np.log10(m200s[posmass])
        log10pdf = duffy_pdf[posmass]
        log10pdf = log10pdf / np.trapz(log10pdf, log10m200s)
        logax.plot(log10m200s, log10pdf, 
                   marker='None', linestyle=linestyle, 
                   color = 'k',
                   linewidth = 2., 
                   label='$M_{200} = $%1.1f x $10^{15} M_{\odot}$' % (true_mass / 1e14))

        logax.axvline(np.log10(true_mass), c=(0.6, 0.6, 0.6), linewidth=2, linestyle=linestyle)



    linax.legend()
    linax.set_xlabel(r'$M_{200} [M_{\odot}]$', fontsize=18)
    linax.set_ylabel(r'Likelihood($M_{200}$)', fontsize=18)
    linax.set_xlim(-0.5e15, 2.5e15)

    linfig.tight_layout()

    linfig.savefig('figures/likelihood_shape.png')
    linfig.savefig('figures/likelihood_shape.eps')
    linfig.savefig('figures/likelihood_shape.pdf')


    logax.legend(loc='upper left')
    logax.set_xlabel(r'$\mathrm{log}_{10}M_{200}$', fontsize=18)
    logax.set_ylabel(r'Likelihood($\mathrm{log}_{10}M_{200}$)', fontsize=18)
    logax.set_xlim(np.log10(np.min(m200s[posmass])), np.log10(2.5e15))
    logfig.tight_layout()
    logfig.savefig('figures/likelihood_logshape.png')
    logfig.savefig('figures/likelihood_logshape.eps')
    logfig.savefig('figures/likelihood_logshape.pdf')

    return figs

##############
    
