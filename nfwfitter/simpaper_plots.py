##############
# Plots for the sims paper
#############

import pylab
import numpy as np
import cPickle

import nfwnoise
import nfwutils
import nfwfit
import varcontainer
import publication_plots as pp
import basicMassCon
import deconvolvedlognorm
import pymc
import fitmodel
import confidenceinterval as ci


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
                   label='$M_{200} = $%1.1f x $10^{14} M_{\odot}$' % (true_mass / 1e14))

        linax.axvline(true_mass, c=(0.6, 0.6, 0.6), linewidth=2, linestyle=linestyle)


        posmass = m200s > 0
        log10m200s = np.log10(m200s[posmass])
        log10pdf = duffy_pdf[posmass]
        log10pdf = log10pdf / np.trapz(log10pdf, log10m200s)
        logax.plot(log10m200s, log10pdf, 
                   marker='None', linestyle=linestyle, 
                   color = 'k',
                   linewidth = 2., 
                   label='$M_{200} = $%1.1f x $10^{14} M_{\odot}$' % (true_mass / 1e14))

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
    

def calcMaxlikeDistro(outputfile = 'maxlikedistro_paperplot_lowsn.pkl', ngalsperarcmin=5):
    '''show that using maximum likelihood estimators does not recover truth'''

    

    figs = []

    #megacamstyle

    ntrials = 2000
    samplesize = 25

    r_mpc_edges = np.linspace(0.75, 3., 13)
    r_mpc = (r_mpc_edges[1:] + r_mpc_edges[:-1])/2.
    zcluster = 0.3 
    duffy = basicMassCon.Duffy(None)
    c200 = 4.0
    beta_s = nfwutils.global_cosmology.beta_s([1.0], zcluster)[0]

    r_arcmin_edges = (r_mpc_edges/nfwutils.global_cosmology.angulardist(zcluster))*(180./np.pi)*60
    bin_areas = np.pi*(r_arcmin_edges[1:]**2 - r_arcmin_edges[:-1]**2)
    ngals_per_bin = ngalsperarcmin * bin_areas
    g_err = 0.25 / np.sqrt(ngals_per_bin)

    m200s = np.arange(1e12, 5e15, 1e13)
    
    config = varcontainer.VarContainer()
    config.massconmodule = 'basicMassCon'
    config.massconrelation = 'constant'
    config.concentration = c200
    model = nfwfit.buildModel(config)
    likelihood = nfwnoise.Likelihood(config)


    unweightedbootstraps_mean = np.zeros(ntrials)
    weightedaverages = np.zeros(ntrials)
    wlssaverage = np.zeros(ntrials)
    linear_fit = np.zeros(ntrials)
    smithaverages = np.zeros(ntrials)
    pdfaverages = np.zeros(ntrials)

    for curtrial in range(ntrials):
        
        true_masses = 10**np.random.uniform(np.log10(5e14), np.log10(2e15), samplesize)
        lens_masses = true_masses * np.exp(0.2*np.random.standard_normal(size=samplesize))

        #compute individual cluster likelihoods
        maxlike_masses = np.zeros(samplesize)
        maxlike_errs = np.zeros((2,samplesize))
        pdfs = np.zeros((samplesize, len(m200s)))
        logpdf = np.zeros(len(m200s))
    
        for curcluster in range(samplesize):
        
            g_obs = nfwnoise.createPerfectProfile(lens_masses[curcluster], c200, 
                                                  zcluster, r_mpc, beta_s) + \
                g_err*np.random.standard_normal(size=len(g_err))

            fitter = nfwfit.NFWFitter(None, model, config)

            maxlikemass, maxlikeerr = fitter.minChisqMethod(r_mpc, g_obs, g_err,
                                                            beta_s, beta_s**2, zcluster,
                                                            guess = [true_masses[curcluster]/model.massScale])

            maxlike_masses[curcluster] = maxlikemass['m200']
            try:
                maxlike_errs[:,curcluster] = maxlikeerr['m200']
            except TypeError:
                maxlike_errs[:,curcluster] = (maxlikeerr['m200'], maxlikeerr['m200'])

            likelihood.bindData(r_mpc, g_obs, g_err, beta_s, zcluster)

            for i, m200 in enumerate(m200s):
                logpdf[i] = likelihood(m200)
            pdf = np.exp(logpdf - np.max(logpdf))
            pdf /= np.trapz(pdf, m200s)
            pdfs[curcluster,:] = pdf

        maxlike_masses *= model.massScale
        maxlike_errs = np.abs(maxlike_errs*model.massScale)

        maxlike_logratio = np.log(maxlike_masses / true_masses)




        #combine using different formulisms
        bootstraps = np.zeros(1000)

        for curboot in range(1000):
            booted = np.random.randint(0, samplesize, samplesize)
            bootstraps[curboot] = np.mean(maxlike_logratio[booted])


        unweightedbootstraps_mean[curtrial] = np.mean(bootstraps)



        maxlike_simpleerr = np.sum(maxlike_errs, axis=0)/(2.*maxlike_masses)

        weight = 1./maxlike_simpleerr**2
        weightedaverages[curtrial] = np.sum(weight*maxlike_logratio)/np.sum(weight)
        
        wlssaverage_weights = 1./(maxlike_simpleerr**2 + (0.2*maxlike_masses)**2)
        wlssaverage[curtrial] = np.sum(wlssaverage_weights*maxlike_logratio)/np.sum(wlssaverage_weights)
        
            
        smithweight = 1./(np.sum(maxlike_errs, axis=0)/2.)**2
        smithaverages[curtrial] = np.sum(smithweight*maxlike_logratio)/np.sum(smithweight)

        
        
        halos = [dict(id=curcluster,
                      true_mass = true_masses[curcluster],
                      masses = m200s,
                      pdf = pdfs[curcluster,:]) for curcluster in range(samplesize)]

        
        pdfmodel = deconvolvedlognorm.buildPDFModel(halos)
        pdfmap = pymc.MAP(pdfmodel)
        pdfmap.fit()
        pdfaverages[curtrial] = pdfmap.logmu.value

        

        linefitter = fitmodel.FitModel(true_masses/1e14, maxlike_masses/1e14,
                                       np.sum(maxlike_errs/1e14,axis=0)/2.,
                                       fitmodel.RatioModel, guess = [1.])
        linefitter.fit()
        linear_fit[curtrial] = linefitter.par_vals['a0']




        

    results = dict(unweightedbootstrapsmean = np.exp(unweightedbootstraps_mean),
                   weightedaverages         = np.exp(weightedaverages),
                   wlssaverage              = np.exp(wlssaverage),
                   smithaverages            = np.exp(smithaverages),
                   pdfaverages              = np.exp(pdfaverages),
                   linearfit                = linear_fit
               )

    with open(outputfile, 'wb') as output:
        cPickle.dump(results, output, -1)


###################

def bootstrapdistro(distro, foo = np.mean, nboots=1000):

    booted = np.zeros(nboots)
    for i in range(nboots):
        booted[i] = foo(distro[np.random.randint(0, len(distro), len(distro))])

    return np.mean(booted), np.std(booted)

def maxlikeComparePlots(inputfile = 'maxlikedistro_paperplot_highsn_large.pkl', outputbase = 'maxlike_comparison'):

    with open(inputfile, 'rb') as input:
        distros = cPickle.load(input)

    fig = pylab.figure()
    ax = fig.add_subplot(1,1,1)

    labels = [('unweightedbootstrapsmean',  'Bootstrap Geometric Mean', '^', 'r'),
              ('weightedaverages',          'Weighted Geometric Mean', 's', 'r'),
              ('wlssaverage',               'Weighted+LSS Geometric Mean', 'p', 'r'),
              ('smithaverages',             'Absolute Err Geometric Mean', 'd', 'r'),
              ('linearfit',                 'Maxlike Linear Fit', 'v', 'k'),
              ('pdfaverages',               'PDF w/intrinsic scatter', 'o', 'k')]

    for key, label, marker, color in labels:
        
        bias, biaserr = bootstrapdistro(distros[key])
        stdev, stdeverr = bootstrapdistro(distros[key], np.std)

        ax.errorbar(bias, stdev, stdeverr, biaserr, label=label, marker=marker, color=color, linestyle='none', markersize=10)


    linearfit_bias = [bootstrapdistro(distros['linearfit'])[0]]
    linearfit_precision = [bootstrapdistro(distros['linearfit'], np.std)[0]]
    for inputfile in 'maxlikedistro_paperplot_lowsn.pkl maxlikedistro_paperplot_extralowsn.pkl'.split():
        with open(inputfile, 'rb') as input:
            distros = cPickle.load(input)
        linearfit_bias.append(bootstrapdistro(distros['linearfit'])[0])
        linearfit_precision.append(bootstrapdistro(distros['linearfit'], np.std)[0])
        

    ax.plot(np.array(linearfit_bias), linearfit_precision, 'k--')
    ax.axvline(1.0, c='k', linestyle=':')
    ax.set_xlim(0.8, 1.1)
    ax.set_ylim(0.0, 0.15)
    ax.set_xscale('log')
    ax.set_xticks([0.8, 0.9, 1.0, 1.1])
    ax.set_xticklabels('0.8 0.9 1.0 1.1'.split())
    ax.set_xlabel('Bias [Lensing / True]', fontsize=18)
    ax.set_ylabel('Precision', fontsize=18)
    ax.legend(numpoints=1, ncol=1, loc='lower left', frameon=False)


        
    

    fig.tight_layout()

    fig.savefig('figures/{}.png'.format(outputbase))
    fig.savefig('figures/{}.eps'.format(outputbase))
    fig.savefig('figures/{}.pdf'.format(outputbase))



    return fig

    
############################

def compareNoiseProfiles(data = None):
    ''' Compare intrinsic noise levels to assumed shape noise'''

    if data is None:
        data = {}

    if 'centers_mpc' not in data:

        config = nfwfit.readConfiguration('/vol/euclid1/euclid1_1/dapple/mxxl_lensing/mxxlsnap54/general-c4-r10-n0_0-xrayNONE/config.sh')

        simreader = nfwfit.buildSimReader(config)

        nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

        fitter = nfwfit.buildFitter(config)

        intrnoise_profiles = []

        for haloid in range(800, 880):

            catbase = '/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap54/halo_54_{}_0'.format(haloid)

            catalog = nfwfit.readSimCatalog(catbase, simreader, config)

            r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens = fitter.prepData(catalog)

            intrnoise_profiles.append(sigma_ghat)

        intrnoise_profiles = np.row_stack(intrnoise_profiles)
        intrnoise = np.mean(intrnoise_profiles, axis=0)
        intrnoise_err = np.std(intrnoise_profiles, axis=0)

        edges_mpc = np.linspace(config.profilemin, config.profilemax, config.nbins+1)
        centers_mpc = (edges_mpc[1:] + edges_mpc[:-1])/2.
        dL = nfwutils.global_cosmology.angulardist(zlens)
        edges_arcmin = (edges_mpc/dL)*(180/np.pi)*60
        bin_areas = np.pi*(edges_arcmin[1:]**2 - edges_arcmin[:-1]**2)

        shapesigma = 0.25
        galdensity = 20. # per sq arc min
        shape_noise = shapesigma / np.sqrt(galdensity*bin_areas)

        data['centers_mpc'] = centers_mpc
        data['intrnoise'] = intrnoise
        data['intrnoise_err'] = intrnoise_err
        data['shape_noise'] = shape_noise

    else:

        centers_mpc = data['centers_mpc']
        intrnoise = data['intrnoise']
        intrnoise_err = data['intrnoise_err']
        shape_noise = data['shape_noise']
    

    fig = pylab.figure()
    ax = pylab.gca()

    ax.plot(centers_mpc, intrnoise/intrnoise[2], label='Intrinsic Noise [80 cluster avg]', marker='None', linestyle='-', linewidth = 2, color = pp.colors[0])
    ax.plot(centers_mpc, shape_noise/shape_noise[2], marker='None', linestyle='-', linewidth = 3, color = pp.colors[1], label='Shape Noise')
    
    ax.legend(loc='upper right')
    ax.set_xlabel('Radius [Mpc]', fontsize=16)
    ax.set_ylabel('Relative Shear Error [Arbit Norm]', fontsize=16)

    fig.tight_layout()

    fig.savefig('figures/relative_shear_noise.png')
    fig.savefig('figures/relative_shear_noise.pdf')
    fig.savefig('figures/relative_shear_noise.eps')


    return fig, data

    
