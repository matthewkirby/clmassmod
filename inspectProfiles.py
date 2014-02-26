###################################
# View shear profiles from bcc simulation, to see if I am having fitting problems
######################################

import pylab, cPickle, numpy as np, glob, os
import nfwfit, nfwmodeltools as tools, nfwutils, fitmodel
import basicBinning

##########################

def profileFromMass(mcrelation, mdelta, delta, r_mpc, simcat ):


    zcluster = simcat.hdu.header['ZLENS']

    concentration = mcrelation(mdelta, zcluster, delta)
    
    r_scale = nfwutils.rscaleConstM(mdelta, concentration, zcluster, delta)
    
    
    nfw_shear_inf = tools.NFWShear(r_mpc, concentration, r_scale, zcluster)
    nfw_kappa_inf = tools.NFWKappa(r_mpc, concentration, r_scale, zcluster)

    beta_s  = np.mean(simcat['beta_s'])
    beta_s2 = np.mean(simcat['beta_s']**2)

    
    g = beta_s*nfw_shear_inf / (1 - ((beta_s2/beta_s)*nfw_kappa_inf) )
                
    return g

###########################################
    


def inspectProfile(simfile, outfile, config):

    simcat = nfwfit_bcc.readSimCatalog(simfile, config)
    zcluster = simcat.hdu.header['ZLENS']

    fitter = nfwfit_bcc.buildFitter(config)
    bootstrapbinner = basicBinning.bootstrapequalbins(config)

    r_mpc, ghat, sigma_ghat = fitter.profileBuilder(simcat, config)

    alt_r_mpc, alt_ghat, alt_sigma_ghat = bootstrapbinner(simcat, config)

    m200s, nfails = cPickle.load(open(outfile))

    fitM200 = np.median(m200s)

    profile_r = np.logspace(np.log10(r_mpc[0]), np.log10(r_mpc[-1]), 250)

    profile_g = profileFromMass(fitter.massconRelation, fitM200, 200, profile_r, simcat)

    reference_g = profileFromMass(fitter.massconRelation, fitM200, 200, r_mpc, simcat)

    chisq = np.sum(((ghat - reference_g)/sigma_ghat)**2)

    fig = pylab.figure()
    pylab.errorbar(alt_r_mpc, alt_ghat, alt_sigma_ghat, fmt='gs')
    pylab.errorbar(r_mpc, ghat, sigma_ghat, fmt='bo')

    pylab.plot(profile_r, profile_g, 'r-')
#    pylab.axhline(1.0, c='r', linewidth=2)
    pylab.xlabel('%f %d' % (chisq, len(r_mpc) - 2))

    return fig



###############################################


def stackProfiles(sims, outdir):

    configfile = '%s/config.sh' % outdir

    config = nfwfit_bcc.readConfiguration(configfile)

    fitter = nfwfit_bcc.buildFitter(config)

    radial_points = []
    amp_points = []

    stackedcat = None

    for simfile in sims:

        base = os.path.basename(simfile)
        root, ext = os.path.splitext(base)

        outfile = '%s/%s.out' % (outdir, root)
        
        simcat = nfwfit_bcc.readSimCatalog(simfile, config)
        if stackedcat is None:
            stackedcat = simcat
        else:
            stackedcat = stackedcat.append(simcat)

    print len(stackedcat)
    r_mpc, ghat, sigma_ghat = fitter.profileBuilder(stackedcat, config)

    return r_mpc, ghat, sigma_ghat, stackedcat
    


##############################################


def stackResiduals(sims, outdir):

    configfile = '%s/config.sh' % outdir

    config = nfwfit_bcc.readConfiguration(configfile)

    fitter = nfwfit_bcc.buildFitter(config)

    radial_points = []
    amp_points = []

    for simfile in sims:

        base = os.path.basename(simfile)
        root, ext = os.path.splitext(base)

        outfile = '%s/%s.out' % (outdir, root)
        

        simcat = nfwfit_bcc.readSimCatalog(simfile, config)

        r_mpc, ghat, sigma_ghat = fitter.profileBuilder(simcat, config)

        m200s, nfails = cPickle.load(open(outfile))

        fitM200 = np.median(m200s)

        trueM200 = simcat.hdu.header['M200C']
        zcluster = simcat.hdu.header['ZLENS']

        concentration = fitter.massconRelation(trueM200, zcluster, 200)
    
        r_scale = nfwutils.rscaleConstM(trueM200, concentration, zcluster, 200)
    

        profile_g = profileFromMass(fitter.massconRelation, fitM200, 200, r_mpc, simcat)

#        residuals = profile_g - ghat

#        scaled_r = r_mpc / r_scale
        ratio = profile_g / ghat

#        scaled_residuals = residuals / (r_scale*nfwutils.deltaC(concentration)*nfwutils.global_cosmology.angulardist(zcluster) * nfwutils.global_cosmology.beta([1e6], zcluster) * nfwutils.global_cosmology.hubble2(zcluster))

        radial_points.extend(r_mpc)
        amp_points.extend(ratio)

    return radial_points, amp_points



######################


def plotLogProbSurface(catalogname, configfile, fig = None, noiselevels=[0., 0.03, 0.07, 0.15]):  

    if fig is None:
        fig = pylab.figure()



    config = nfwfit.readConfiguration(configfile)

    simreader = nfwfit.buildSimReader(config)

    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    catalog = nfwfit.readSimCatalog(catalogname, simreader, config)

    fitter = nfwfit.buildFitter(config)

    r_mpc, ghat, sigma_ghat, beta_s, beta_s2, zlens = fitter.prepData(catalog, config)

    print sigma_ghat

    fitter.model.setData(beta_s, beta_s2, zlens)

    guess = fitter.model.guess()

    massgrid = np.arange(1e14, 1e15, 2.5e13)/fitter.model.massScale
    masscenters = (massgrid[1:] + massgrid[:-1])/2.

    for noise in noiselevels:

        noisy_g = ghat + noise*np.random.standard_normal(len(ghat))
        noisy_sigma = np.sqrt(sigma_ghat**2 + noise**2)


        if len(guess) == 2:

            concengrid = np.arange(1., 15., 0.25)
            concencenters = (concengrid[1:] + concengrid[:-1])/2.

            chisqgrid = np.zeros((len(masscenters), len(concencenters)))

            for i in range(len(masscenters)):
                for j in range(len(concencenters)):

                    chisqgrid[i,j] = fitmodel.ChiSqStat(noisy_g, 
                                                        noisy_sigma, 
                                                        fitter.model(r_mpc, masscenters[i], concencenters[j]))

            probgrid = np.exp(-0.5*(chisqgrid - np.max(chisqgrid)))
            massprobgrid = np.sum(probgrid, axis=1) / np.sum(probgrid)
            logprob = np.log(massprobgrid)



        else:

            chisqgrid = np.zeros(len(masscenters))

            for i in range(len(massgrid)):

                chisqgrid[i] = fitmodel.ChiSqStat(noisy_g, noisy_sigma, fitter.model(r_mpc, masscenters[i]))

            logprob = -0.5*chisqgrid




        ax = pylab.gca()
        ax.plot(masscenters, np.exp(logprob - np.max(logprob)), label=noise)

        print 'Max: {0}'.format(masscenters[np.argmax(logprob)])

    pylab.legend()
    return fig, massgrid, concengrid, chisqgrid, logprob

    



    

    



    
