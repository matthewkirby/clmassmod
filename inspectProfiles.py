###################################
# View shear profiles from bcc simulation, to see if I am having fitting problems
######################################

import pylab, cPickle, numpy as np, glob, os
import nfwfit_bcc, nfwmodeltools as tools, nfwutils
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


def plotChisqSurface(simfile, config, fig = None):  


    catalog = nfwfit_bcc.readSimCatalog(simfile, config)

    fitter = nfwfit_bcc.buildFitter(config)

    fitdata = fitter.prepFit(catalog, config)

    r_mpc, ghat, sigma_ghat, beta_s, ebta_s2, z_lens = fitdata

    model = fitter.createNFWModel(beta_s, ebta_s2, z_lens, 1e14)

    masses = np.arange(1e14, 1e15, 1e13)/1e14

    chisqs = np.zeros_like(masses)

    def chisq(m200):
        g_pred = model(r_mpc, m200)
        return  np.sum(((g_pred - ghat)/sigma_ghat)**2)

    for i in range(len(masses)):
        chisqs[i] = chisq(masses[i])

    if fig is None:
        fig = pylab.figure()
    ax = pylab.gca()
    ax.plot(masses, chisqs)

    return fig

    



    

    



    
