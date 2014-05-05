# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Current Status of N-Body Cluster WL Calibration

# <markdowncell>

# Multiple studies using ray-traced n-body simulations have already looked at the bias in weak lensing (WL) measurements of cluster. My goal is to apply exactly the same methods to simulations as have been used on existing WL cluster studies. I will then predict the expected bias for those studies and methods, while also optimizing those methods for LSST.
# 
# Currently, I am using an NFW-halo model fit using $\chi^2$ statistics from shear profiles. I also have access to Henk Hoekstra's mass-aperture code. 
# 
# I am using the BCC simulations and the Becker & Kratsov simulations.

# <headingcell level=2>

# BCC Simulation Results with NFW Halo Fitting

# <markdowncell>

# I can use the BCC simulations to measure mass bias as a function of mass and redshift. Currently only $M_{200}$ values are available. I am assuming that the ray-tracing and halo-mass measurements in the simulation are valid over all masses and redshifts.
# 
# I select galaxies behind each halo to avoid dealing with cluster and foreground contamination. The background is characterized by $<\beta>$ and $<\beta^2>$. Note that cluster and foreground contamination can be included, at least partially, to test correction procedures. Note also that the BK simulations have a single plane of "background" at z=1.
# 
# Results approximately agree with Joerg Dietrich for NFW fits with a Duffy mass-concentration relation, giving me confidence that I am interpreting the simulations correctly.

# <headingcell level=2>

# Setup

# <codecell>

from pylab import *
import cPickle, astropy.io.fits as pyfits, ldac
from mpl_toolkits.mplot3d import Axes3D

# <codecell>

def loadData(simfile, h=None, useM200=True):

    results = cPickle.load(open(simfile))


    medians = results['measured_m200s']
    actuals = results['true_m200s']




    ratio = medians/actuals

    clean = medians > 0
    
    
    data = {'medians' : medians[clean],
            'sigmas' : 0.01*np.ones_like(medians), 
            'actuals' : actuals[clean],
            'redshifts' : results['redshifts'][clean],
            'ratio' : ratio[clean],
            'logratio' : np.log(ratio[clean])}
    
    return data

# <codecell>

#c4_inner500kpc = loadData('massive_c4_inner500kpc.pkl', h=0.72) # 500kpc - 3Mpc, c=4
#c4_inner750kpc = loadData('massive_c4_inner750kpc.pkl', h=0.72) # 750kpc - 3Mpc, c=4
#duffy_inner500kpc = loadData('massive_duffy_inner500kpc.pkl', h=0.72) # 500kpc - 3Mpc, duffy08 m-c relation
#henk_massapp = loadData('henk_massive.pkl') # Hoekstra's mass aperture 
#becker_c4_inner750kpc = loadData('becker_wtg_baseline_inner750.pkl')  #B&K simulations, c=4, 750kpc-3Mpc; h is already factored out

# <codecell>

def summary1D(simdata, selection=None):
    if selection is None:
        selection = np.ones_like(simdata['ratio']) == 1
    print 'Median Ratio: ', np.median(simdata['ratio'][selection])
    print 'Mean Log Ratio: ', np.exp(np.mean(simdata['logratio'][selection]))
    print 'Log Scatter: ', np.std(simdata['logratio'][selection])
    hist(simdata['logratio'][selection], bins=50)
    axvline(0.0, c='k')
    axvline(np.median(simdata['logratio'][selection]), c='r')
    

# <codecell>

def summary2DRedshift(simdata):
    medians, sigmas, actuals, redshifts, ratio, logratio = simdata['medians'], simdata['sigmas'], simdata['actuals'], simdata['redshifts'], simdata['ratio'], simdata['logratio']
    redshiftbin = []
    ratiobin = []
    errup = []
    errdown = []
    median_err = []
    redshiftsorted = np.argsort(redshifts)
    for i in range(0, len(redshifts), 100):
        maxtake = min(len(redshifts), i+100)
        redshiftbin.append(np.median(redshifts[redshiftsorted][i:maxtake]))
        localratios = np.sort(logratio[redshiftsorted][i:maxtake])
        localMedian = localratios[int(0.5*len(localratios))]
        r68up = localratios[int(0.84*len(localratios))] - localMedian
        r68down = localMedian - localratios[int(0.16*len(localratios))]
        ratiobin.append(np.median(localratios))
        errup.append(r68up)
        errdown.append(r68down)
        median_err.append(bootstrapMedianErr(localratios))
#    subplot(1,2,1)
#    hexbin(redshifts, logratio, gridsize=50, bins='log')
#    errorbar(redshiftbin, ratiobin, np.row_stack([errdown, errup]), fmt='ko', markersize=5, linewidth=3)
#    axhline(0.0, c='r', linewidth=3)
#    xlabel('Redshift')
#    ylabel('Log 10 Ratio')
#    subplot(1,2,2)
#    plot(redshiftbin, np.exp(np.array(ratiobin)))

    return redshiftbin, np.exp(ratiobin), median_err

# <codecell>

def summary2DMass(simdata, selection = None, axisrange=None):
    if selection is None:
        selection = np.ones_like(simdata['ratio']) == 1
    medians, sigmas, actuals, redshifts, ratio, logratio = [x[selection] for x in [simdata['medians'], simdata['sigmas'], simdata['actuals'], simdata['redshifts'], simdata['ratio'], simdata['logratio']]]

#    ngals = len(medians) / 6 
    ngals = 600
    

    massbin = []
    ratiobin = []
    errup = []
    errdown = []
    median_err = []
    masssorted = np.argsort(actuals)
    sortedlog10actuals = np.log10(actuals)[masssorted]
    
    i=0
    while i < len(actuals):

#        if len(actuals) - i < 10:
#            maxtake = len(actuals)
#        else:
#
        maxtake = min(len(actuals), i+ngals)
            

#            binbound = sortedlog10actuals[i:maxtake] - sortedlog10actuals[i] < 0.2
#            maxtake = np.arange(maxtake)[binbound][-1] + i
#            


        massbin.append(np.median(np.log10(actuals[masssorted][i:maxtake])))
        localratios = np.sort(logratio[masssorted][i:maxtake])
        localMedian = localratios[int(0.5*len(localratios))]
        r68up = localratios[int(0.84*len(localratios))] - localMedian
        r68down = localMedian - localratios[int(0.16*len(localratios))]
        ratiobin.append(np.median(localratios))
        errup.append(r68up)
        errdown.append(r68down)
        median_err.append(bootstrapMedianErr(localratios))

        i = maxtake

#    subplot(1,2,1)
#    hexbin(np.log10(actuals), logratio, gridsize=50, bins='log', extent=axisrange)
#    xlabel('True Mass (M200)')
#    ylabel('Log10 Ratio')
#    errorbar(massbin, ratiobin, np.row_stack([errdown, errup]), fmt='ko', markersize=5, linewidth=3)
#    axhline(0.0, c='r', linewidth=3)
#    subplot(1,2,2)
#    plot(massbin, np.exp(np.array(ratiobin)))

    return massbin, np.exp(ratiobin), median_err
#
# <codecell>

def summary3D(simdata):
    medians, sigmas, actuals, redshifts, ratio, logratio = simdata['medians'], simdata['sigmas'], simdata['actuals'], simdata['redshifts'], simdata['ratio'], simdata['logratio']
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(redshifts, actuals, logratio)
    ax.set_xlabel('Redshift')
    ax.set_ylabel('True Mass M200')
    ax.set_zlabel('Log10 Ratio')

# <codecell>

def bootstrapMedianErr(aset, nboot=1000):
    
    bootvals = np.zeros(nboot)
    for i in range(nboot):
        curboot = np.random.randint(0, len(aset), len(aset))
        bootvals[i] = np.median(aset[curboot])
        
    return np.std(bootvals)

# <codecell>

def scatterSummary(simdata, selection = None):

    if selection is None:
        selection = np.ones_like(simdata['ratio']) == 1
    medians, sigmas, actuals, redshifts, ratio, logratio = [x[selection] for x in [simdata['medians'], simdata['sigmas'], simdata['actuals'], simdata['redshifts'], simdata['ratio'], simdata['logratio']]]
    massbin = []
    massbin_min = []
    massbin_max = []
    intrinsic_scatter = []
    logratiodists = []
    sigmadists = []
    median_err = []
    ratiobin = []
    errup = []
    errdown = []
    masssorted = np.argsort(actuals)
    for i in range(0, len(actuals), 600):
        maxtake = min(len(actuals), i+200)
        massbin.append(np.median(np.log10(actuals[masssorted][i:maxtake])))
        massbin_min.append(np.min(np.log10(actuals[masssorted][i:maxtake])))
        massbin_max.append(np.max(np.log10(actuals[masssorted][i:maxtake])))
        
        sigmadists.append((sigmas/medians)[masssorted][i:maxtake])
        
        localratios = np.sort(logratio[masssorted][i:maxtake])
        logratiodists.append(localratios)
     
        localMedian = localratios[int(0.5*len(localratios))]
        r68up = localratios[int(0.84*len(localratios))] - localMedian
        r68down = localMedian - localratios[int(0.16*len(localratios))]
        ratiobin.append(np.median(localratios))
        errup.append(r68up)
        errdown.append(r68down)

        
        intrinsic_scatter.append(np.std(localratios))
        
        median_err.append(bootstrapMedianErr(localratios))
        
    sym_distro = (np.array(errup) + np.array(errdown))/2.
        
    
    return massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists

# <codecell>

bright_and_modz = lambda simdata: np.logical_and(np.log10(simdata['actuals']) > 14.0, simdata['redshifts'] < 0.6)

# <codecell>

gauss = lambda x, mu, sig: np.exp(-0.5*(x-mu)**2/sig**2)/np.sqrt(2*np.pi*sig**2)

# <headingcell level=2>

# What Halos are Available?

# <codecell>
#
#print "Number of Halos:", len(c4_inner500kpc['redshifts'])
#hexbin(c4_inner500kpc['redshifts'], np.log10(c4_inner500kpc['actuals']), gridsize=50)
#xlabel('Redshift')
#ylabel('Log10 Mass')
#
## <codecell>
#
#subplot(1,2,1)
#hist(c4_inner500kpc['redshifts'], bins=50)
#xlabel('Redshift')
#subplot(1,2,2)
#hist(c4_inner500kpc['actuals'], bins=50, log=True)
#xlabel('Mass $M_{200}$')
#
## <markdowncell>
#
## I've selected ~2000 halos from the BCC that are above a mass of $M_{200} = 10^{14.5} M_{\odot}$, spanning $0.2 < z < 1.5$.
#
## <codecell>
#
#print "Number of Halos:", len(becker_c4_inner750kpc['redshifts'])
#hist(becker_c4_inner750kpc['actuals'], bins=50, log=True)
#xlabel('Mass $M_{500}$')
#
## <markdowncell>
#
## The B&K simulations have two redshift slices, at z=0.25, and z=0.5. For the moment, I've plotted just the mass distro from z=0.25. There are only ~730 clusters. 
#
## <headingcell level=2>
#
## Baseline Ratio
#
## <markdowncell>
#
## There is clearly an overall bias reported in the BCC simulations. There is also clearly a mass and redshift dependence to the bias.
#
## <codecell>
#
#summary1D(c4_inner500kpc)
#
## <codecell>
#
#summary1D(c4_inner500kpc, selection = bright_and_modz(c4_inner500kpc))
#
## <markdowncell>
#
## This is best seen in 2D:
#
## <codecell>
#
#summary2DRedshift(c4_inner500kpc)
#
## <codecell>
#
#summary2DMass(c4_inner500kpc)
#
## <codecell>
#
#summary2DMass(c4_inner500kpc, selection = bright_and_modz(c4_inner500kpc), axisrange=(14.5, 15.2, -1, 1))
#
## <markdowncell>
#
## Notice that elimination of low mass and high redshift halos cleans up the noise considerably. Most of the low-scattering tail at low mass is at high redshift. 
## 
## Also keep in mind that this is M200, not M500, and fitting in to 500kpc. 
#
## <headingcell level=2>
#
## Comparing Mass-Concentration Relations
#
## <headingcell level=3>
#
## Summary Plots for Duffy08 Mass-Concentration Relation
#
## <codecell>
#
#summary1D(duffy_inner500kpc)
#
## <codecell>
#
#summary1D(duffy_inner500kpc, selection = bright_and_modz(duffy_inner500kpc))
#
## <codecell>
#
#summary2DRedshift(duffy_inner500kpc)
#
## <codecell>
#
#summary2DMass(duffy_inner500kpc)
#
## <codecell>
#
#summary2DMass(duffy_inner500kpc, bright_and_modz(duffy_inner500kpc), axisrange=(14.5, 15.2, -1, 1))
#
## <markdowncell>
#
## The Duffy08 Mass-con relation appears to perform slightly better than c=4. From above, for comparison is the c=4 plot again:
#
## <codecell>
#
#summary2DMass(c4_inner500kpc, selection = bright_and_modz(c4_inner500kpc), axisrange=(14.5, 15.2, -1, 1))
#
## <headingcell level=3>
#
## Comparing Inner Fit Radii
#
## <markdowncell>
#
## If r=500kpc is too aggressive, does fitting only to r=750kpc do better?
#
## <codecell>
#
#summary2DMass(c4_inner750kpc, bright_and_modz(c4_inner750kpc), axisrange=(14.5, 15.2, -1, 1))
#
## <markdowncell>
#
## Yes, the situation appears to improve, but the bias is still slightly larger (in aggregate) for c=4. However, c=4 does not appear to induce as steep a mass-dependent trend when fitting to 750kpc.
#
## <headingcell level=3>
#
## Mass Apertures and the BK simulations
#
## <codecell>
#
#summary2DMass(henk_massapp)
#
## <codecell>
#
#summary2DMass(becker_c4_inner750kpc)
#
## <markdowncell>
#
## I am still struggling to understand the inputs and outputs for the BK simulations and Henk's mass aperture code. Results currently are obviously off.
#
## <headingcell level=2>
#
## Questions & Tasks
#
## <markdowncell>
#
## * Bootstrap error bars for summary bias plots
## * Get BK sims working
## * Get Henk's Mass-App working
## * What are the centers defined in the halo catalog?
## * Measure P(M_true | M_meas)
## * Measure mass bias at M500.
## * Match mass & redshift distributions of WtG, CCCP, SPT-Magellan, SPT-HST
## * What fraction of the scatter is due to measurement uncertainty, versus intrinsic scatter?
## * What is the scatter shape in mass bins, and how well do we measure those shape parameters?
## * Where are we resolution limited -- what are the limits of the current simulation?
## * What happens when we add foreground contamination? Cluster contamination?
#
## <headingcell level=2>

# Scatter

# <codecell>

#massbin, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = scatterSummary(c4_inner500kpc, selection = bright_and_modz(c4_inner500kpc))
#plot(massbin, sym_distro, label='Central 68%')
#plot(massbin, intrinsic_scatter, label='Intrinsic Scatter - Stddev')
#plot(massbin, median_err, label='Median Err')
#plot(massbin, [np.median(x) for x in sigmadists], label='Typical Errorbar')
#legend()
#
# <markdowncell>

# Looks like errors on individual cluster fits are insignificant compared to the intrinsic scatter. So at least for this mass range, we don't need more sample points per cluster.

# <codecell>

#nplots = len(logratiodists)
#for i in range(nplots):
#    subplot(nplots, 1, i+1)
#    bins, edges, patches = hist(logratiodists[i]/intrinsic_scatter[i], bins=50, normed=True, range=(-5,5))
#    plot(edges, gauss(edges, 0., 1.0), 'r-')
#    print intrinsic_scatter[i], sym_distro[i]
#
## <markdowncell>

# Looks like our error distributions are kind of non-gaussian, with extended tails. That would account for the difference between the 68% range and the standard error for the width of the distributions. That also means we will need to characterize the extended tails if we want to do this accurately.

# <codecell>


