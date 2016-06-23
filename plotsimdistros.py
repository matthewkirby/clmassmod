import publication_plots as pp
import pylab
import numpy as np
import cPickle
import deconvolvedlognorm as dln
import pymc
import load_chains, os, glob
import confidenceinterval as ci
import readtxtfile
import nfwutils
import rundln
import astropy.io.ascii as asciireader
import sys



#############

def bootstrapMean(sample, nboots=1000):

    nsamples = len(sample)

    bootedmeans = np.zeros(nboots)
    for i in range(nboots):
        curboot = np.random.randint(0, nsamples, nsamples)
        bootedmeans[i] = np.mean(sample[curboot])

    return bootedmeans


############

def bootstrapMedian(sample, nboots=1000):

    nsamples = len(sample)

    bootedmedians = np.zeros(nboots)
    for i in range(nboots):
        curboot = np.random.randint(0, nsamples, nsamples)
        bootedmedians[i] = np.median(sample[curboot])

    return bootedmedians




#############


def bootstrapStd(sample, nboots=1000):

    nsamples = len(sample)

    bootedstds = np.zeros(nboots)
    for i in range(nboots):
        curboot = np.random.randint(0, nsamples, nsamples)
        bootedstds[i] = np.std(sample[curboot])

    return bootedstds

#############

def createErrorbars(samples):

    nsamples = len(samples)
    
    centers = np.zeros(nsamples)
    errs = np.zeros((2,nsamples))

    for i, sample in enumerate(samples):

        sortedsample = np.sort(sample)
        nentries = len(sample)
        
        low = sortedsample[int(0.16*nentries)]
        med = sortedsample[int(0.5*nentries)]
        high = sortedsample[int(0.84*nentries)]

        centers[i] = med
        errs[0,i] = med - low
        errs[1,i] = high - med

    return centers, errs


#############

c = [(.9,.6,0), (.35, .7, .9), (0,.6,.5), 
     (0, .45, .7), (.8, .4, 0), (.8, .6, .7), (0.95, 0.9, 0.25)]

############

def fitLogNormDistro(truemass, measuredmass, measuredmasserr, massedges, meanax, stdax, colorindex):

    log10massedges = np.log10(massedges)

    nbins = len(log10massedges) - 1

    ylows = []
    yhighs = []
    xpoints = []    

    ystdlows = []
    ystdhighs = []

    for i in range(nbins):

        inbin = np.logical_and(truemass >= massedges[i],
                               truemass < massedges[i+1])


        if len(truemass[inbin]) < 25:
            continue

        xpoints.append(massedges[i])
        xpoints.append(massedges[i+1])

        print len(measuredmass[inbin]), len(measuredmasserr[inbin]), len(truemass[inbin])
        parts = None
        for i in range(20):
            try:
                parts = dln.buildModel(measuredmass[inbin], measuredmasserr[inbin], truemass[inbin])
                break
            except pymc.ZeroProbability:
                continue
        if parts is None:
            raise pymc.ZeroProbability
        (logmu, logmuerr), (logsigma, logsigmaerr) = dln.runFit(parts)
        

        mu_low = np.exp(logmu[0] - logmuerr[0,0])
        mu_high = np.exp(logmu[0] + logmuerr[0,0])
        std_low = np.exp(logsigma[0] - logsigmaerr[0,0])
        std_high = np.exp(logsigma[0] + logsigmaerr[0,0])
        
        ylows.append( mu_low)
        ylows.append( mu_low)
        yhighs.append(mu_high)
        yhighs.append(mu_high)

        ystdlows.append(std_low)
        ystdlows.append(std_low)
        ystdhighs.append(std_high)
        ystdhighs.append(std_high)
                     


    meanax.fill_between(xpoints, ylows, yhighs, alpha=0.8, color = c[colorindex], hatch = None)
    stdax.fill_between(xpoints, ystdlows, ystdhighs, alpha=0.8, color = c[colorindex], hatch = None)
    patch = pylab.Rectangle((0, 0), 1, 1, fc=c[colorindex], alpha=0.8, hatch = None)

    return patch
        

################

def gatherChainFiles(chaindir, delta, binnum = None):

    if binnum is None:

        #may 2015 style
        chainfiles = glob.glob('%s/rundln*.%d.[0-9]' % (chaindir, delta))

        if len(chainfiles) == 0:
            chainfiles = glob.glob('%s/dln_*.%d.chain.0' % (chaindir, delta))

        if len(chainfiles) == 0:
            chainfiles = glob.glob('%s/dln_*.chain.0' % chaindir)

    else:
        chainfiles = glob.glob('%s/rundln*.%d.%d' % (chaindir, delta, binnum))
        


#    print chainfiles


    return sorted(chainfiles)

################

def weightedaverage(means, errs):

    weights = 1./errs**2
    totalweight = np.sum(weights)
    mu = np.sum(weights*means)/totalweight
    sig = np.sqrt(1./totalweight)

    return mu, sig

def precomputedLogNormDistro(chaindir, delta, meanax, stdax, colorindex, alpha=0.8, biaslabel = True, xoffset = 1.0, binnum=None):



    ylows = []
    yhighs = []
    xpoints = []    

    ave = []
    aveerr = []
    stdev = []
    stdeverr = []

    ystdlows = []
    ystdhighs = []

    chainfiles = gatherChainFiles(chaindir, delta, binnum = binnum)

    assert(len(chainfiles) > 0)

    for chainfile in chainfiles:

        fileroot = chainfile.split('.chain')[0]


            

        chain = load_chains.loadChains([chainfile], trim=True)

        print chainfile, len(chain['logmu'])
        if len(chain['logmu'][0,:]) < 5000:
            print 'Skipping'
            continue

        split = int((chain['logmu'].shape[1] + 1000)/2.)
        splitlen = split - 1000
        c1mean = np.mean(chain['logmu'][0,1000:split])
        c1err = np.std(chain['logmu'][0,1000:split])
        c2mean = np.mean(chain['logmu'][0,split:])
        c2err = np.std(chain['logmu'][0,split:])
        assert(np.abs(c1mean - c2mean)/np.sqrt(c1err**2 + c2err**2) < 3.)


        massbinlow, massbinhigh = [x[0] for x in readtxtfile.readtxtfile('%s.massrange' % fileroot)]

        xpoints.append(massbinlow)
        xpoints.append(massbinhigh)

        mu, muerr = ci.maxDensityConfidenceRegion(np.exp(chain['logmu'][0,1000::3]))
        sig, sigerr = ci.maxDensityConfidenceRegion(np.exp(chain['logsigma'][0,1000::3]))

        ave.append(np.mean(np.exp(chain['logmu'][0,1000::3])))
        aveerr.append(np.std(np.exp(chain['logmu'][0,1000::3])))
        stdev.append(np.mean(np.exp(chain['logsigma'][0,1000::3])))
        stdeverr.append(np.std(np.exp(chain['logsigma'][0,1000::3])))


        
        mu_low =  mu - muerr[0]
        mu_high = mu + muerr[1]
        std_low = sig - sigerr[0]
        std_high =sig + sigerr[1]
        
        ylows.append( mu_low)
        ylows.append( mu_low)
        yhighs.append(mu_high)
        yhighs.append(mu_high)

        ystdlows.append(std_low)
        ystdlows.append(std_low)
        ystdhighs.append(std_high)
        ystdhighs.append(std_high)

        x_center = xoffset*(massbinlow + massbinhigh)/2.
        mu_center = (mu_high + mu_low)/2.
        mu_err = (mu_high - mu_low)/2.
        std_center = (std_high + std_low)/2.
        std_err = (std_high - std_low)/2.

        print mu_center, mu_err


        meanax.errorbar([x_center], [mu_center], [mu_err], [[x_center - massbinlow], [massbinhigh - x_center]], color = c[colorindex], marker='None', linestyle='None', elinewidth=2.)
        stdax.errorbar([x_center], [std_center], [std_err], [[x_center - massbinlow], [massbinhigh - x_center]], color = c[colorindex], marker='None', linestyle='None', elinewidth=2.)


#        meanax.fill_between([massbinlow, massbinhigh], 
#                            [mu_low, mu_low], 
#                            [mu_high, mu_high], 
#                            alpha=alpha, color = c[colorindex], hatch = None)
#        stdax.fill_between([massbinlow, massbinhigh],
#                           [std_low, std_low],
#                           [std_high, std_high],
#                           alpha=alpha, color = c[colorindex], hatch = None)
# 


                     
    if len(xpoints) == 0:
        return None

    ave = np.array(ave)
    aveerr = np.array(aveerr)
    stdev = np.array(stdev)
    stdeverr = np.array(stdeverr)

    summary = weightedaverage(ave, aveerr), weightedaverage(stdev, stdeverr)

    if biaslabel is True:
        meanax.text(2.5e14, 0.75 + float(colorindex)/10., '%1.2f +/- %1.2f' % (summary[0][0], summary[0][1]))

    patch = pylab.Rectangle((0, 0), 1, 1, fc=c[colorindex], alpha=alpha, hatch = None)

    return patch, summary
        





################    


def plotLogNormDistro(truemass, measuredmass, massedges, meanax, nongaussax, stdax, label, colorindex, useLog = True):

    log10massedges = np.log10(massedges)

    log10centers = (log10massedges[:-1] + log10massedges[1:])/2.
    nbins = len(log10centers)




    ratio = measuredmass / truemass

    centers = []
    medians = []
    nongausses = []
    stds = []
    ylows = []
    yhighs = []
    xpoints = []    

    for i in range(len(log10centers)):

        inbin = ratio[np.logical_and(truemass >= massedges[i],
                               truemass < massedges[i+1])]

        if len(inbin) < 25:
            continue

        xpoints.append(massedges[i])
        xpoints.append(massedges[i+1])

        centers.append(10**(log10centers[i]))

        if (inbin < 0).any() and useLog is True:
            print 'ILLEGAL'


        if useLog:
            logratio = np.log(inbin)

            medians.append(bootstrapMean(logratio))
            nongausses.append(np.median(logratio) - np.mean(logratio))
            stds.append(bootstrapStd(logratio))

        else:

            medians.append(bootstrapMean(inbin))
            nongausses.append(np.median(inbin) - np.mean(inbin))
            stds.append(bootstrapStd(inbin))

    centers = np.array(centers)



    mediancenter, medianerrs = createErrorbars(medians)

    print mediancenter, medianerrs

    if useLog is True:
        for i in range(len(mediancenter)):
            ylows.append( np.exp(mediancenter[i] - medianerrs[0,i]))
            ylows.append( np.exp(mediancenter[i] - medianerrs[0,i]))
            yhighs.append(np.exp(mediancenter[i] + medianerrs[1,i]))
            yhighs.append(np.exp(mediancenter[i] + medianerrs[1,i]))
    else:
        for i in range(len(mediancenter)):
            ylows.append( mediancenter[i] - medianerrs[0,i])
            ylows.append( mediancenter[i] - medianerrs[0,i])
            yhighs.append(mediancenter[i] + medianerrs[1,i])
            yhighs.append(mediancenter[i] + medianerrs[1,i])


#    meanax.errorbar(centers + offset, meancenter, meanerrs, **plotargs)
    print len(xpoints), 
    meanax.fill_between(xpoints, ylows, yhighs, alpha=0.8, color = c[colorindex], label = label, hatch = None)
    patch = pylab.Rectangle((0, 0), 1, 1, fc=c[colorindex], alpha=0.8, hatch = None)
                 


    print nongausses

#    nongaussax.plot(centers, nongausses, marker='o', **plotargs)


    stdcenter, stderrs = createErrorbars(stds)
    stdax.errorbar(centers, stdcenter, stderrs, label = label, color = c[colorindex])

    return patch


################################


def plotRadiusMXXL():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)
    stdsfig = pylab.figure()
    stdsax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    
    radialranges = [5,6,8,9]
    radialnames = ['0.5 - 1.5', '0.5 - 2.5', '0.75 - 1.5', '0.75 - 2.5']
    offsets = np.arange(-1.5e13, 2.0e13, 1e13)



    for i, radrange in enumerate(radialranges):
        with open('run7consolidated/mxxlsnap41.c4-r%d-n0_0_corenone.pkl' % radrange, 'rb') as input:

            consol = cPickle.load(input)

            plotLogNormDistro(consol['true_m500s'], 
                              consol['measured_m500s'],
                              massedges,
                              meansax,
                              stdsax,
                              offset = offsets[i],
                              label = radialnames[i],
                              linestyle='None',
                              linewidth=2.)

    meansax.set_xscale('log')
    meansax.set_xlabel('Mass', fontsize=16)
    meansax.set_ylabel('Mean Log-Bias', fontsize=16)
    meansax.legend(loc='lower left')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('radiusmxxl_mean.png')

    stdsax.set_xscale('log')
    stdsax.set_xlabel('Mass', fontsize=16)
    stdsax.set_ylabel('Standard Deviation Log-Bias', fontsize=16)
    stdsax.legend()
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('radiusmxxl_std.png')

    return meansfig, stdsfig


#####################################

 
def plotNoiseMXXL():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    

#    chaindirs = [#'mxxl_imperial/mxxlsnap41/mcmc_linear-c4-r5-n0_0-corenone-linearbins12',
#                 '/users/dapple/astro/mxxlsims/mxxl_imperial/mxxlsnap41/mcmc_linear-c4-r5-n2_2-corenone-lineargaussbins12',
#                 '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-testprofile',
#                 #'mxxl_imperial/mxxlsnap41/mcmc_linear-c4-r5-n6_4-corenone-lineargaussbins12',
#                 #'mxxl_imperial/mxxlsnap41/mcmc_linear-c4-r5-n4_3-corenone-lineargaussbins12',
#        
#    ]

#    noisenames = [#'No Noise', 
#                  '20 gals/sq. arcmin $\sigma_e = 0.33$',
#                  'hst noise control',
#                  #'10 gals/sq. arcmin $\sigma_e = 0.4$',
#                  #'4 gals/sq. arcmin $\sigma_e = 0.5$'
#    ]
#


    chaindirs = ['/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-simple',
                 '/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-lownoise',
                 '/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0000-5748',
                 '/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-highnoise']

    noisenames = ['Simple Profile',
                  'Low Noise',
                  'Actual Noise',
                  'High Noise']



    patches = []
    labels = []


    for i in range(len(chaindirs)-1, -1,-1):

        chaindir = chaindirs[i]

        print chaindir


        label = noisenames[i]

        patch = precomputedLogNormDistro(chaindir, 
                                         massedges,
                                         meansax,
                                         stdax,
                                         colorindex = i)

        if patch is None:
            continue

        patches.append(patch)
        labels.append(label)


    meansax.set_xscale('log')
    meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
    meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    meansax.set_xlim(2e14, 1.3e15)
    meansax.set_ylim(0.65, 1.2)
    meansax.set_xticks([1e15])
    meansax.set_xticklabels(['10'])
    meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    meansax.legend(patches[::-1], labels[::-1], loc='upper left')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('hstnoisemxxl_logmean_control.png')

    stdax.set_xscale('log')
    stdax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
    stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    stdax.set_xlim(2e14, 1.3e15)
#    stdax.set_ylim(0.85, 1.10)
    stdax.set_xticks([1e15])
    stdax.set_xticklabels(['10'])
    stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    stdax.legend(patches[::-1], labels[::-1], loc='upper left')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('hstnoisemxxl_logstd_control.png')


    return meansfig, stdsfig



############################


def plotShearErrEstimateMXXL():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    
    radrange = 5
    radialname = ['0.5 - 1.5']
    noiseranges = ['2_2', '4_3']
    noisenames = ['ng=20  $\sigma_e = 0.33$',
                  'ng=4 $\sigma_e = 0.5$']
    errests = ['', '-gaussianshearerr']
    errestnames = ['bootstrap', 'gaussian scaled']


    patches = []
    labels = []


    for i, errest in enumerate(errests):
        for j, noiserange in enumerate(noiseranges):

            consolfile = 'mxxl_imperial/rundirs/run8consolidated/mxxlsnap41.c4-r%d-n%s-corenone-linearbins12%s.pkl' % (radrange, noiserange, errest)
            print consolfile

            with open(consolfile, 'rb') as input:

                consol = cPickle.load(input)

                label = '%s; %s' % (noisenames[j], errestnames[i])

                patch = fitLogNormDistro(consol['true_m200s'], 
                                         consol['measured_m200s'],
                                         consol['measured_m200errs'],
                                         massedges,
                                         meansax,
                                         stdax,
                                         colorindex = i)

                patches.append(patch)
                labels.append(label)


    meansax.set_xscale('log')
    meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
    meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    meansax.set_xlim(2e14, 1.3e15)
#    meansax.set_ylim(0.85, 1.10)
    meansax.set_xticks([1e15])
    meansax.set_xticklabels(['10'])
    meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    meansax.legend(patches, labels, loc='upper left')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('noisemxxl_logmean.png')

    stdax.set_xscale('log')
    stdax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
    stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    stdax.set_xlim(2e14, 1.3e15)
#    stdax.set_ylim(0.85, 1.10)
    stdax.set_xticks([1e15])
    stdax.set_xticklabels(['10'])
    stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    stdax.legend(patches, labels, loc='upper left')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('noisemxxl_logstd.png')


    return meansfig, stdsfig



############################



def plotBinningMXXL():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    
    radialrange = 5
    radialname = '0.5 - 1.5'
    noiserange = '6_4'
    noisename = ['10 gals/sq. arcmin $\sigma_e = 0.4$']
    binnings = ['linearbins6', 'linearbins12', 'logbins6']
    binningnames = ['linear 6 bins', 'linear 12 bins', 'log 6 bins']


    patches = []
    labels = []


    for i, binning in enumerate(binnings):


            consolfile = 'mxxl_imperial/rundirs/run9consolidated/mxxlsnap41.c4-r%d-n%s-corenone-%s.pkl' % (radialrange, noiserange,binning)
            print consolfile

            with open(consolfile, 'rb') as input:

                consol = cPickle.load(input)

                label = binningnames[i]

                patch = fitLogNormDistro(consol['true_m200s'], 
                                         consol['measured_m200s'],
                                         consol['measured_m200errs'],
                                         massedges,
                                         meansax,
                                         stdax,
                                         colorindex = i)

                patches.append(patch)
                labels.append(label)


    meansax.set_xscale('log')
    meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
    meansax.set_title(r'10 Galaxies/ sq arcmin; $\sigma_e = 0.4$')
    meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    meansax.set_xlim(2e14, 1.3e15)
#    meansax.set_ylim(0.85, 1.10)
    meansax.set_xticks([1e15])
    meansax.set_xticklabels(['10'])
    meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    meansax.legend(patches, labels, loc='upper left')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('binningmxxl_logmean.png')

    stdax.set_xscale('log')
    stdax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
    stdax.set_title(r'10 Galaxies/ sq arcmin; $\sigma_e = 0.4$')
    stdax.set_xlim(2e14, 1.3e15)
#    stdax.set_ylim(0.85, 1.10)
    stdax.set_xticks([1e15])
    stdax.set_xticklabels(['10'])
    stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    stdax.legend(patches, labels, loc='lower left')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('binningmxxl_logstd.png')


    return meansfig, stdsfig



############################



def plotCoreMXXL():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)
    nongaussfig = pylab.figure()
    nongaussax = nongaussfig.add_subplot(1,1,1)
    stdsfig = pylab.figure()
    stdsax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(5e15), 12)
    
    radialranges = [5]
    radialnames = ['0.50 - 1.5']
    coreranges = ['none', '0', '5']
    corenames = ['Exact Centering', r"Typical Miscentering: $\theta_c = 0.25'$",
                 r"Worst Case: $\theta_c = 1.5'$"]

    patches = []
    labels = []


    for i, radrange in enumerate(radialranges):
        for j, corerange in enumerate(coreranges):

            consolfile = 'mxxl_imperial/rundirs/run7consolidated/mxxlsnap41.c4-r%d-n0_0_core%s.pkl' % (radrange, corerange)
            print consolfile

            with open(consolfile, 'rb') as input:

                consol = cPickle.load(input)
                
                label = corenames[j]
                

                patch = plotLogNormDistro(consol['true_m200s'], 
                                  consol['measured_m200s'],
                                  massedges,
                                  meansax,
                                  nongaussax,
                                  stdsax,
                                  label = label,
                                  colorindex = j)

                patches.append(patch)
                labels.append(label)


  
    nongaussax.set_xscale('log')
    nongaussax.set_xlabel('Mass', fontsize=16)
    nongaussax.set_ylabel('Median - Mean', fontsize=16)
    nongaussax.axhline(0.0, c='k', linewidth=2, linestyle='--')
    nongaussax.legend(loc='upper left')
    nongaussfig.canvas.draw()
    nongaussfig.tight_layout()
    nongaussfig.savefig('coremxxl_nongauss.png')



    stdsax.set_xscale('log')
    stdsax.set_xlabel('Mass', fontsize=16)
    stdsax.set_ylabel('Standard Deviation Log-Bias', fontsize=16)
    stdsax.set_ybound(0.15, 0.32)
    stdsax.legend(loc='lower left')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('coremxxl_std.png')


    return meansfig, nongaussfig, stdsfig


#######################################


def plotCoreBK11():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)
    nongaussfig = pylab.figure()
    nongaussax = nongaussfig.add_subplot(1,1,1)
    stdsfig = pylab.figure()
    stdsax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(1e14), np.log10(1e15), 6)
    
    radialranges = [9,6 ]
    radialnames = ['0.75 - 2.5', '0.5 - 2.5']
    coreranges = ['none', '0', '5']
    corenames = ['none', '0.25', '1.5']
    offsets = np.linspace(-2.2e13, 2.2e13, 6)

    colors = ['b', 'g', 'r']
    alphas = [1.0, 0.5]


    for i, radrange in enumerate(radialranges):
        for j, corerange in enumerate(coreranges):

            consolfile = '../rundirs/run7consolidated/bk11snap141..c4-r%d-n4_3_core%s.pkl' % (radrange, corerange)
            print consolfile

            with open(consolfile, 'rb') as input:

                consol = cPickle.load(input)

                

                plotLogNormDistro(consol['true_m500s'], 
                                  consol['measured_m500s'],
                                  massedges,
                                  meansax,
                                  nongaussax,
                                  stdsax,
                                  offset = offsets[3*i+j],
                                  label = '%s %s' % (radialnames[i], corenames[j]),
                                  linestyle='None',
                                  linewidth=2.,
                                  color=colors[j],
                                  alpha = alphas[i])

    meansax.set_xscale('log')
    meansax.set_xlabel('Mass', fontsize=16)
    meansax.set_ylabel('Mean Log-Bias', fontsize=16)
    meansax.axhline(0.0, c='k', linewidth=2, linestyle='--')
    meansax.legend(loc='upper right')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('corebk11_mean.png')

    nongaussax.set_xscale('log')
    nongaussax.set_xlabel('Mass', fontsize=16)
    nongaussax.set_ylabel('Median - Mean', fontsize=16)
    nongaussax.axhline(0.0, c='k', linewidth=2, linestyle='--')
    nongaussax.legend(loc='upper left')
    nongaussfig.canvas.draw()
    nongaussfig.tight_layout()
    nongaussfig.savefig('corebk11_nongauss.png')



    stdsax.set_xscale('log')
    stdsax.set_xlabel('Mass', fontsize=16)
    stdsax.set_ylabel('Standard Deviation Log-Bias', fontsize=16)
    stdsax.set_ybound(0.15, 0.32)
    stdsax.legend(loc='upper right')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('corebk11_std.png')


    return meansfig, nongaussfig, stdsfig


#######################################

    
    
def plotNoiseBK11():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)
    nongaussfig = pylab.figure()
    nongaussax = nongaussfig.add_subplot(1,1,1)
    stdsfig = pylab.figure()
    stdsax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(1e14), np.log10(1e15), 6)
    
    radialranges = [8]
    radialnames = ['0.75 - 1.5']
    noiseranges = ['0_0', '3_2', '4_3']
    noisenames = ['NoNoise', '20-0.33', '7-0.5']
    offsets = np.linspace(-1.5e13, 1.5e13, 3)

    colors = 'b g r'.split()
    alphas = [1.0, 0.5]


    for i, radrange in enumerate(radialranges):
        for j, noiserange in enumerate(noiseranges):

            consolfile = '../rundirs/run7consolidated/bk11snap141..c4-r%d-n%s_corenone.pkl' % (radrange, noiserange)
            print consolfile

            with open(consolfile, 'rb') as input:

                consol = cPickle.load(input)

                

                plotLogNormDistro(consol['true_m500s'], 
                                  consol['measured_m500s'],
                                  massedges,
                                  meansax,
                                  nongaussax,
                                  stdsax,
                                  offset = offsets[3*i+j],
                                  label = '%s %s' % (radialnames[i], noisenames[j]),
                                  linestyle='None',
                                  linewidth=2.,
                                  color = colors[j],
                                  alpha = alphas[i])

    meansax.set_xscale('log')
    meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    meansax.set_ylabel(r'Mean Bias in $M_{200}$', fontsize=16)
    meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    meansax.set_xlim(2e14, 1.3e15)
    meansax.set_ylim(0.7, 1.2)
    meansax.set_xticks([1e15])
    meansax.set_xticklabels(['10'])
    meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    meansax.legend(patches, labels, loc='upper left')
    meansax.set_title(r'Fit Range: $0.5 < r < 1.5$ Mpc')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('coremxxl_mean_r6.png')

    nongaussax.set_xscale('log')
    nongaussax.set_xlabel('Mass', fontsize=16)
    nongaussax.set_ylabel('Median - Mean', fontsize=16)
    nongaussax.axhline(0.0, c='k', linewidth=2, linestyle='--')
    nongaussax.legend(loc='upper left')
    nongaussfig.canvas.draw()
    nongaussfig.tight_layout()
    nongaussfig.savefig('noisebk11_nongauss.png')


    stdsax.set_xscale('log')
    stdsax.set_xlabel('Mass', fontsize=16)
    stdsax.set_ylabel('Standard Deviation Log-Bias', fontsize=16)
    stdsax.legend()
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('coremxxl_std_r6.png')

    return meansfig, nongaussfig, stdsfig



############################


def plotHSTNoiseNoOffset():



    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)

####
# By Core Radius
##    
#    chaingroups = [['/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2337-5942',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2331-5051',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0533-5005'],
#                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2342-5411',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2106-5844',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0615-5746'],
#                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0000-5748',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2040-5725',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0546-5345'],
#                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0102-4915',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2341-5119'],
#                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2359-5009',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0559-5249']]
#
    

#    clustergroups = [['J2337-5942',
#                   'J2331-5051',
#                   'J0533-5005'],
#                   ['J2342-5411',
#                   'J2106-5844',
#                   'J0615-5746'],
#                   ['J0000-5748',
#                   'J2040-5725',
#                   'J0546-5345'],
#                   ['J0102-4915',
#                   'J2341-5119'],
#                   ['J2359-5009',
#                    'J0559-5249']]
#
#    groupnames = ['0a',
#                  '0b',
#                  '1',
#                  '2',
#                  '3']
#
















#########
# By Redshift
##

    chaingroups = [['/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2331-5051',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0559-5249',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0000-5748'],
                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2359-5009',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2337-5942',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0102-4915'],
                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0533-5005',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2040-5725',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0615-5746'],
                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2341-5119',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0546-5345'],
                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2342-5411',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2106-5844']]


    clustergroups = [['SPT-CLJ2331-5051',
                      'SPT-CLJ0559-5249',
                      'SPT-CLJ0000-5748'],
                   [  'SPT-CLJ2359-5009',
                      'SPT-CLJ2337-5942',
                      'SPT-CLJ0102-4915'],
                   [  'SPT-CLJ0533-5005',
                      'SPT-CLJ2040-5725',
                      'SPT-CLJ0615-5746'],
                   [  'SPT-CLJ2341-5119',
                      'SPT-CLJ0546-5345'],
                   [  'SPT-CLJ2342-5411',
                      'SPT-CLJ2106-5844']]

    groupnames = ['0',
                  '1',
                  '2',
                  '3',
                  '4']


######
# By Number of Bins
#





#
#    chaingroups = [['/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2106-5844',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0546-5345',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0559-5249'],
#                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2342-5411',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0615-5746',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2040-5725'],
#                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0533-5005',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2331-5051',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2337-5942'],
#                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0000-5748',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ0102-4915'],
#                 [  '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2341-5119',
#                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2359-5009']]
#
#
#    clustergroups = [['SPT-CLJ2106-5844',
#                      'SPT-CLJ0546-5345',
#                      'SPT-CLJ0559-5249'],
#                   [  'SPT-CLJ2342-5411',
#                      'SPT-CLJ0615-5746',
#                      'SPT-CLJ2040-5725'],
#                   [  'SPT-CLJ0533-5005',
#                      'SPT-CLJ2331-5051',
#                      'SPT-CLJ2337-5942'],
#                   [  'SPT-CLJ0000-5748',
#                      'SPT-CLJ0102-4915'],
#                   [  'SPT-CLJ2341-5119',
#                      'SPT-CLJ2359-5009']]
#
#    groupnames = ['0',
#                  '1',
#                  '2',
#                  '3',
#                  '4']
#
#





    meansfigs = []
    stdsfigs = []

    for curgroup in range(len(groupnames)):

        chaindirs = chaingroups[curgroup]
        clusternames = clustergroups[curgroup]

        meansfig = pylab.figure()
        meansfigs.append(meansfigs)
        meansax = meansfig.add_subplot(1,1,1)

        stdsfig = pylab.figure()
        stdsfigs.append(stdsfig)
        stdax = stdsfig.add_subplot(1,1,1)



        patches = []
        labels = []


        for i in range(len(clusternames)):

            chaindir = chaindirs[i]

            print chaindir


            label = clusternames[i]

            patch = precomputedLogNormDistro(chaindir, 
                                             massedges,
                                             meansax,
                                             stdax,
                                             colorindex = i%4)

            if patch is None:
                continue

            patches.append(patch)
            labels.append(label)

        meansax.set_title(groupnames[curgroup])
        meansax.set_xscale('log')
        meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
        meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
        meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
        meansax.set_xlim(2e14, 1.3e15)
        meansax.set_ylim(0.65, 1.2)
        meansax.set_xticks([1e15])
        meansax.set_xticklabels(['10'])
        meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
        meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        meansax.legend(patches[::-1], labels[::-1], loc='upper left')
        meansfig.canvas.draw()
        meansfig.tight_layout()
        meansfig.savefig('hstnoisemxxl_logmean_corenone_%s.png' % groupnames[curgroup] )

        stdax.set_title(groupnames[curgroup])
        stdax.set_xscale('log')
        stdax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
        stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
        stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
        stdax.set_xlim(2e14, 1.3e15)
    #    stdax.set_ylim(0.85, 1.10)
        stdax.set_xticks([1e15])
        stdax.set_xticklabels(['10'])
        stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
        stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        stdax.legend(patches[::-1], labels[::-1], loc='upper left')
        stdsfig.canvas.draw()
        stdsfig.tight_layout()
        stdsfig.savefig('hstnoisemxxl_logstd_corenone_%s.png' % groupnames[curgroup])


    return meansfigs, stdsfigs


############################

def plotSplitNoise():



    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    

    chaingroups = [['/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/lineargaussbins-c4-r5-splitdensity-fakeSPT2106',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2106-5844'],
                   ['/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/lineargaussbins-c4-r5-splitdensity-fakeSPT2331',
                    '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-corenone-SPT-CLJ2331-5051']]
    clustergroups = [['Simulated SPT2106 Noise',
                      'Actual SPT2106 Shear Noise'],
                   [  'Simulated SPT2331 Shear Noise',
                      'Actual SPT2331 Shear Noise']]

    groupnames = ['spt2106',
                  'spt2331']





    meansfigs = []
    stdsfigs = []

    for curgroup in range(len(groupnames)):

        print curgroup

        chaindirs = chaingroups[curgroup]
        clusternames = clustergroups[curgroup]

        meansfig = pylab.figure()
        meansfigs.append(meansfigs)
        meansax = meansfig.add_subplot(1,1,1)

        stdsfig = pylab.figure()
        stdsfigs.append(stdsfig)
        stdax = stdsfig.add_subplot(1,1,1)



        patches = []
        labels = []


        for i in range(len(clusternames)):

            chaindir = chaindirs[i]

            print chaindir


            label = clusternames[i]

            patch = precomputedLogNormDistro(chaindir, 
                                             massedges,
                                             meansax,
                                             stdax,
                                             colorindex = i%4)

            if patch is None:
                continue

            patches.append(patch)
            labels.append(label)

        meansax.set_title(groupnames[curgroup])
        meansax.set_xscale('log')
        meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
        meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
        meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
        meansax.set_xlim(2e14, 1.3e15)
        meansax.set_ylim(0.65, 1.2)
        meansax.set_xticks([1e15])
        meansax.set_xticklabels(['10'])
        meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
        meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        meansax.legend(patches[::-1], labels[::-1], loc='upper left')
        meansfig.canvas.draw()
        meansfig.tight_layout()
#        meansfig.savefig('hstnoisemxxl_logmean_corenone_%s.png' % groupnames[curgroup] )

        stdax.set_title(groupnames[curgroup])
        stdax.set_xscale('log')
        stdax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
        stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
        stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
        stdax.set_xlim(2e14, 1.3e15)
    #    stdax.set_ylim(0.85, 1.10)
        stdax.set_xticks([1e15])
        stdax.set_xticklabels(['10'])
        stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
        stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        stdax.legend(patches[::-1], labels[::-1], loc='upper left')
        stdsfig.canvas.draw()
        stdsfig.tight_layout()
 #       stdsfig.savefig('hstnoisemxxl_logstd_corenone_%s.png' % groupnames[curgroup])


    return meansfigs, stdsfigs


        

############################    


def plotHSTNoiseSZOffset():


    delta = 500
    snap = 54
    mcrelation='c4'




    chaingroups = [['/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core0-SPT-CLJ2337-5942',
                    '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core0-SPT-CLJ2331-5051',
                    '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core0-SPT-CLJ0533-5005'],
                   ['/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core0-SPT-CLJ2342-5411',
                    '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core0-SPT-CLJ2106-5844',
                    '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core0-SPT-CLJ0615-5746'],
                   ['/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core1-SPT-CLJ0000-5748',
                    '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core1-SPT-CLJ2040-5725',
                    '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core1-SPT-CLJ0546-5345'],
                   ['/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core2-SPT-CLJ0102-4915',
                    '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core2-SPT-CLJ2341-5119'],
                   ['/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core3-SPT-CLJ2359-5009',
                    '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/hstnoisebins-%s-r5-core3-SPT-CLJ0559-5249']]

    

    clustergroups = [['J2337-5942',
                   'J2331-5051',
                   'J0533-5005'],
                   ['J2342-5411',
                   'J2106-5844',
                   'J0615-5746'],
                   ['J0000-5748',
                   'J2040-5725',
                   'J0546-5345'],
                   ['J0102-4915',
                   'J2341-5119'],
                   ['J2359-5009',
                    'J0559-5249']]

    groupnames = ['core0a',
                  'core0b',
                  'core1',
                  'core2',
                  'core3']






    meansfigs = []
    stdsfigs = []

    for curgroup in range(len(groupnames)):

        chaindirs = chaingroups[curgroup]
        clusternames = clustergroups[curgroup]

        meansfig = pylab.figure()
        meansfigs.append(meansfigs)
        meansax = meansfig.add_subplot(1,1,1)

        stdsfig = pylab.figure()
        stdsfigs.append(stdsfig)
        stdax = stdsfig.add_subplot(1,1,1)



        patches = []
        labels = []


        for i in range(len(clusternames)):

            chaindir = chaindirs[i] % (snap, mcrelation)

            print chaindir


            label = clusternames[i]

            patch = precomputedLogNormDistro(chaindir, 
                                             delta,
                                             meansax,
                                             stdax,
                                             colorindex = i%4)

            if patch is None:
                continue

            patches.append(patch)
            labels.append(label)

        meansax.set_title(groupnames[curgroup])
        meansax.set_xscale('log')
        meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
        meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
        meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
        meansax.set_xlim(2e14, 1.3e15)
        meansax.set_ylim(0.5, 1.05)
        meansax.set_xticks([1e15])
        meansax.set_xticklabels(['10'])
        meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
        meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        meansax.legend(patches[::-1], labels[::-1], loc='upper left')
        meansfig.canvas.draw()
        meansfig.tight_layout()
        meansfig.savefig('hstnoisemxxl_szsimcenter_logmean_%s.snap%d.delta%d.%s.png' % (groupnames[curgroup], snap, delta, mcrelation) )

        stdax.set_title(groupnames[curgroup])
        stdax.set_xscale('log')
        stdax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
        stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
#        stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
        stdax.set_xlim(2e14, 1.3e15)
    #    stdax.set_ylim(0.85, 1.10)
        stdax.set_xticks([1e15])
        stdax.set_xticklabels(['10'])
        stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
        stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        stdax.legend(patches[::-1], labels[::-1], loc='upper left')
        stdsfig.canvas.draw()
        stdsfig.tight_layout()
        stdsfig.savefig('hstnoisemxxl_szsimcenter_logstd_%s.snap%d.delta%d.%s.png' % (groupnames[curgroup], snap, delta, mcrelation) )



    return meansfigs, stdsfigs


    
###############################################


def plotHSTNoiseXrayOffset():

    deltas = [200]

    snaps = [41]

    centerings = 'NONE WTG CCCP SPTHST'.split()
#    centerings = 'xrayNONE corenone'.split()


    
    configtemplate = '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap%(snap)d/hstnoisebins-c4-r5-xray%(centering)s-%(cluster)s'
#    configtemplate = '/users/dapple/euclid1raid1/mxxl_lensing/mxxlsnap%(snap)d/hstnoisebins-c4-r5-%(centering)s-%(cluster)s'


    clusters = ['SPT-CLJ2331-5051',
                'SPT-CLJ0559-5249',
                'SPT-CLJ0000-5748']
#                'SPT-CLJ2359-5009',
#                'SPT-CLJ2337-5942',
#                'SPT-CLJ0102-4915',
#                'SPT-CLJ0533-5005',
#                'SPT-CLJ2040-5725',
#                'SPT-CLJ0615-5746',
#                'SPT-CLJ2341-5119',
#                'SPT-CLJ0546-5345',
#                'SPT-CLJ2342-5411',
#                'SPT-CLJ2106-5844']
#

    meansfigs = []
    stdsfigs = []

    for delta in deltas:

        for snap in snaps:

            simtype = 'mxxlsnap%d' % snap
            massedges = rundln.defineMassEdges(simtype, delta)
#            massedge_sets = [rundln.defineMassEdges(simtype, delta),
#                             np.logspace(np.log10(2e14), np.log10(1e15), 7)]


            for cluster in clusters:
                
                chaindirs = [configtemplate % dict(snap = snap, centering=x, cluster=cluster) for x in centerings]


                meansfig = pylab.figure()
                meansfigs.append(meansfigs)
                meansax = meansfig.add_subplot(1,1,1)

                stdsfig = pylab.figure()
                stdsfigs.append(stdsfig)
                stdax = stdsfig.add_subplot(1,1,1)



                patches = []
                labels = []


                for i in range(len(centerings)):

                    chaindir = chaindirs[i]

                    print chaindir

                    label = centerings[i]

#                    massedges = massedge_sets[i]

                    patch = precomputedLogNormDistro(chaindir, 
                                                     delta,
                                                     massedges,
                                                     meansax,
                                                     stdax,
                                                     colorindex = i%4)

                    if patch is None:
                        continue

                    patches.append(patch)
                    labels.append(label)

                meansax.set_title(cluster)
                meansax.set_xscale('log')
                meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
                meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
                meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
                meansax.set_xlim(2e14, 6e15)
                meansax.set_ylim(0.5, 1.05)
                meansax.set_xticks([1e15])
                meansax.set_xticklabels(['10'])
                meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 
                                    2e15, 3e15, 4e15, 5e15, 6e15], minor=True)
                meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', 
                                         '20', '', '40', '', '60'], minor=True)
                meansax.legend(patches[::-1], labels[::-1], loc='upper left')
                meansfig.canvas.draw()
                meansfig.tight_layout()
                meansfig.savefig('hst_sim_plots/%s-snap%d-delta%d.logmean.png' % (cluster, snap, delta))

                stdax.set_title(cluster)
                stdax.set_xscale('log')
                stdax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
                stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
                stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
                stdax.set_xlim(2e14, 1.6e15)
                #    stdax.set_ylim(0.85, 1.10)
                stdax.set_xticks([1e15])
                stdax.set_xticklabels(['10'])
                stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 
                                    2e15, 3e15, 4e15, 5e15, 6e15], minor=True)
                stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', 
                                         '20', '', '40', '', '60'], minor=True)
                stdax.legend(patches[::-1], labels[::-1], loc='upper left')
                stdsfig.canvas.draw()
                stdsfig.tight_layout()
                stdsfig.savefig('hst_sim_plots/%s-snap%d-delta%d.logstd.png' % (cluster, snap, delta))


    return meansfigs, stdsfigs


    
###############################################


 
def plotNoiseGradient():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)
    meansax.axhline(1.0, c='k', linewidth=1, linestyle='--')

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    


    chaindirs = [ '/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-risingnoise',
                  '/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-flatnoise',
                 '/vol/euclid1/euclid1_raid1/dapple/mxxl_lensing/mxxlsnap41/hstnoisebins-c4-r5-droppingnoise']

    noisenames = ['Dense Core Sampling',
                  'Flat Noise Profile',
                  'Sparse Core Sampling']


    patches = []
    labels = []

    xoffsets = [0.98, 1.0, 1.02]

    for i in range(len(chaindirs)-1, -1,-1):

        chaindir = chaindirs[i]

        print chaindir


        label = noisenames[i]

        patch, summary = precomputedLogNormDistro(chaindir, 
                                                  200,
                                                  meansax,
                                                  stdax,
                                                  colorindex = i, 
                                                  biaslabel=False,
                                                  xoffset = xoffsets[i])

        if patch is None:
            continue

        patches.append(patch)
        labels.append(label)


    meansax.set_xscale('log')
    meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)

    meansax.set_xlim(2e14, 1.3e15)
    meansax.set_ylim(0.70, 1.1)
    meansax.set_xticks([1e15])
    meansax.set_xticklabels(['10'])
    meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    meansax.legend(patches[::-1], labels[::-1], loc='lower left')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('figures/hstnoisemxxl_logmean_noiseprofiles.png')

    stdax.set_xscale('log')
    stdax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
#    stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    stdax.set_xlim(2e14, 1.3e15)
    stdax.set_ylim(0.0, 0.5)
    stdax.set_xticks([1e15])
    stdax.set_xticklabels(['10'])
    stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    stdax.legend(patches[::-1], labels[::-1], loc='lower left')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('figures/hstnoisemxxl_logstd_noiseprofiles.png')


    return meansfig, stdsfig





#############################################


def plotMegacamSnapComp():

    config = 'mega-c4-r6-sigma0.25-corenone'


    clusters = ['SPT-CLJ0307-6225',
                'SPT-CLJ2138-6008',
                'SPT-CLJ0426-5455',
                'SPT-CLJ0317-5935',
                'SPT-CLJ2022-6324',
                'SPT-CLJ0551-5709',
                'SPT-CLJ2145-5644',
                'SPT-CLJ2136-5726',
                'SPT-CLJ2031-5638',
                'SPT-CLJ0346-5438',
                'SPT-CLJ0509-5342',
                'SPT-CLJ0254-5857',
                'SPT-CLJ2332-5358',
                'SPT-CLJ0234-5831',
                'SPT-CLJ0240-5946',
                'SPT-CLJ2032-5627',
                'SPT-CLJ2355-5056',
                'SPT-CLJ0516-5430',
                'SPT-CLJ0348-4514']


#    clusters = ['SPT-CLJ0234-5831',
#               'SPT-CLJ0240-5946',
#               'SPT-CLJ0254-5857',
#               'SPT-CLJ0307-6225',
#               'SPT-CLJ0317-5935',
#               'SPT-CLJ0346-5438',
#               'SPT-CLJ0348-4514',
#               'SPT-CLJ0426-5455',
#               'SPT-CLJ0509-5342',
#               'SPT-CLJ0516-5430',
#               'SPT-CLJ0551-5709',
#               'SPT-CLJ2022-6324',
#               'SPT-CLJ2031-5638',
#               'SPT-CLJ2032-5627',
#               'SPT-CLJ2136-5726',
#               'SPT-CLJ2138-6008',
#               'SPT-CLJ2145-5644',
#               'SPT-CLJ2332-5358',
#               'SPT-CLJ2355-5056']
#
    snapdirs = ['/users/dapple/euclid1raid1/bk11_lensing/snap124/intlength400',
                '/users/dapple/euclid1raid1/bk11_lensing/snap141/intlength400']


    snapnames = ['Snap 124',
                 'Snap 141']

    
    datafile = readtxtfile.readtxtfile('configfiles/megacam_siminput.list')
    redshiftlookup = {}
    densitylookup = {}
    for line in datafile:
        redshiftlookup[line[0]] = float(line[1])
        densitylookup[line[0]] = float(line[2])

    redshifts = np.array([redshiftlookup[x] for x in clusters])
    densities = np.array([densitylookup[x] for x in clusters])
    Dls = np.array([nfwutils.global_cosmology.angulardist(x) for x in redshifts])
    conversion = ((60*180)/(np.pi*Dls))**2
    effdensities = densities*conversion
    

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    for cursnap in range(len(snapdirs)):

        chaindirs = ['%s/%s-%s' % (snapdirs[cursnap], config, x) for x in clusters]

        patches = []
        labels = []

        biasmean = []
        biaserrs = []
        stdmean = []
        stderr = []


        for i in range(len(clusters)):

            chaindir = chaindirs[i]

            print chaindir

            chainfile = '%s/dln_0.chain.0' % (chaindir)
            chain = load_chains.loadChains([chainfile], trim=True)
            print chainfile, len(chain['logmu'])
            if len(chain['logmu'][0,:]) < 5000:
                print 'Skipping'
                continue

            mu, muerr = ci.maxDensityConfidenceRegion(np.exp(chain['logmu'][0,1000::3]))
            biasmean.append(mu)
            biaserrs.append(muerr)
            sig, sigerr = ci.maxDensityConfidenceRegion(np.exp(chain['logsigma'][0,1000::3]))
            stdmean.append(sig)
            stderr.append(sigerr)


        biaserrs = np.array(biaserrs).T
        stderr = np.array(stderr).T

#        meansax.errorbar(redshifts, biasmean, biaserrs, label=snapnames[cursnap], linestyle='none', c=c[cursnap])
#        stdax.errorbar(redshifts, stdmean, stderr, label=snapnames[cursnap], linestyle='none', c=c[cursnap])

        meansax.errorbar(effdensities, biasmean, biaserrs, label=snapnames[cursnap], linestyle='none', c=c[cursnap])
        stdax.errorbar(effdensities, stdmean, stderr, label=snapnames[cursnap], linestyle='none', c=c[cursnap])





#        meansax.set_xscale('log')
#        meansax.set_xlabel('Cluster Redshift', fontsize=16)
        meansax.set_xlabel(r'Galaxy Density [gals/mpc$^2$]', fontsize=16)
        meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
        meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
#        meansax.set_xlim(0.2, 0.7)
        meansax.set_ylim(0.7, 1.1)
#        meansax.set_xticks([1e15])
#        meansax.set_xticklabels(['10'])
#        meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
#        meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        meansax.legend(loc='upper left')
        meansfig.canvas.draw()
        meansfig.tight_layout()
        meansfig.savefig('megacam_snapcomp_effdensity_logmean.png')


#        stdax.set_xscale('log')
#        stdax.set_xlabel('Cluster Redshift', fontsize=16)
        stdax.set_xlabel(r'Galaxy Density [gals/mpc$^2$]', fontsize=16)
        stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
#        stdax.set_xlim(2e14, 1.3e15)
    #    stdax.set_ylim(0.85, 1.10)
#        stdax.set_xticks([1e15])
#        stdax.set_xticklabels(['10'])
#        stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
#        stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        stdax.legend(loc='upper left')
        stdsfig.canvas.draw()
        stdsfig.tight_layout()
        stdsfig.savefig('megacam_snapcomp_effdensity_logstd.png')


    return meansfig, stdsfig



def plotMegacamSnapCompOffset():

    config = 'mega-c4-r6-sigma0.25-core'

    clusters = ['SPT-CLJ0234-5831',
               'SPT-CLJ0240-5946',
               'SPT-CLJ0254-5857',
               'SPT-CLJ0307-6225',
               'SPT-CLJ0317-5935',
               'SPT-CLJ0346-5438',
               'SPT-CLJ0348-4514',
               'SPT-CLJ0426-5455',
               'SPT-CLJ0509-5342',
               'SPT-CLJ0516-5430',
               'SPT-CLJ0551-5709',
               'SPT-CLJ2022-6324',
               'SPT-CLJ2031-5638',
               'SPT-CLJ2032-5627',
               'SPT-CLJ2136-5726',
               'SPT-CLJ2138-6008',
               'SPT-CLJ2145-5644',
               'SPT-CLJ2332-5358',
               'SPT-CLJ2355-5056']

    snapdirs = ['/users/dapple/euclid1raid1/bk11_lensing/snap124/intlength400',
                '/users/dapple/euclid1raid1/bk11_lensing/snap141/intlength400']


    snapnames = ['Snap 124',
                 'Snap 141']

    
    datafile = readtxtfile.readtxtfile('configfiles/megacam_siminput.list')
    redshiftlookup = {}
    corelookup = {}
    for line in datafile:
        redshiftlookup[line[0]] = float(line[1])
        corelookup[line[0]] = int(line[-1])

    redshifts = np.array([redshiftlookup[x] for x in clusters])
    cores = np.array([corelookup[x] for x in clusters])


    

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    for cursnap in range(len(snapdirs)):

        chaindirs = ['%s/%s%d-%s' % (snapdirs[cursnap], config,cores[i], clusters[i]) \
                     for i in range(len(clusters))]

        patches = []
        labels = []

        biasmean = []
        biaserrs = []
        stdmean = []
        stderr = []


        for i in range(len(clusters)):

            chaindir = chaindirs[i]

            print chaindir

            chainfile = '%s/dln_0.chain.0' % (chaindir)
            chain = load_chains.loadChains([chainfile], trim=True)
            print chainfile, len(chain['logmu'])
            if len(chain['logmu'][0,:]) < 5000:
                print 'Skipping'
                continue

            mu, muerr = ci.maxDensityConfidenceRegion(np.exp(chain['logmu'][0,1000::3]))
            biasmean.append(mu)
            biaserrs.append(muerr)
            sig, sigerr = ci.maxDensityConfidenceRegion(np.exp(chain['logsigma'][0,1000::3]))
            stdmean.append(sig)
            stderr.append(sigerr)


        biaserrs = np.array(biaserrs).T
        stderr = np.array(stderr).T

        meansax.errorbar(redshifts, biasmean, biaserrs, label=snapnames[cursnap], linestyle='none', c=c[cursnap])
        stdax.errorbar(redshifts, stdmean, stderr, label=snapnames[cursnap], linestyle='none', c=c[cursnap])





#        meansax.set_xscale('log')
        meansax.set_xlabel('Cluster Redshift', fontsize=16)
        meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
        meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
#        meansax.set_xlim(0.2, 0.7)
        meansax.set_ylim(0.7, 1.1)
#        meansax.set_xticks([1e15])
#        meansax.set_xticklabels(['10'])
#        meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
#        meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        meansax.legend(loc='upper left')
        meansfig.canvas.draw()
        meansfig.tight_layout()
        meansfig.savefig('megacam_snapcomp_core_logmean.png')


#        stdax.set_xscale('log')
        stdax.set_xlabel('Cluster Redshift', fontsize=16)
        stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
#        stdax.set_xlim(2e14, 1.3e15)
    #    stdax.set_ylim(0.85, 1.10)
#        stdax.set_xticks([1e15])
#        stdax.set_xticklabels(['10'])
#        stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
#        stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        stdax.legend(loc='upper left')
        stdsfig.canvas.draw()
        stdsfig.tight_layout()
        stdsfig.savefig('megacam_snapcomp_core_logstd.png')


    return meansfig, stdsfig



def plotMegacamRangeComp():

    rs = [6,9]

    configs = ['mega-c4-r%d-sigma0.25-core' % x for x in rs]

    clusters = ['SPT-CLJ0234-5831',
               'SPT-CLJ0240-5946',
               'SPT-CLJ0254-5857',
               'SPT-CLJ0307-6225',
               'SPT-CLJ0317-5935',
               'SPT-CLJ0346-5438',
               'SPT-CLJ0348-4514',
               'SPT-CLJ0426-5455',
               'SPT-CLJ0509-5342',
               'SPT-CLJ0516-5430',
               'SPT-CLJ0551-5709',
               'SPT-CLJ2022-6324',
               'SPT-CLJ2031-5638',
               'SPT-CLJ2032-5627',
               'SPT-CLJ2136-5726',
               'SPT-CLJ2138-6008',
               'SPT-CLJ2145-5644',
               'SPT-CLJ2332-5358',
               'SPT-CLJ2355-5056']

    snapdir = '/users/dapple/euclid1raid1/bk11_lensing/snap141/intlength400'


    rangenames = ['0.5 - 2.5 Mpc',
                  '0.75 - 2.5 Mpc']

    
    datafile = readtxtfile.readtxtfile('configfiles/megacam_siminput.list')
    redshiftlookup = {}
    corelookup = {}
    for line in datafile:
        redshiftlookup[line[0]] = float(line[1])
        corelookup[line[0]] = int(line[-1])

    redshifts = np.array([redshiftlookup[x] for x in clusters])
    cores = np.array([corelookup[x] for x in clusters])


    

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    for currange in range(len(rs)):

        chaindirs = ['%s/%s%d-%s' % (snapdir, configs[currange], cores[i], clusters[i]) \
                     for i in range(len(clusters))]

        patches = []
        labels = []

        biasmean = []
        biaserrs = []
        stdmean = []
        stderr = []


        for i in range(len(clusters)):

            chaindir = chaindirs[i]

            print chaindir

            chainfile = '%s/dln_0.chain.0' % (chaindir)
            chain = load_chains.loadChains([chainfile], trim=True)
            print chainfile, len(chain['logmu'])
            if len(chain['logmu'][0,:]) < 5000:
                print 'Skipping'
                continue

            mu, muerr = ci.maxDensityConfidenceRegion(np.exp(chain['logmu'][0,1000::3]))
            biasmean.append(mu)
            biaserrs.append(muerr)
            sig, sigerr = ci.maxDensityConfidenceRegion(np.exp(chain['logsigma'][0,1000::3]))
            stdmean.append(sig)
            stderr.append(sigerr)


        biaserrs = np.array(biaserrs).T
        stderr = np.array(stderr).T

        meansax.errorbar(redshifts, biasmean, biaserrs, label=rangenames[currange], linestyle='none', c=c[currange])
        stdax.errorbar(redshifts, stdmean, stderr, label=rangenames[currange], linestyle='none', c=c[currange])





#        meansax.set_xscale('log')
        meansax.set_xlabel('Cluster Redshift', fontsize=16)
        meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
        meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
#        meansax.set_xlim(0.2, 0.7)
        meansax.set_ylim(0.7, 1.1)
#        meansax.set_xticks([1e15])
#        meansax.set_xticklabels(['10'])
#        meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
#        meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        meansax.legend(loc='upper left')
        meansfig.canvas.draw()
        meansfig.tight_layout()
        meansfig.savefig('megacam_rangecomp_core_logmean.png')


#        stdax.set_xscale('log')
        stdax.set_xlabel('Cluster Redshift', fontsize=16)
        stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
#        stdax.set_xlim(2e14, 1.3e15)
    #    stdax.set_ylim(0.85, 1.10)
#        stdax.set_xticks([1e15])
#        stdax.set_xticklabels(['10'])
#        stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
#        stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        stdax.legend(loc='upper left')
        stdsfig.canvas.draw()
        stdsfig.tight_layout()
        stdsfig.savefig('megacam_rangecomp_core_logstd.png')


    return meansfig, stdsfig




########################

def plotMegacamMiscenterCompNoOffset():

    miscenterings = ['core%d', 'szanalytic']
    miscenternames = ['Hydro Sims',
                      'Analytic']

    configs = ['mega-diemer15-r%d-sigma0.25-corenone' % x for x in rs]

#    clusters = ['SPT-CLJ0234-5831',
#               'SPT-CLJ0240-5946',
#               'SPT-CLJ0254-5857',
#               'SPT-CLJ0307-6225',
#               'SPT-CLJ0317-5935',
#               'SPT-CLJ0346-5438',
#               'SPT-CLJ0348-4514',
#               'SPT-CLJ0426-5455',
#               'SPT-CLJ0509-5342',
#               'SPT-CLJ0516-5430',
#               'SPT-CLJ0551-5709',
#               'SPT-CLJ2022-6324',
#               'SPT-CLJ2031-5638',
#               'SPT-CLJ2032-5627',
#               'SPT-CLJ2136-5726',
#               'SPT-CLJ2138-6008',
#               'SPT-CLJ2145-5644',
#               'SPT-CLJ2332-5358',
#               'SPT-CLJ2355-5056']
#
    clusters = ['SPT-CLJ0234-5831',
               'SPT-CLJ0346-5438',
               'SPT-CLJ0426-5455',
               'SPT-CLJ2032-5627',
               'SPT-CLJ2136-5726',
               'SPT-CLJ2138-6008']


    snapdir = '/users/dapple/euclid1_2/rundlns/bk11snap141/'




    datafile = readtxtfile.readtxtfile('configfiles/megacam_siminput.list')
    redshiftlookup = {}
    for line in datafile:
        redshiftlookup[line[0]] = float(line[1])

    redshifts = np.array([redshiftlookup[x] for x in clusters])


    

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    for currange in range(len(rs)):

        chaindirs = ['%s/%s-%s' % (snapdir, configs[currange], clusters[i]) \
                     for i in range(len(clusters))]

        patches = []
        labels = []

        biasmean = []
        biaserrs = []
        stdmean = []
        stderr = []


        for i in range(len(clusters)):

            chaindir = chaindirs[i]

            print chaindir

            chainfile = '%s/dln_0.chain.0' % (chaindir)
            chain = load_chains.loadChains([chainfile], trim=True)
            print chainfile, len(chain['logmu'])
            if len(chain['logmu'][0,:]) < 5000:
                print 'Skipping'
                continue

            mu, muerr = ci.maxDensityConfidenceRegion(np.exp(chain['logmu'][0,1000::3]))
            biasmean.append(mu)
            biaserrs.append(muerr)
            sig, sigerr = ci.maxDensityConfidenceRegion(np.exp(chain['logsigma'][0,1000::3]))
            stdmean.append(sig)
            stderr.append(sigerr)


        biaserrs = np.array(biaserrs).T
        stderr = np.array(stderr).T

        meansax.errorbar(redshifts, biasmean, biaserrs, label=rangenames[currange], linestyle='none', c=c[currange])
        stdax.errorbar(redshifts, stdmean, stderr, label=rangenames[currange], linestyle='none', c=c[currange])





#        meansax.set_xscale('log')
        meansax.set_xlabel('Cluster Redshift', fontsize=16)
        meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
        meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
#        meansax.set_xlim(0.2, 0.7)
        meansax.set_ylim(0.7, 1.1)
#        meansax.set_xticks([1e15])
#        meansax.set_xticklabels(['10'])
#        meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
#        meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        meansax.legend(loc='upper left')
        meansfig.canvas.draw()
        meansfig.tight_layout()
        meansfig.savefig('megacam_rangecomp_nocore_logmean.png')


#        stdax.set_xscale('log')
        stdax.set_xlabel('Cluster Redshift', fontsize=16)
        stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
#        stdax.set_xlim(2e14, 1.3e15)
    #    stdax.set_ylim(0.85, 1.10)
#        stdax.set_xticks([1e15])
#        stdax.set_xticklabels(['10'])
#        stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
#        stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
        stdax.legend(loc='upper left')
        stdsfig.canvas.draw()
        stdsfig.tight_layout()
        stdsfig.savefig('megacam_rangecomp_nocore_logstd.png')


    return meansfig, stdsfig


##################################



def plotMiscenteringCompvSomething(clusters, something, somethinglabel, massbins, prefix, suffix, outdir, ylim=(0.8, 1.1), rs='r5', mc='diemer15', delta=500, snap='bk11snap141', miscentering='xrayNONE'):

    snapnum = int(snap.split('snap')[1])


    offsets = [-0.005,0.005]
    markers = ['o', 's']
    
    snapdir = '/vol/euclid1/euclid1_2/dapple/rundlns/{}'.format(snap)
    config = '{}-{}-{}'.format(prefix, mc, rs)    

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)

    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    for curbini, massbin in enumerate(massbins):

                            
        chaindirs = ['%s/%s-%s-%s-%s' % (snapdir, config, miscentering, 
                                         clusters[i], suffix) \
                     for i in range(len(clusters))]

        patches = []
        labels = []
        
        biasmean = []
        biaserrs = []
        stdmean = []
        stderr = []


        for i in range(len(clusters)):

            chaindir = chaindirs[i]
    
            print chaindir


            chainfile = '%s/rundln%d.%d.%d.chain.0' % (chaindir, snapnum, delta, massbin)
            chain = load_chains.loadChains([chainfile], trim=True)
            print chainfile, len(chain['logmu'])
            if len(chain['logmu'][0,:]) < 5000:
                print 'Skipping'
                continue

            mu, muerr = ci.maxDensityConfidenceRegion(np.exp(chain['logmu'][0,1000::3]))
            biasmean.append(mu)
            biaserrs.append(muerr)
            sig, sigerr = ci.maxDensityConfidenceRegion(np.exp(chain['logsigma'][0,1000::3]))
            stdmean.append(sig)
            stderr.append(sigerr)


        biaserrs = np.array(biaserrs).T
        stderr = np.array(stderr).T

                            

        meansax.errorbar(something + offsets[curbini], biasmean, biaserrs, label='Mass Bin {}'.format(massbin), linestyle='none', c=c[curbini], linewidth=3, marker=markers[curbini])
        stdax.errorbar(something + offsets[curbini], stdmean, stderr, label='Mass Bin {}'.format(massbin), linestyle='none', c=c[curbini], linewidth=3, marker=markers[curbini])

        averr = (biaserrs[0,:] + biaserrs[1,:])/2.
        weights = 1./averr**2
        avebias = np.sum(biasmean*weights)/np.sum(weights)
        meansax.axhline(avebias, c=c[curbini], linewidth=3, linestyle=':')





        
    meansax.set_title('%d %s %s %s %s' % (delta, rs, mc, snap, miscentering))
    meansax.set_xlabel(somethinglabel, fontsize=20)
    meansax.set_ylabel(r'Mean Bias in $Ln(M_{%d})$' % delta, fontsize=20)
    meansax.axhline(1.0, c='k', linewidth=2, linestyle='--')
    meansax.set_ylim(ylim[0], ylim[1])

    meansax.legend(loc='lower right', numpoints=1)
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('%s/%s_%d_%s_%s_%s_%s_%s_logmean.png' % (outdir, prefix, delta, rs, mc, snap, miscentering, somethinglabel))
    meansfig.savefig('%s/%s_%d_%s_%s_%s_%s_%s_logmean.pdf' % (outdir, prefix, delta, rs, mc, snap, miscentering, somethinglabel))
    meansfig.savefig('%s/%s_%d_%s_%s_%s_%s_%s_logmean.eps' % (outdir, prefix, delta, rs, mc, snap, miscentering, somethinglabel))

    stdax.set_title('%d %s %s %s %s' % (delta, rs, mc, snap, miscentering))
    stdax.set_xlabel(somethinglabel, fontsize=16)
    stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
    stdax.legend(loc='upper left')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('%s/%s_%d_%s_%s_%s_%s_%s_logstd.png' % (outdir, prefix, delta, rs, mc, snap, miscentering, somethinglabel))
    stdsfig.savefig('%s/%s_%d_%s_%s_%s_%s_%s_logstd.pdf' % (outdir, prefix, delta, rs, mc, snap, miscentering, somethinglabel))
    stdsfig.savefig('%s/%s_%d_%s_%s_%s_%s_%s_logstd.eps' % (outdir, prefix, delta, rs, mc, snap, miscentering, somethinglabel))



                        

    
    return meansfig, stdsfig

###############################################################

def plotComparison(chaindirs, labels, delta, ylim=(0.8,1.1)):

    alphas = [0.8, 0.3]

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)
        
    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)

    patches = []
    legendlabels = []

    for curi, chaindir, label in zip(range(len(labels)), chaindirs, labels):

        print chaindir

        patch,summary = precomputedLogNormDistro(chaindir, 
                                                 delta,
                                                 meansax,
                                                 stdax,
                                                 colorindex = curi,
                                                 alpha = alphas[curi],
                                                 biaslabel = False)


        if patch is None:
            continue


        patches.append(patch)
        legendlabels.append(label)

    meansax.set_xscale('log')
    meansax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
    meansax.set_ylabel(r'Mean Bias in $Ln(M_{%d})$' % delta, fontsize=16)
    meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    meansax.set_xlim(1e14, 4e15)
    meansax.set_ylim(ylim[0], ylim[1])
    meansax.set_xticks([1e14, 1e15])
    meansax.set_xticklabels(['1', '10'])
    meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
    meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
    meansax.legend(patches[::-1], legendlabels[::-1], loc='upper left')
    meansfig.canvas.draw()
    meansfig.tight_layout()



    stdax.set_xscale('log')
    stdax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
    stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)

    stdax.set_xlim(1e14, 4e15)

    stdax.set_xticks([1e14, 1e15])
    stdax.set_xticklabels(['1', '10'])
    stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
    stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
    stdax.legend(patches[::-1], legendlabels[::-1], loc='upper left')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()


    return meansfig, stdsfig

            

    
###

def plotMiscenteringComp(clusters, prefix, suffix, outdir, ylim=(0.8, 1.1), rs='r5', mc='diemer15'):

    delta = 500

    snaps = ['bk11snap141']


    miscenterings = ['xrayNONE', 'szanalytic', 'szmag']


    names = ['Perfect Centers',
             'Song+12 Miscentering',
             'Hydro Miscentering']


    
    meansfigs = []
    stdsfigs = []

    for curcluster, cluster in enumerate(clusters):

        
        meansfig = pylab.figure()
        meansax = meansfig.add_subplot(1,1,1)
        
        stdsfig = pylab.figure()
        stdax = stdsfig.add_subplot(1,1,1)

        patches = []
        labels = []



        for cursnap, snap in enumerate(snaps):

            snapdir = '/vol/euclid1/euclid1_2/dapple/rundlns/%s' % snap
            config = '%s-%s-%s' % (prefix, mc, rs)

            for curcenter in range(len(miscenterings)):

                chaindir = '%s/%s-%s-%s-%s' % (snapdir, 
                                               config, 
                                               miscenterings[curcenter],
                                               cluster,
                                               suffix) 

                print chaindir

                label = names[curcenter]
                
                patch,summary = precomputedLogNormDistro(chaindir, 
                                                 delta,
                                                 meansax,
                                                 stdax,
                                                 colorindex = curcenter%3,
                                                         alpha=1.0,
                                                 biaslabel = False)


                if patch is None:
                    continue
                    
                if cursnap == 0:
                    patches.append(patch)
                    labels.append(label)
#            


        meansax.set_title(cluster)
        meansax.set_xscale('log')
        meansax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
        meansax.set_ylabel(r'Mean Bias in $Ln(M_{%d})$' % delta, fontsize=16)
        meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
        meansax.set_xlim(1e14, 4e15)
        meansax.set_ylim(ylim[0], ylim[1])
        meansax.set_xticks([1e14, 1e15])
        meansax.set_xticklabels(['1', '10'])
        meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
        meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
        meansax.legend(patches[::-1], labels[::-1], loc='upper left')
        meansfig.canvas.draw()
        meansfig.tight_layout()
        meansfig.savefig('%s/%s_szcentercomp_logmean_%s.delta%d.%s.png' % (outdir, prefix, cluster, delta, mc) )

        stdax.set_title(cluster)
        stdax.set_xscale('log')
        stdax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
        stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
#        stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
        stdax.set_xlim(1e14, 4e15)

        stdax.set_xticks([1e14, 1e15])
        stdax.set_xticklabels(['1', '10'])
        stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
        stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
        stdax.legend(patches[::-1], labels[::-1], loc='upper left')
        stdsfig.canvas.draw()
        stdsfig.tight_layout()
        stdsfig.savefig('%s/%s_szcentercomp_logstd_%s.delta%d.%s.png' % (outdir, prefix, cluster, delta, mc) )


        meansfigs.append(meansfig)
        stdsfigs.append(stdsfig)

    return meansfigs, stdsfigs







    


###############################################################


def dumpConfigList():

    centers = 'xrayNONE xrayXVP szxvptcenter core%d xraylensingpeak szlensingpeak'.split()
    mcs = 'c4 duffy'.split()
    rss = 'r5 r16'.split()

    datafile = asciireader.read('sptdat')
    nametranslator = {}
    for i in range(len(datafile)):
        nametranslator[datafile['name'][i]] = datafile['altname'][i]

    clusters = datafile['name']


    corefileindex = readtxtfile.readtxtfile('shearprofiles/coresizeindex.list')
    corelookup = {}
    for line in corefileindex:
        corelookup[line[0]] = int(line[1])
    cores = np.array([corelookup[nametranslator[x]] for x in clusters])
    config = 'hstnoisebins-{mc}-{rs}-{curcenter}-{clustername}'

    with open('hstsummary_all_configs.list', 'w') as output:

        for rs in rss:
            for center in centers:
                for mc in mcs:
                    for curcluster, clustername in enumerate(clusters):

                        curcenter = center
                        if center == 'core%d':
                            curcenter = center % cores[curcluster]


                        curaltname = nametranslator[clustername]
                        curconfig = config.format(mc = mc, rs = rs, 
                                                  curcenter = curcenter, 
                                                  clustername = curaltname)

                        output.write('{}\n'.format(curconfig))

                



###############################################################


def plotHST_MXXL_BK11_Summary(outputdir, binnum = None):


    
    deltas = [200, 500, 2500]

#    rss = 'r5 r16'.split()
    rss = ['r5']

    mcs = 'c4 diemer15'.split()

    
    centers = 'xrayNONE xraymag szmag szanalytic'.split()



    mxxlsnaps = [41, 54]
    mxxlredshifts = ['z=1.0', 'z=0.25']
    bk11snaps = [124, 141]
    bk11redshifts = ['z=0.5', 'z=0.25']

    config = 'hstnoisebins-{mc}-{rs}-{curcenter}-{clustername}-feb2016'


    datafile = asciireader.read('sptdat')
    nametranslator = {}
    for i in range(len(datafile)):
        nametranslator[datafile['name'][i]] = datafile['altname'][i]

    clusters = datafile['name']


    corefileindex = readtxtfile.readtxtfile('shearprofiles/coresizeindex.list')
    corelookup = {}
    for line in corefileindex:
        corelookup[line[0]] = int(line[1])
    cores = np.array([corelookup[nametranslator[x]] for x in clusters])



    meansfigs = []
    stdsfigs = []

    with open('{}/hstbiassummary'.format(outputdir), 'w') as output:

        output.write('cluster zcluster core sim rad mc delta center b b_err sig sig_err\n')

        for delta in deltas:

            for rs in rss:

                for mc in mcs:

                    print 'MC: ', mc

                    for center in centers:

                        print 'CENTER: ', center

                        prefix = '%s %s %d %s' % (rs, mc, delta, center)
                        output.write('# %s\n' % prefix )

                        makePretty = False


                        for curcluster, clustername in enumerate(clusters):

                            curaltname = nametranslator[clustername]

                            clusterinfo = '%s %1.2f %1.2f' % (curaltname,
                                                              datafile['z_l'][curcluster],
                                                              cores[curcluster])

                            print 'CLUSTER: ', clusterinfo

                            curcenter = center

                            curconfig = config.format(mc = mc, rs = rs, 
                                                      curcenter = curcenter, 
                                                      clustername = curaltname)

                            meansfig = pylab.figure()
                            meansax = meansfig.add_subplot(1,1,1)

                            stdsfig = pylab.figure()
                            stdax = stdsfig.add_subplot(1,1,1)

                            patches = []
                            labels = []


                            #first bk11

                            curcolor = 0

                            for bk11snap, bk11redshift in zip(bk11snaps, bk11redshifts):

                                chaindir = '/users/dapple/euclid1_2/rundlns/bk11snap%d/%s' % (bk11snap, curconfig)

                                try:
                                    patch, summary = precomputedLogNormDistro(chaindir, 
                                                                              delta,
                                                                              meansax,
                                                                              stdax,
                                                                              colorindex = curcolor,
                                                                              biaslabel = False,
                                                                              binnum = binnum)

                                    (avebias, errbias), (avestd, errstd) = summary


                                    output.write('%s BK11%d %s %f %f %f %f\n' % (clusterinfo,
                                                                                 bk11snap, prefix, 
                                                                                 avebias, errbias, 
                                                                                 avestd, errstd))



                                    if patch is None:
                                        print 'Error. Skipped'
                                        continue


                                    patches.append(patch)
                                    labels.append('BK11 %s' % bk11redshift)

                                    makePretty = True


                                    
                                except AssertionError, IOError:
                                    print 'Skipping BK11 %s' % clustername


                                curcolor += 1


                            #then mxxl

                            curcolor = 2
                            
                            for mxxlsnap, mxxlredshift in zip(mxxlsnaps, mxxlredshifts):

                                chaindir = '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/%s' % (mxxlsnap, curconfig)
                                try:
                                    patch, summary = precomputedLogNormDistro(chaindir, 
                                                                              delta,
                                                                              meansax,
                                                                              stdax,
                                                                              colorindex = curcolor,
                                                                              biaslabel = False,
                                                                              binnum = binnum)

                                    (avebias, errbias), (avestd, errstd) = summary



                                    output.write('%s MXXL%d %s %f %f %f %f\n' % (clusterinfo,
                                                                                 mxxlsnap, prefix, 
                                                                                 avebias, errbias, 
                                                                                 avestd, errstd))



                                    if patch is None:
                                        print 'Error. Skipped'
                                        continue


                                    patches.append(patch)
                                    labels.append('MXXL %s' % mxxlredshift)

                                    makePretty = True

                                except AssertionError:
                                    print 'Skipping MXXL %s' % clustername

                                curcolor += 1


                            if makePretty:

                                meansax.set_title('%s' % (clusterinfo))
                                meansax.set_xscale('log')
                                meansax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
                                meansax.set_ylabel(r'Mean Bias in $Ln(M_{%d})$' % delta, fontsize=16)
                                meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
                                meansax.set_xlim(3e14, 4e15)
                                meansax.set_ylim(0.5, 1.3)
                                meansax.set_xticks([1e15])
                                meansax.set_xticklabels(['10'])
                                meansax.set_xticks([3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
                                meansax.set_xticklabels(['', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
                                if mc == 'c4':
                                    meansax.legend(patches[::-1], labels[::-1], loc='upper right')
                                else:
                                    meansax.legend(patches[::-1], labels[::-1], loc='lower left')
                                meansfig.canvas.draw()
                                meansfig.tight_layout()
                                meansfig.savefig('%s/hst_mxxlbk11_comp_logmean_%s.%s.delta%d.%s.%s.png' % (outputdir, clustername, rs, delta, mc, curcenter) )

                                stdax.set_title('%s' % (clusterinfo))
                                stdax.set_xscale('log')
                                stdax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
                                stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
                        #        stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
                                stdax.set_xlim(3e14, 4e15)
                        #        stdax.set_ylim(0.5, 1.05)
                                stdax.set_xticks([1e15])
                                stdax.set_xticklabels(['10'])
                                stdax.set_xticks([3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
                                stdax.set_xticklabels(['', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
                                stdax.legend(patches[::-1], labels[::-1], loc='upper left')
                                stdsfig.canvas.draw()
                                stdsfig.tight_layout()
                                stdsfig.savefig('%s/hst_mxxlbk11_comp_logstd_%s.%s.delta%d.%s.%s.png' % (outputdir, clustername, rs, delta, mc, curcenter) )


                                meansfigs.append(meansfig)
                                stdsfigs.append(stdsfig)




    return meansfigs, stdsfigs





########################


def plotWTG_MXXL_BK11_Summary():

    bk11_snap_ranges = {124 : {500 : 1e14*np.array([1.5, 6.4]),
                               200 : 1e14*np.array([4, 10])},
                        141 : {500 : 1e14*np.array([2, 9]),
                               200 : 1e14*np.array([4, 1.6])}}

    
#    deltas = [200, 500, 2500]
    deltas = [2500, 500]
#    rss = 'r6 r7 r9 r10 r17 r18'.split()
    rss = ['r10']
#    mcs = 'c4 diemer15'.split()
    mcs = ['c4']
    
    centers = 'xrayNONE xrayWTG'.split()


#    noiselevels='n2_4 n3_4 n5_5'.split()
    noiselevels=['n2_4']

    mxxlsnaps = [41, 54]
    mxxlredshifts = ['z=1.0', 'z=0.25']
    bk11snaps = [141, 124]
    bk11redshifts = ['z=0.5', 'z=0.25']



    config = 'general-{mc}-{rs}-{noiselevel}-{curcenter}'

    meansfigs = []
    stdsfigs = []

    with open('wtgbiassummary', 'w') as output:

        output.write('sim rad mc delta center noise b b_err sig sig_err\n')

        for delta in deltas:

            for rs in rss:

                for mc in mcs:

                    print 'MC: ', mc

                    for curcenter in centers:
                        
                        print 'CENTER: ', curcenter

                        for noiselevel in noiselevels:

                            curcolor = -1

                            print 'Noise: ', noiselevel

                            prefix = '%s %s %d %s %s' % (rs, mc, delta, curcenter, noiselevel)
                            output.write('# %s\n' % prefix )

                            makePretty = False


                            curconfig = config.format(mc = mc, rs = rs, 
                                                      noiselevel = noiselevel,
                                                      curcenter = curcenter)

                            meansfig = pylab.figure()
                            meansax = meansfig.add_subplot(1,1,1)

                            stdsfig = pylab.figure()
                            stdax = stdsfig.add_subplot(1,1,1)

                            patches = []
                            labels = []
                            
                            if delta != 2500:

                                #first bk11

                                for bk11snap,bk11redshift in zip(bk11snaps, bk11redshifts):
                                    curcolor += 1

                                    chaindir = '/users/dapple/euclid1_2/rundlns/bk11snap%d/%s' % (bk11snap, curconfig)
                                    chainfile = '%s/rundln%d.%d.0.chain.0' % (chaindir, bk11snap, delta)
                                    try:
                                        chain = load_chains.loadChains([chainfile], trim=True)
                                        print chainfile, len(chain['logmu'])
                                        if len(chain['logmu'][0,:]) < 5000:
                                            print 'Skipping'
                                            continue

                                        mu, muerr = ci.maxDensityConfidenceRegion(np.exp(chain['logmu'][0,1000::3]))
                                        sig, sigerr = ci.maxDensityConfidenceRegion(np.exp(chain['logsigma'][0,1000::3]))



                                        meansax.fill_between(bk11_snap_ranges[bk11snap][delta], 
                                                             mu-muerr[0], mu+muerr[1], 
                                                             facecolor=c[curcolor],alpha=0.8)
                                        stdax.fill_between(bk11_snap_ranges[bk11snap][delta], 
                                                           sig-sigerr[0], sig+sigerr[1], 
                                                           facecolor=c[curcolor], alpha=0.8)

                                        patch = pylab.Rectangle((0, 0), 1, 1, fc=c[curcolor], alpha=0.8, hatch = None)

                                        patches.append(patch)
                                        labels.append('BK11 %s' % bk11redshift)


                                        output.write('BK11%d %s %f %f %f %f\n' % (bk11snap, prefix, 
                                                                                  mu,np.mean(muerr),
                                                                                  sig, np.mean(sigerr)))

                                        makePretty = True


                                    except IOError:
                                        print 'Skipping BK11 %s' % prefix



                            #then mxxl
                            
                            for mxxlsnap, mxxlredshift in zip(mxxlsnaps, mxxlredshifts):
                                curcolor +=1
                            
                                chaindir = '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/%s' % (mxxlsnap, curconfig)
                                try:
                                    patch, summary = precomputedLogNormDistro(chaindir, 
                                                                              delta,
                                                                              meansax,
                                                                              stdax,
                                                                              colorindex = curcolor,
                                                                              biaslabel = False)

                                    (avebias, errbias), (avestd, errstd) = summary



                                    output.write('MXXL%d %s %f %f %f %f\n' % (mxxlsnap, prefix, 
                                                                              avebias, errbias, 
                                                                              avestd, errstd))



                                    if patch is None:
                                        print 'Error. Skipped'
                                        continue


                                    patches.append(patch)
                                    labels.append('MXXL %s' % mxxlredshift)

                                    makePretty = True

                                except AssertionError:
                                    print 'Skipping MXXL %s' % clustername


                            if makePretty:

                                title = '%s %s %d %s %s' % (rs, mc, delta, curcenter, noiselevel.replace('_', '-'))

                                meansax.set_title(title)
                                meansax.set_xscale('log')
                                meansax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
                                meansax.set_ylabel(r'Mean Bias in $Ln(M_{%d})$' % delta, fontsize=16)
                                meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
                                meansax.set_xlim(1e14, 4e15)
                                meansax.set_ylim(0.7, 1.15)
                                meansax.set_xticks([1e14, 1e15])
                                meansax.set_xticklabels(['1', '10'])
                                meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
                                meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
#                                if mc == 'c4':
#                                    meansax.legend(patches[::-1], labels[::-1], loc='upper right')
#                                else:
                                meansax.legend(patches[::-1], labels[::-1], loc='lower right')
                                meansfig.canvas.draw()
                                meansfig.tight_layout()
                                meansfig.savefig('wtg_sim_plots/wtg_mxxlbk11_comp_logmean.%s.delta%d.%s.%s.%s.png' % (rs, delta, mc, curcenter, noiselevel) )

                                stdax.set_title(title)
                                stdax.set_xscale('log')
                                stdax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
                                stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
                        #        stdax.axhline(1.0, c='k', linewidth=3, linestyle='--')
                                stdax.set_xlim(1e14, 4e15)
                        #        stdax.set_ylim(0.5, 1.05)
                                stdax.set_xticks([1e14, 1e15])
                                stdax.set_xticklabels(['1', '10'])
                                stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
                                stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
                                stdax.legend(patches[::-1], labels[::-1], loc='upper left')
                                stdsfig.canvas.draw()
                                stdsfig.tight_layout()
                                stdsfig.savefig('hst_sim_plots/hst_mxxlbk11_comp_logstd.%s.delta%d.%s.%s.%s.png' % (rs, delta, mc, curcenter, noiselevel) )


                                meansfigs.append(meansfig)
                                stdsfigs.append(stdsfig)




    return meansfigs, stdsfigs





########################

def plotNoiseComp():
    '''Compare bias & scatter results for 4 noise levels in the MXXL snap41 simulation. Noise levels are different gals/per arcmin density.'''

    delta = 200    
    rs = 'r10'
    mc = 'c4'
    center = 'xrayNONE'
    
    noiselevels = ['n5_5', 'n3_4', 'n2_4', 'n0_0']
    noiselabels = ['1 gal / arcmin$^2$',
                   '7 gals / arcmin$^2$',
                   '20 gals / arcmin$^2$',
                   'Intrinsic Noise Only']

    xoffsets = [0.97, 0.99, 1.01, 1.03]


                   

    mxxlsnap = 54

    config = 'general-{mc}-{rs}-{noiselevel}-{curcenter}'

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)
    stdsfig = pylab.figure()
    stdax = stdsfig.add_subplot(1,1,1)
                            

    curcolor = 0
    patches = []
    labels = []

    for noiselevel, noiselabel in zip(noiselevels, noiselabels):
        print 'Noise: ', noiselevel

        curconfig = config.format(mc = mc, rs = rs, 
                                  noiselevel = noiselevel,
                                  curcenter = center)
        
        

                            

        


                            
        chaindir = '/vol/euclid1/euclid1_2/dapple/rundlns/mxxlsnap%d/%s' % (mxxlsnap, curconfig)
        try:
            patch, summary = precomputedLogNormDistro(chaindir, 
                                                      delta,
                                                      meansax,
                                                      stdax,
                                                      colorindex = curcolor,
                                                      biaslabel = False,
                                                      xoffset = xoffsets[curcolor])

            (avebias, errbias), (avestd, errstd) = summary

            if patch is None:
                print 'Error. Skipped'
                continue


            patches.append(patch)
            labels.append(noiselabel)

            makePretty = True
            
        except AssertionError:
            print 'Chains empty. Skipped.'

        curcolor += 1



            
    meansax.set_xscale('log')
    meansax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
    meansax.set_ylabel(r'Mean Bias in $Ln(M_{%d})$' % delta, fontsize=16)
    meansax.axhline(1.0, c='k', linewidth=1, linestyle='--')
    meansax.set_xlim(1e14, 4e15)
    meansax.set_ylim(0.7, 1.2)
    meansax.set_xticks([1e14, 1e15])
    meansax.set_xticklabels(['1', '10'])
    meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
    meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
    meansax.legend(patches[::-1], labels[::-1], loc='lower left')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('figures/bias_differingnoiselevels.png')
    meansfig.savefig('figures/bias_differingnoiselevels.eps')
    meansfig.savefig('figures/bias_differingnoiselevels.pdf')


    stdax.set_xscale('log')
    stdax.set_xlabel(r'Mass $M_{%d} [10^{14} M_{\odot}]$' % delta, fontsize=16)
    stdax.set_ylabel(r'Noise Magnitude $\sigma$', fontsize=16)
    stdax.set_xlim(1e14, 4e15)
    stdax.set_xticks([1e14, 1e15])
    stdax.set_xticklabels(['1', '10'])
    stdax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 2e15, 3e15, 4e15], minor=True)
    stdax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '20', '', '40'], minor=True)
    stdax.legend(patches[::-1], labels[::-1], loc='upper left')
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('figures/scatter_differingnoiselevels.png')
    stdsfig.savefig('figures/scatter_differingnoiselevels.eps')
    stdsfig.savefig('figures/scatter_differingnoiselevels.pdf')




    return meansfig, stdsfig





if __name__ == '__main__':

    import matplotlib
    matplotlib.use('agg')

    outputdir = sys.argv[1]

    binnum = None
    if len(sys.argv) >= 3:
        binnum = int(sys.argv[2])
        print 'Bin Number:', binnum

    plotHST_MXXL_BK11_Summary(outputdir, binnum = binnum)
