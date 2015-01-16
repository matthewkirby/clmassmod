import publication_plots as pp
import pylab
import numpy as np
import cPickle

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

c = [(.9,.6,0), (.35, .7, .9), (0,.6,.5)]


def plotLogNormDistro(truemass, measuredmass, massedges, meanax, nongaussax, stdax, label, colorindex):

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

        if (inbin < 0).any():
            print 'ILLEGAL'

        logratio = np.log(inbin)

        medians.append(bootstrapMean(logratio))
        nongausses.append(np.median(logratio) - np.mean(logratio))
        stds.append(bootstrapStd(logratio))

    centers = np.array(centers)



    mediancenter, medianerrs = createErrorbars(medians)

    print mediancenter, medianerrs

    for i in range(len(mediancenter)):
        ylows.append( np.exp(mediancenter[i] - medianerrs[0,i]))
        ylows.append( np.exp(mediancenter[i] - medianerrs[0,i]))
        yhighs.append(np.exp(mediancenter[i] + medianerrs[1,i]))
        yhighs.append(np.exp(mediancenter[i] + medianerrs[1,i]))

#    meanax.errorbar(centers + offset, meancenter, meanerrs, **plotargs)
    print len(xpoints), 
    meanax.fill_between(xpoints, ylows, yhighs, alpha=0.8, color = c[colorindex], label = label, hatch = None)
    patch = pylab.Rectangle((0, 0), 1, 1, fc=c[colorindex], alpha=0.8, hatch = None)
                 


    print nongausses

    nongaussax.plot(centers+offset, nongausses, marker='o', **plotargs)


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
    meansax.legend()
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
    nongaussfig = pylab.figure()
    nongaussax = nongaussfig.add_subplot(1,1,1)
    stdsfig = pylab.figure()
    stdsax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    
    radialranges = [6]
    radialnames = ['0.5 - 1.5', '0.75 - 2.5']
    noiseranges = ['0_0', '4_3']
    noisenames = ['No Noise', '7 gals/sq. arcmin $\sigma_e = 0.5$']
    offsets = np.linspace(-1.5e13, 1.5e13, 6)

    patches = []
    labels = []


    for i, radrange in enumerate(radialranges):
        for j, noiserange in enumerate(noiseranges):

            consolfile = 'mxxl_imperial/rundirs/run7consolidated/mxxlsnap41.c4-r%d-n%s_corenone.pkl' % (radrange, noiserange)
            print consolfile

            with open(consolfile, 'rb') as input:

                consol = cPickle.load(input)

                label = noisenames[j]

                plotLogNormDistro(consol['true_m200s'], 
                                  consol['measured_m200s'],
                                  massedges,
                                  meansax,
                                  nongaussax,
                                  stdsax,
                                  label = '%s %s' % (radialnames[i], noisenames[j]),
                                  colorindex = j)


    meansax.set_xscale('log')
    meansax.set_xlabel(r'Mass $M_{200} [10^{14} M_{\odot}]$', fontsize=16)
    meansax.set_ylabel(r'Mean Bias in $Ln(M_{200})$', fontsize=16)
    meansax.axhline(1.0, c='k', linewidth=3, linestyle='--')
    meansax.set_xlim(2e14, 1.3e15)
    meansax.set_ylim(0.85, 1.10)
    meansax.set_xticks([1e15])
    meansax.set_xticklabels(['10'])
    meansax.set_xticks([2e14, 3e14, 4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 11e14, 12e14, 13e14], minor=True)
    meansax.set_xticklabels(['2', '', '4', '', '6', '', '8', '', '', '12', ''], minor=True)
    meansax.legend(patches, labels, loc='upper left')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('noisemxxl_median.png')

    nongaussax.set_xscale('log')
    nongaussax.set_xlabel('Mass', fontsize=16)
    nongaussax.set_ylabel('Median - Mean', fontsize=16)
    nongaussax.axhline(0.0, c='k', linewidth=2, linestyle='--')
    nongaussax.legend(loc='upper left')
    nongaussfig.canvas.draw()
    nongaussfig.tight_layout()
    nongaussfig.savefig('noisemxxl_nongauss.png')


    stdsax.set_xscale('log')
    stdsax.set_xlabel('Mass', fontsize=16)
    stdsax.set_ylabel('Standard Deviation Log-Bias', fontsize=16)
    stdsax.legend()
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('noisemxxl_std.png')

    return meansfig, nongaussfig, stdsfig



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
        

    

    
