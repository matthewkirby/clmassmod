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


def plotLogNormDistro(truemass, measuredmass, massedges, meanax, nongaussax, stdax, offset = 0., **plotargs):

    log10massedges = np.log10(massedges)

    log10centers = (log10massedges[:-1] + log10massedges[1:])/2.
    nbins = len(log10centers)

    ratio = measuredmass / truemass

    centers = []
    means = []
    nongausses = []
    stds = []
    

    for i in range(len(log10centers)):

        inbin = ratio[np.logical_and(truemass >= massedges[i],
                               truemass < massedges[i+1])]

        if len(inbin) < 25:
            continue

        centers.append(10**(log10centers[i]))

        if (inbin < 0).any():
            print 'ILLEGAL'

        logratio = np.log(inbin)

        means.append(bootstrapMean(logratio))
        nongausses.append(np.median(logratio) - np.mean(logratio))
        stds.append(bootstrapStd(logratio))

    centers = np.array(centers)



    meancenter, meanerrs = createErrorbars(means)
    meanax.errorbar(centers + offset, meancenter, meanerrs, **plotargs)

    print nongausses

    nongaussax.plot(centers+offset, nongausses, marker='o', **plotargs)


    stdcenter, stderrs = createErrorbars(stds)
    stdax.errorbar(centers + offset, stdcenter, stderrs, **plotargs)


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
    
    radialranges = [8]
    radialnames = ['0.5 - 1.5', '0.75 - 2.5']
    noiseranges = ['0_0', '3_2', '4_3']
    noisenames = ['NoNoise', '20-0.33', '7-0.5']
    offsets = np.linspace(-1.5e13, 1.5e13, 3)

    colors = 'b g r'.split()
    alphas = [1.0, 0.5]


    for i, radrange in enumerate(radialranges):
        for j, noiserange in enumerate(noiseranges):

            consolfile = '../rundirs/run7consolidated/mxxlsnap41.c4-r%d-n%s_corenone.pkl' % (radrange, noiserange)
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
    meansax.set_xlabel('Mass', fontsize=16)
    meansax.set_ylabel('Mean Log-Bias', fontsize=16)
    meansax.legend()
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('noisemxxl_mean.png')

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

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    
    radialranges = [8, 5]
    radialnames = ['0.75 - 1.5', '0.5-1.5']
    coreranges = ['none', '0', '5']
    corenames = ['none', '0.25', '1.5']
    offsets = np.linspace(-2.2e13, 2.2e13, 6)

    colors = ['b', 'g', 'r']
    alphas = [1.0, 0.5]


    for i, radrange in enumerate(radialranges):
        for j, corerange in enumerate(coreranges):

            consolfile = '../rundirs/run7consolidated/mxxlsnap41.c4-r%d-n4_3_core%s.pkl' % (radrange, corerange)
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
    meansax.legend(loc='lower left')
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('coremxxl_mean.png')

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
    meansax.set_xlabel('Mass', fontsize=16)
    meansax.set_ylabel('Mean Log-Bias', fontsize=16)
    meansax.legend()
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('noisebk11_mean.png')

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
    stdsfig.savefig('noisebk11_std.png')

    return meansfig, nongaussfig, stdsfig



############################
        

    

    
