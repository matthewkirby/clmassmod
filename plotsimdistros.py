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


def plotLogNormDistro(truemass, measuredmass, massedges, meanax, stdax, offset = 0., **plotargs):

    log10massedges = np.log10(massedges)

    log10centers = (log10massedges[:-1] + log10massedges[1:])/2.
    nbins = len(log10centers)

    ratio = measuredmass / truemass

    centers = []
    means = []
    stds = []
    

    for i in range(len(log10centers)):

        inbin = ratio[np.logical_and(truemass >= massedges[i],
                               truemass < massedges[i+1])]

        if len(inbin) < 25:
            continue

        centers.append(10**(log10centers[i]))

        logratio = np.log(inbin)

        means.append(bootstrapMean(logratio))
        stds.append(bootstrapStd(logratio))

    centers = np.array(centers)



    meancenter, meanerrs = createErrorbars(means)
    meanax.errorbar(centers + offset, meancenter, meanerrs, **plotargs)

    print meancenter


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
    stdsfig = pylab.figure()
    stdsax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    
    radialranges = [5,9]
    radialnames = ['0.5 - 1.5', '0.75 - 2.5']
    noiseranges = ['0_0', '3_2', '4_3']
    noisenames = ['NoNoise', '20-0.33', '7-0.5']
    offsets = np.linspace(-1.5e13, 1.5e13, 6)


    for i, radrange in enumerate(radialranges):
        for j, noiserange in enumerate(noiseranges):

            consolfile = 'run7consolidated/mxxlsnap41.c4-r%d-n%s_corenone.pkl' % (radrange, noiserange)
            print consolfile

            with open(consolfile, 'rb') as input:

                consol = cPickle.load(input)

                

                plotLogNormDistro(consol['true_m500s'], 
                                  consol['measured_m500s'],
                                  massedges,
                                  meansax,
                                  stdsax,
                                  offset = offsets[3*i+j],
                                  label = '%s %s' % (radialnames[i], noisenames[j]),
                                  linestyle='None',
                                  linewidth=2.)

    meansax.set_xscale('log')
    meansax.set_xlabel('Mass', fontsize=16)
    meansax.set_ylabel('Mean Log-Bias', fontsize=16)
    meansax.legend()
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('noisemxxl_mean.png')

    stdsax.set_xscale('log')
    stdsax.set_xlabel('Mass', fontsize=16)
    stdsax.set_ylabel('Standard Deviation Log-Bias', fontsize=16)
    stdsax.legend()
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('noisemxxl_std.png')

    return meansfig, stdsfig



############################


def plotCoreMXXL():

    meansfig = pylab.figure()
    meansax = meansfig.add_subplot(1,1,1)
    stdsfig = pylab.figure()
    stdsax = stdsfig.add_subplot(1,1,1)

    massedges = np.logspace(np.log10(2e14), np.log10(1e15), 7)
    
    radialranges = [9]
    radialnames = ['0.75 - 2.5']
    coreranges = ['0', '5', '11', 'none']
    offsets = np.linspace(-2.2e13, 2.2e13, 8)


    for i, radrange in enumerate(radialranges):
        for j, corerange in enumerate(coreranges):

            consolfile = 'run7consolidated/mxxlsnap41.c4-r%d-n0_0_core%s.pkl' % (radrange, corerange)
            print consolfile

            with open(consolfile, 'rb') as input:

                consol = cPickle.load(input)

                

                plotLogNormDistro(consol['true_m500s'], 
                                  consol['measured_m500s'],
                                  massedges,
                                  meansax,
                                  stdsax,
                                  offset = offsets[4*i+j],
                                  label = '%s %s' % (radialnames[i], corerange),
                                  linestyle='None',
                                  linewidth=2.)

    meansax.set_xscale('log')
    meansax.set_xlabel('Mass', fontsize=16)
    meansax.set_ylabel('Mean Log-Bias', fontsize=16)
    meansax.legend()
    meansfig.canvas.draw()
    meansfig.tight_layout()
    meansfig.savefig('coremxxl_mean_r9.png')

    stdsax.set_xscale('log')
    stdsax.set_xlabel('Mass', fontsize=16)
    stdsax.set_ylabel('Standard Deviation Log-Bias', fontsize=16)
    stdsax.legend()
    stdsfig.canvas.draw()
    stdsfig.tight_layout()
    stdsfig.savefig('coremxxl_std_r9.png')

    return meansfig, stdsfig

    
    
        

    

    
