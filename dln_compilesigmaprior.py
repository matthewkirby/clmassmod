#############
''' Use a no-noise chain to measure a prior on logsigma for low noise chains.'''

import sys, glob, os.path
import cPickle
import numpy as np
import load_chains
import pylab

#############


def run(inputchaindir, outputdir):

    massbinchainfiles = glob.glob('{}/*.chain.0'.format(inputchaindir))

    sigmapriors = {}

    for chainfile in massbinchainfiles:

        chainbase = os.path.basename(chainfile)

        massbin = int(chainbase.split('.')[2])
        
        chain = load_chains.loadChains([chainfile])
        
        logsigmas = chain['logsigma'][0,1000::3]
        mu = np.mean(logsigmas)
        std = np.std(logsigmas)
        xs = np.linspace(np.min(logsigmas), np.max(logsigmas), 1000)
        model = np.exp(-0.5*((xs - mu)/std)**2)/(std*np.sqrt(2*np.pi))

        sigmapriors[massbin] = (mu, std)
        
        fig = pylab.figure()
        ax = pylab.gca()
        ax.hist(logsigmas, bins=50, normed=True)
        ax.plot(xs, model, 'r-')
        ax.set_title(chainfile)
        fig.savefig('{}/{}.logsigma.png'.format(outputdir, chainbase))

    cPickle.dump(sigmapriors, open('{}/sigmapriors.pkl'.format(outputdir), 'wb'), -1)


############

if __name__ == '__main__':

    inputchaindir = sys.argv[1]
    outputdir = sys.argv[2]

    run(inputchaindir, outputdir)
