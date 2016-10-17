#!/usr/bin/env python

#########

import rundln, nfwfit, nfwutils
import deconvolvedlognorm as dln
import sys, pymc

#########

__samples__ = 50000
__ngauss__ = 5
__singlecore__ = False


##########

def run(simtype, chaindir, outfile, delta, massbin=0,
        samples = __samples__,
        ngauss = __ngauss__,
        singlecore = __singlecore__):

    outfile = '{}.gauss{}'.format(outfile, ngauss)

    config = nfwfit.readConfiguration('%s/config.sh' % chaindir)
    simreader = nfwfit.buildSimReader(config)
    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    massedges = rundln.defineMassEdges(simtype, delta)
        

    halos = dln.loadPDFs(chaindir, simtype, simreader, delta, massedges, massbin)


    if len(halos) < 10:
        sys.exit(0)

    success = False
    for i in range(20):

        try:
        
            parts = dln.buildGaussMixture1DModel(halos, ngauss, modeltype='ratio')
            success=True
            break

        except (AssertionError, pymc.ZeroProbability) as e:
            
            continue

    assert(success is True)

    with open('%s.massrange' % outfile, 'w') as output:
        output.write('%f\n%f\n' % (massedges[massbin], massedges[massbin+1]))
    dln.sample(parts, outfile, samples, singlecore=singlecore)


if __name__ == '__main__':

    massbin = 0

    simtype=sys.argv[1]
    chaindir=sys.argv[2]
    outfile=sys.argv[3]
    delta=int(sys.argv[4])
    if len(sys.argv) > 5:
        massbin=int(sys.argv[5])

    print 'Called with:', dict(simtype=simtype, chaindir=chaindir, 
                               outfile=outfile, delta=delta, 
                               massbin=massbin)

    
    run(simtype, chaindir, outfile, delta, massbin)
        
