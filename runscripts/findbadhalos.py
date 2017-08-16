#!/usr/bin/env python
#####################
# Identify halos that have bad pdfs, as defined by nfwfitter.nfwfit.verifyfit. Rerun them.
######

import sys, glob, os.path
import nfwfitter.nfwfit as nfwfit
import readtxtfile

#######

def verifyOutput(outfile):

    with open(pdffile, 'rb') as input:
        fitvals = cPickle.load(input)

    #verify pdfs
    isGood = nfwfit.verifyfit(None, None, fitvals, None, raiseException = False)
    
    #if any fail, initiate rerun

    if not isGood:

        print 'Reprocessing {}'.format(outfile)

        configdir, outputbase = os.path.split(outfile)
        haloid, outext = os.path.splitext(outputbase)

        nfwfit.runNFWFit('{}/../{}'.format(configdir, haloid),
                         '{}/config.py'.format(configdir),
                         outfile)



def scanConfigs(simdir, configs):

    for config in configs:

        configdir = '{}/{}'.format(simdir, config)
        outfiles = glob.glob('{}/*.out')
        for outfile in outfiles:
            verifyOutput(outfile)
        



if __name__ == '__main__':

    simdir = sys.argv[1]
    configlist = sys.argv[2]

    configs = [x[0] for x in readtxtfile.readtxtfile(configlist)]

    scanConfigs(simdir, configs)
