#pp!/usr/bin/env python
#####################
# Identify halos that have bad pdfs, as defined by nfwfitter.nfwfit.verifyfit. Rerun them.
######

import sys, glob, os.path
import nfwfitter.nfwfit as nfwfit
import cPickle
#######

def verifyOutput(outfile, simdatadir):
    print outfile
    
    with open(outfile, 'rb') as input:
        fitvals = cPickle.load(input)

    pdfscanner = nfwfit.PDFScanner()
        
    #verify pdfs
    isGood = pdfscanner.verifyfit(None, None, fitvals, None, raiseException = False)
    
    #if any fail, initiate rerun

    if not isGood:

        print 'Reprocessing {}'.format(outfile)
        
        configdir, outputbase = os.path.split(outfile)
        haloid, outext = os.path.splitext(outputbase)

        nfwfit.runNFWFit('{}/{}'.format(simdatadir, haloid),
                         '{}/config.py'.format(configdir),
                         outfile)



def scanConfigs(outputdir, configs, simdatadir):

    for config in configs:

        configdir = '{}/{}'.format(outputdir, config)
        outfiles = glob.glob('{}/*.out'.format(configdir))
        for outfile in outfiles:
            verifyOutput(outfile, simdatadir)
        



if __name__ == '__main__':

    outputdir = sys.argv[1]
    configlist = sys.argv[2]
    simdatadir = sys.argv[3]
    
    configs = [x.strip() for x in open(configlist,'r')]
    print configs
    scanConfigs(outputdir, configs, simdatadir)
