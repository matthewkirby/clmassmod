#!/usr/bin/env python
######################
# This program takes one input and multiple configuration files
#  and runs multiple nfwfits, saving the results in different locations
######################

import sys, json, os, shutil, glob
import nfwfit
import utils.simutils as simutils


def runMultiConfigs(jobparams, jobname=''):

    inputfiles = jobparams['inputfiles']
    outputExt = jobparams['outputExt']
    workbase = jobparams['workbase']
    doTransfer = workbase is not None



    inputname = jobparams['catname']
    outbasename = jobparams['outbasename']

    if doTransfer:

        inputbase = os.path.basename(inputname)

        if jobname != '':
            jobname = '-' + jobname

        workdir = '{0}/{1}{2}'.format(workbase, inputbase, jobname)
        print 'WORKDIR: ' + workdir


        if not os.path.exists(workdir):
            os.mkdir(workdir)

        for inputfile in inputfiles:
            shutil.copy(inputfile, workdir)


        inputname = '{0}/{1}'.format(workdir, inputbase)

        

    simreader = None



    for configfile in jobparams['configurations']:

        outdir = os.path.dirname(configfile)

        outputname = '{0}/{1}{2}'.format(outdir, outbasename, outputExt)

        print configfile, outputname

        if simreader is None:

            config, simreader = nfwfit.preloadNFWFit(configfile)

        else:

            config = simutils.readConfiguration(configfile)

        nfwfit.runNFWFit_Preloaded(simreader, inputname, config, outputname)





    if doTransfer:
        
        shutil.rmtree(workdir)

########

def loadJobfile(jobfile):

    with open(jobfile) as input:
        jobparams = json.load(input)

    return jobparams
    
#########

def runAll(jobfile):

    jobparams = loadJobfile(jobfile)
    runMultiConfigs(jobparams)

#########

if __name__ == '__main__':

    jobfile = sys.argv[1]

    runAll(jobfile)

