#!/usr/bin/env python
######################
# This program takes one input and multiple configuration files
#  and runs multiple nfwfits, saving the results in different locations
######################

import sys, json, os, shutil, glob
import nfwfit


def runMultiConfigs(jobparams):

    inputfiles = jobparams['inputfiles']
    workdir = jobparams['workdir']
    doTransfer = workdir is not None

    try:

        inputname = jobparams['catname']
        outbasename = jobparams['outbasename']

        if doTransfer:
            if not os.path.exists(workdir):
                os.makedirs(workdir)

            for inputfile in inputfiles:
                shutil.copy(inputfile, workdir)

            inputbase = os.path.basename(inputname)
            inputname = '{0}/{1}'.format(workdir, catname)



        for configfile in jobparams['configurations']:

            outdir = os.path.dirname(configfile)

            outputname = '{0}/{1}'.format(outdir, outbasename)

            nfwfit.runNFWFit(catname, configfile, outputname)

    finally:

        if doTransfer:
            workfiles = glob.glob('{0}/*'.format(workdir))
            for workfile in workfiles:
                os.remove(workfile)

########

def loadJobfile(jobfile):

    with open(jobfile) as input:
        jobparams = json.loads(input)

    return jobparams
    
#########

def runAll(jobfile):

    jobparams = loadJobfile(jobfile)
    runMultiConfigs(jobparams)

#########

if __name__ == '__main__':

    jobfile = sys.argv[1]

    runAll(jobfile)

