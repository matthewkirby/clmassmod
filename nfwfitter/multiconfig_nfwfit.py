#!/usr/bin/env python
######################
# This program takes one input and multiple configuration files
#  and runs multiple nfwfits, saving the results in different locations
######################

import sys, os, shutil, glob
import nfwfit
import utils.simutils as simutils
from configparser import ConfigParser


'''
Take the configuration object and kick off a run for each config file

Parameters
----------
configparams : ConfigParser object
               Contains all of the details about the runs the user wishes to use

Returns
-------
None

Notes
-----
MK 2/9: 
    -Unsure what the doTransfer stuff is for. 
    -Commenting this out and seeing if the code works without it.
    -At some point down the line it can be deleted if a purpose is not discovered.
'''

def runMultiConfigs(jobparams, jobname=''):

    outputExt = jobparams['outputExt']
    outbasename = jobparams['outbasename']
    simreader = None

    #inputfiles = jobparams['inputfiles']
    #workbase = jobparams['workbase']
    #doTransfer = workbase is not None
    #if doTransfer:
    #    inputbase = os.path.basename(inputname)
    #
    #    if jobname != '':
    #        jobname = '-' + jobname
    #
    #    workdir = '{0}/{1}{2}'.format(workbase, inputbase, jobname)
    #    print 'WORKDIR: ' + workdir
    #
    #
    #    if not os.path.exists(workdir):
    #        os.mkdir(workdir)
    #
    #    for inputfile in inputfiles:
    #        shutil.copy(inputfile, workdir)
    #
    #
    #    inputname = '{0}/{1}'.format(workdir, inputbase)

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


    
    
'''
Run the NFW fit for each set of configurations

The user is permitted to pass in any amount of configurations. NFWfit will then be run on each
set of configurations.

Parameters
----------
configfile : string
             Name of the config file
             
Returns
-------
None
'''

def runAll(configfile):
    configparams = ConfigParser()
    configparams.read(configfile)
    runMultiConfigs(configparams)


if __name__ == '__main__':

    configfile = sys.argv[1]

    runAll(configfile)

