#!/usr/bin/env python
#####################
# This program creates a job file that multiconfig_nfwfit uses to run multiple fits
#  with different configurations
#####################


import sys, os, json, argparse, glob

####################

## add preconfiged functions for slac, condor running

####################

def createJobParams(catalogname, confignames, outputExt = '.out', workbase = None, stripCatExt = True):

    inputfiles = glob.glob('{0}*'.format(catalogname))

    catalogbase = os.path.basename(catalogname)
    if stripCatExt:
        catalogbase, ext = os.path.splitext(catalogbase)

    workdir = None


    jobparams = dict(inputfiles = inputfiles, workbase = workbase, 
                    catname = catalogname, outbasename = catalogbase, 
                    configurations = confignames)

    return jobparams


####################



def writeJobfile(jobparams, jobfile):

    with open(jobfile, 'w') as output:

        json.dumps(jobparams, output)

###################


def main(argv = sys.argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('inputname')
