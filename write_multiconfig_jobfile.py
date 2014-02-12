#!/usr/bin/env python
#####################
# This program creates a job file that multiconfig_nfwfit uses to run multiple fits
#  with different configurations
#####################


import sys, os, json, argparse

####################

## add preconfiged functions for slac, condor running

####################

def createJobParams(catalogname, inputfiles, confignames, outputExt, workdir = None, stripCatExt = True):

    

    catalogbase = os.path.basename(catalogname)
    if stripCatExt:
        catalogbase, ext = os.path.splitext(catalogbase)


    jobparams = dir(inputfiles = inputfiles, workdir = workdir, 
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

    parser.add_argument('inputname'
