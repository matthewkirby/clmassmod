#!/usr/bin/env python
#####################
# This program creates a job file that multiconfig_nfwfit uses to run multiple fits
#  with different configurations
#####################


import sys, os, json, argparse, glob

####################

## add preconfiged functions for slac, condor running

####################

def setupCondor_MXXL(configs, jobdir, jobname, simdir = '/home/dapple/braid1/mxxl/snap41', toInclude = range(6300)):
    
    if not os.path.exists(jobdir):
        os.mkdir(jobdir)

    configfiles = ['{0}/{1}/config.sh'.format(simdir, config) for config in configs]

    input_extensions = 'convergence_map shear_1_map shear_2_map'.split()

    for i, id in enumerate(toInclude):

        catname = '{0}/halo_cid{1}'.format(simdir, id)

        inputfiles = ['{0}.{1}'.format(catname, x) for x in input_extensions]
        
        jobparams = createJobParams(catname,
                                    configfiles,
                                    inputfiles = inputfiles,
                                    workbase = './',
                                    stripCatExt = False)
        writeJobfile(jobparams, '{0}/{1}.{2}.job'.format(jobdir, jobname, i))

    condorfile = '''executable = /home/dapple/braid1/mxxl/mxxlsims/nfwfit_condorwrapper.sh
universe = vanilla
Error = {jobdir}/{jobname}.$(Process).stderr
Output = {jobdir}/{jobname}.$(Process).stdout
Log = {jobdir}/{jobname}.$(Process).batch.log
Arguments = {jobdir}/{jobname}.$(Process).job
queue {njobs}
'''.format(jobdir = jobdir, jobname = jobname, njobs = len(toInclude))

    with open('{0}/{1}.submit'.format(jobdir, jobname), 'w') as output:
        output.write(condorfile)

        

    
    


####################

def createJobParams(catalogname, confignames, inputfiles, outputExt = '.out', workbase = None, stripCatExt = True):


    catalogbase = os.path.basename(catalogname)
    if stripCatExt:
        catalogbase, ext = os.path.splitext(catalogbase)

    workdir = None


    jobparams = dict(inputfiles = inputfiles, outputExt = outputExt, workbase = workbase, 
                    catname = catalogname, outbasename = catalogbase, 
                    configurations = confignames)

    return jobparams


####################



def writeJobfile(jobparams, jobfile):

    with open(jobfile, 'w') as output:

        json.dump(jobparams, output)

###################


def main(argv = sys.argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('inputname')
