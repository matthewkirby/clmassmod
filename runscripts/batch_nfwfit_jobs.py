#Take nfwfit job files, and batch them up into bash scripts to lessen disk I/O

import os

import nfwfitter.rundln as rundln

aifa_batch_header = '''#!/bin/bash

export HOME=/home/dapple
export PATH=/vol/aibn218/data1/dapple/anaconda/bin:/home/dapple/bin:/vol/software/software/astro/theli/THELI//theli/gui/:/vol/software/software/astro/theli/THELI//theli/bin/Linux_64/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
export LD_LIBRARY_PATH=/home/dapple/lib
export SRCLOC=/vol/euclid1/euclid1_raid1/dapple/mxxlsims
'''

midway_batch_header = '''#!/bin/bash

#SBATCH --job-name={jobname}
#SBATCH --output={jobdir}/{jobname}_{runner}.out
#SBATCH --partition=kicp
#SBATCH --account=kicp
#SBATCH --nodes=1
#SBATCH --time {time}

echo $SLURM_JOB_ID starting execution `date` on `hostname`

'''



def batchNFWFitJobs(jobs, outputdir, nrunners=128, batch_header = aifa_batch_header):
    '''Maximum allowed time on kicp midway is 48:00:00.'''
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    njobs = len(jobs)

    for currunner in range(nrunners):
        curjobs = jobs[currunner::nrunners]
        
        cmds_to_run = []
        for job in curjobs:
            jobbase, jobext = os.path.splitext(job)

            cmds_to_run.append('python nfwfitter/multiconfig_nfwfit.py {jobfile} 1>{jobbase}.stdout 2>{jobbase}.stderr\n'.format(jobfile=job, jobbase=jobbase))

        with open('{}/nfwfitbatch_{}.sh'.format(outputdir, currunner), 'w') as output:
            output.write(batch_header.format(jobdir = outputdir, jobname = 'batchnfwfit', runner=currunner, time='48:00:00'))
            for cmd in cmds_to_run:
                output.write(cmd)

    condorfile = '''executable = /bin/bash
universe = vanilla
Error = {outputdir}/nfwfitbatch_$(Process).stderr
Output = {outputdir}/nfwfitbatch_$(Process).stdout
Log = {outputdir}/nfwfitbatch_$(Process).batch.log
Arguments = {outputdir}/nfwfitbatch_$(Process).sh
queue {nrunners}
'''.format(outputdir = outputdir, nrunners = nrunners)


    with open('{}/nfwfitbatch.submit'.format(outputdir), 'w') as output:
        output.write(condorfile)


#############


def buildDLNArgsets(chainbase, configs, simtype, delta, massmodel, outdir):
    '''chainbase : where the MCMC chains are stored
    configs : list from text file like run25mxxl54
    simtype : e.g. mxxlsnap54
    delta : 200 or 500 in critical
    massmodel : function name of model, e.g. buildMCMCModel_massconcentration
    outdir : where the output goes  
    Note - the length of argsets is the number of mass bins * the number of configurations
    '''
    #number of bins

    massbins = rundln.defineMassEdges(simtype, delta)
    nbins = len(massbins)

    #loop over configs and bins to build arglist

    argsets = []
    for config in configs:
        chaindir = '{}/{}'.format(chainbase, config)
        for massbin in range(nbins):
            workdir = '{outdir}/{config}'.format(outdir = outdir, config = config)
            if not os.path.exists(workdir):
                os.mkdir(workdir)
            outfile = '{workdir}/rundln{simtype}.{massmodel}.{delta}.{massbin}'.format(workdir = workdir,
                                                                                       simtype = simtype,
                                                                                       massmodel = massmodel,
                                                                                       delta = delta,
                                                                                       massbin = massbin)
            argsets.append([simtype, chaindir, outfile, delta, massmodel, massbin])

    return argsets


#############
                               

def batchRunDLNJobs(argsets, outputdir, nrunners=64, prefix='rundlnbatch', batch_header = midway_batch_header):
    '''nrunners should be 64 in Midway (number of simultaneous jobs)'''
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    nsets = len(argsets)
    nperrunner = int(float(nsets)/nrunners)
    timeest = max(60*nperrunner,60)

    for currunner in range(nrunners):
        curargsets = argsets[currunner::nrunners]
        
        cmds_to_run = []
        for argset in curargsets:

            cmds_to_run.append('python nfwfitter/rundln.py {argset}\n'.format(argset = ' '.join(map(str, argset))))

        with open('{}/{}_{}.sh'.format(outputdir, prefix, currunner), 'w') as output:
            output.write(batch_header.format(jobdir = outputdir, jobname = 'batchrundln', runner=currunner, time=timeest))
            for cmd in cmds_to_run:
                output.write(cmd)

    condorfile = '''executable = /bin/bash
universe = vanilla
Error = {outputdir}/{prefix}_$(Process).stderr
Output = {outputdir}/{prefix}_$(Process).stdout
Log = {outputdir}/{prefix}_$(Process).batch.log
Arguments = {outputdir}/{prefix}_$(Process).sh
queue {nrunners}
'''.format(outputdir = outputdir, nrunners = nrunners, prefix=prefix)


    with open('{}/{}.submit'.format(outputdir, prefix), 'w') as output:
        output.write(condorfile)
    
    
############

def checkRunDlnOutput(bashfiles):

    argsets_to_rerun = []

    for bashfile in bashfiles:
        
        with open(bashfile) as input:
            lines = input.readlines()

        for line in lines:
            try:
                tokens = line.split()
                outfile = tokens[4]
                if tokens[2] == 'starting' :
                    continue

                if not os.path.exists(outfile):

                    argsets_to_rerun.append(tokens[2:])

            except IndexError:
                pass

    return argsets_to_rerun





