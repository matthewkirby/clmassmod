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

    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    njobs = len(jobs)
    nperrunner = int(float(njobs)/nrunners)

    for currunner in range(nrunners):
        curjobs = jobs[currunner::nrunners]
        
        cmds_to_run = []
        for job in curjobs:
            jobbase, jobext = os.path.splitext(job)

            cmds_to_run.append('python nfwfitter/multiconfig_nfwfit.py {jobfile} 1>{jobbase}.stdout 2>{jobbase}.stderr\n'.format(jobfile=job, jobbase=jobbase))

        with open('{}/nfwfitbatch_{}.sh'.format(outputdir, currunner), 'w') as output:
            output.write(batch_header.format(jobdir = outputdir, jobname = 'batchnfwfit', runner=currunner, time='02:00:00'))
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


def buildDLNArgsets(chainbase, configs, simtype, delta, outdir):

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
            outfile = '{workdir}/rundln{simtype}.{delta}.{massbin}'.format(workdir = workdir,
                                                                           simtype = simtype,
                                                                           delta = delta,
                                                                           massbin = massbin)
            argsets.append([simtype, chaindir, outfile, delta, 'pdf', massbin])

    return argsets


#############
                               

def batchRunDLNJobs(argsets, outputdir, nrunners=128, prefix='rundlnbatch', batch_header = midway_batch_header):

    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    nsets = len(argsets)
    nperrunner = int(float(nsets)/nrunners)
    timeest = 10*nperrunner

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

                if not os.path.exists(outfile):

                    argsets_to_rerun.append(tokens[2:])

            except IndexError:
                pass

    return argsets_to_rerun





