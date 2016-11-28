#Take nfwfit job files, and batch them up into bash scripts to lessen disk I/O

import os

aifa_batch_header = '''#!/bin/bash

export HOME=/home/dapple
export PATH=/vol/aibn218/data1/dapple/anaconda/bin:/home/dapple/bin:/vol/software/software/astro/theli/THELI//theli/gui/:/vol/software/software/astro/theli/THELI//theli/bin/Linux_64/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
export LD_LIBRARY_PATH=/home/dapple/lib
export SRCLOC=/vol/euclid1/euclid1_raid1/dapple/mxxlsims
'''

midway_batch_header = '''#!/bin/bash

export HOME=/home/dapple
export PATH=/home/dapple/anaconda2/bin:/software/slurm-current-el6-x86_64/bin:/software/gnuplot-4.6-el6-x86_64/bin:/software/xpdf-3.04-el6-x86_64/bin:/software/tkdiff-4.2
export LD_LIBRARY_PATH=/home/dapple/lib
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

            cmds_to_run.append('python /vol/euclid1/euclid1_raid1/dapple/mxxlsims/multiconfig_nfwfit.py {jobfile} 1>{jobbase}.stdout 2>{jobbase}.stderr\n'.format(jobfile=job, jobbase=jobbase))

        with open('{}/nfwfitbatch_{}.sh'.format(outputdir, currunner), 'w') as output:
            output.write(batch_header)
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


#############
                               

def batchRunDLNJobs(argsets, outputdir, nrunners=128, prefix='rundlnbatch'):

    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    nsets = len(argsets)
    nperrunner = int(float(nsets)/nrunners)

    for currunner in range(nrunners):
        curargsets = argsets[currunner::nrunners]
        
        cmds_to_run = []
        for argset in curargsets:

            cmds_to_run.append('python /vol/euclid1/euclid1_raid1/dapple/mxxlsims/rundln.py {argset}\n'.format(argset = ' '.join(argset)))

        with open('{}/{}_{}.sh'.format(outputdir, prefix, currunner), 'w') as output:
            output.write(batch_header)
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





