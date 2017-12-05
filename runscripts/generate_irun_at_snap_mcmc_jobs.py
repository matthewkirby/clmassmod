import batch_nfwfit_jobs as bnj
from nfwfitter.write_multiconfig_jobfile import setupCondor_MXXL
from glob import glob

import sys

if len(sys.argv) != 3 :
    print "usage: python [run number] [snapshot]"
    print "e.g. python generate.py 30 41"
    sys.exit()
irun = sys.argv[1]
snapnum = sys.argv[2]

#configs = [ line.strip() for line in open('../nfwfitter/data/run29mxxl41','r') ]
configs = [ line.strip() for line in open('../nfwfitter/data/run'+irun+'mxxl'+snapnum,'r') ]

#setupCondor_MXXL(configs, '/project/kicp/avestruz/storage/rundirs/mxxlrun29/','mxxlrun29',simdir='/project/kicp/dapple/mxxlsims/mxxlsnap41/',outputdir='/project/kicp/avestruz/storage/mxxlsims/mxxlsnap41/')

setupCondor_MXXL(configs, '/project/kicp/avestruz/storage/rundirs/mxxlrun'+irun+'/','mxxlrun'+irun,simdir='/project/kicp/dapple/mxxlsims/mxxlsnap'+snapnum+'/',outputdir='/project/kicp/avestruz/storage/mxxlsims/mxxlsnap'+snapnum+'/')


bnj.batchNFWFitJobs(glob('/project/kicp/avestruz/storage/rundirs/mxxlrun'+irun+'/*.job'), '/project/kicp/avestruz/storage/rundirs/mxxlrun'+irun+'/', nrunners=64, batch_header = bnj.midway_batch_header)
