import batch_nfwfit_jobs as bnj
from nfwfitter.write_multiconfig_jobfile import setupCondor_MXXL
from glob import glob

configs = [ line.strip() for line in open('../nfwfitter/data/run27mxxl41','r') ]
setupCondor_MXXL(configs, '/project/kicp/avestruz/storage/rundirs/mxxlrun27/','mxxlrun27',simdir='/project/kicp/dapple/mxxlsims/mxxlsnap41/',outputdir='/project/kicp/avestruz/storage/mxxlsims/mxxlsnap41/')



bnj.batchNFWFitJobs(glob('/project/kicp/avestruz/storage/rundirs/mxxlrun27/*.job'), '/project/kicp/avestruz/storage/rundirs/mxxlrun27/', nrunners=64, batch_header = bnj.midway_batch_header)
