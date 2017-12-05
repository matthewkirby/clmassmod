'''Submit the rundlns'''

import batch_nfwfit_jobs as bnj
configs = [ line.strip() for line in open('../nfwfitter/data/run27mxxl41','r') ]
argsets = bnj.buildDLNArgsets('/project/kicp/avestruz/storage/mxxlsims/mxxlsnap41/', configs, 'mxxlsnap41', 200, 'buildMCMCModel_massconcentration','/project/kicp/avestruz/storage/rundlns/')

# Length of argsets is the number of mass bins * configurations
bnj.batchRunDLNJobs(argsets, '/project/kicp/avestruz/storage/rundirs/mxxlrun27/',nrunners=64)

