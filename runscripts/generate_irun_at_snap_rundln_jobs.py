'''Submit the rundlns'''
import batch_nfwfit_jobs as bnj
import sys

if len(sys.argv) < 3 :
    print "Usage: python generate.py irun snapnum"
    print "e.g. python generate.py 30 54"
    sys.exit()

irun = sys.argv[1]
snapnum = sys.argv[2]
    
#configs = [ line.strip() for line in open('../nfwfitter/data/run29mxxl41','r') ]
configs = [ line.strip() for line in open('../nfwfitter/data/run'+irun+'mxxl'+snapnum,'r') ]
#argsets = bnj.buildDLNArgsets('/project/kicp/avestruz/storage/mxxlsims/mxxlsnap41/', configs, 'mxxlsnap41', 200, 'buildMCMCModel_massconcentration','/project/kicp/avestruz/storage/rundlns/')
argsets = bnj.buildDLNArgsets('/project/kicp/avestruz/storage/mxxlsims/mxxlsnap'+snapnum+'/', configs, 'mxxlsnap'+snapnum, 200, 'buildMCMCModel_massconcentration','/project/kicp/avestruz/storage/rundlns/')


# Length of argsets is the number of mass bins * configurations
#bnj.batchRunDLNJobs(argsets, '/project/kicp/avestruz/storage/rundirs/mxxlrun29/',nrunners=64)
bnj.batchRunDLNJobs(argsets, '/project/kicp/avestruz/storage/rundirs/mxxlrun'+irun+'/',nrunners=64)

