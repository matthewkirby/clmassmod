#/usr/bin/env python
######################
# Utility that burns and thins individual cluster posteriors, and compiles them into one file
######################

import sys, cPickle
import measurebias

answerfile = sys.argv[1]
chaindir = sys.argv[2]
burn = int(sys.argv[3])
thin = int(sys.argv[4])

clusters = measurebias.loadClusterData(answerfile, chaindir, burn, thin)

with open('%s/clusters.pkl' % chaindir, 'wb') as output:
    cPickle.dump(clusters, output)




    
