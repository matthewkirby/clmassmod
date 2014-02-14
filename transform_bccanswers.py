#!/usr/bin/env python
#######################
# This file compiles the truth tables for BCC clusters
#######################

import h5py, h5pyutils, sys, os, glob, cPickle
import readtxtfile

#######################



clusterlist = sys.argv[1]
outfile = sys.argv[2]

clusterfiles = [x[0] for x in readtxtfile.readtxtfile(clusterlist)]

clusterstofind = {}

clusterinfo = {}

for i, clusterfile in enumerate(clusterfiles):

    root, ext = os.path.splitext(clusterfile)
    clusterid = float(root.split('_')[1])
    clusterstofind[clusterid] = None


for halofile in glob.glob('/users/dapple/localdisk/halos/*.fit'):
    
    halocat = ldac.LDACCat(pyfits.open(halofile)[1])

    for halo_index, id in enumerate(halocat['HALOID']):

        if id in clusterstofind:

            clusterinfo[id] = dict(m500 = cat['M500'][halo_index], 
                                   m200 = cat['M200C'][halo_index],
                                   concen = cat['R200'][halo_index]/cat['RS'][halo_index],
                                   redshift = cat['Z'][halo_index])

            del clusterstofind[id]


assert(len(clusterstofind == 0))


with open(outfile, 'wb') as output:

    cPickle.dump(clusterinfo, output, -1)





    

    


    
