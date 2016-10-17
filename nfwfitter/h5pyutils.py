###############################
#From Joerg Dietrich 2013
#Modified by Douglas Applegate 2013



import h5py
import numpy as np

def getClusterProperties(f):
    ra = f['Cluster']['Cluster']['ra']
    dec = f['Cluster']['Cluster']['dec']
    z = f['Cluster']['Cluster']['z']
    mass = f['Cluster']['Cluster']['m200']
    return ra, dec, z, mass

def getShearCat(f):
    data = f['Shear']['Shear Catalog'][...]
    data = data.view(np.recarray)
    return data
