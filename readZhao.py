#######################
# Read in m-c relations from zhao
#######################

import os
import numpy as np
import scipy.interpolate as interp

#
#  row 1:  ispec,nbin1,zobs,omegam0,omegaL0,sigma8,tilt,omegab0,
#                omeganu0,T_cmb,N_nu,h/note ('---' will take place of
#                those parameters that needn't be input)
#
#        row >1: ziz,Miz,c,rhohiz,R,V,Ms,rhos,Rs,Vs,Miz_200c,c_200c,
#                rhohiz_200c,Miz_200m,c_200m,rhohiz_200m,uniage_iz


if os.environ.get('SRCLOC') is None:
    datdir='./zhaodat'
else:
    datdir='%s/zhaodat' % os.environ.get('SRCLOC')



def readZhao(filename):

    cat = {}

    with open(filename) as input:

        lines = input.readlines()

        ispec,nbin1,zobs,omegam0,omegaL0,sigma8,tilt,omegab0,omeganu0,T_cmb,N_nu,h = lines[0].split()
        cat['zobs'] = float(zobs)
        cat['omegam0'] = float(omegam0)
        cat['omegaL0'] = float(omegaL0)
        cat['sigma8'] = float(sigma8)
        cat['h'] = float(h)

        cat['ziz'] = np.array([float(x.split()[0]) for x in lines[1:]])
        cat['rscale'] = np.array([float(x.split()[8]) for x in lines[1:]])
        cat['M200c'] = np.array([float(x.split()[10]) for x in lines[1:]])
        cat['c200c'] = np.array([float(x.split()[11]) for x in lines[1:]])

    return cat

#####

__MASS_SCALING__ = 1e15


######

class MCFunction(object):

    def __init__(self, interpmc, scaling = __MASS_SCALING__):

        self.interpmc = interpmc
        self.scaling = scaling

    ####

    def __call__(self, m, z):

        return self.interpmc([[z, m/self.scaling]])

######


def createInterp(zi, filenames, scaling = __MASS_SCALING__):

    mcs = [readZhao(afile) for afile in filenames]

    nzs = len(zi)
    nmasses = len(mcs[0]['M200c'])


    xvals = np.zeros((nzs*nmasses,2))
    yvals = np.zeros(nzs*nmasses)

    for i in range(nzs):
        xvals[nmasses*i:nmasses*(i+1),0] = zi[i]
        xvals[nmasses*i:nmasses*(i+1),1] = mcs[i]['M200c']/scaling
        yvals[nmasses*i:nmasses*(i+1)] = mcs[i]['c200c']


    interpmc = MCFunction(interp.LinearNDInterpolator(xvals,yvals), scaling)

    return interpmc

#######
    

def readBCC():

    zs = '00 10 20 25 50 75 100'.split()
    zi = [0.0, 0.1, 0.2, 0.25, 0.5, 0.75, 1.0]
    mcrelation = createInterp(zi, ['%s/bcc_cosmo_mc_z%s.dat' % (datdir,x) for x in zs])

    return mcrelation

#######

def readBK11():

    zs = '25 05'.split()
    zi = [0.25, 0.5]
    mcrelation = createInterp(zi, ['%s/bk11_cosmo_mc_z%s.dat' % (datdir,x) for x in zs])

    return mcrelation
    

####### 

def readMXXL():

    zs = '25 100'.split()
    zi = [0.2425, 0.99]
    mcrelation = createInterp(zi, ['%s/mxxl_cosmo_mc_z%s.dat' % (datdir,x) for x in zs])

    return mcrelation


#######

availableSims = dict(readBCC = readBCC, readBK11 = readBK11, readMXXL = readMXXL)

                            
class ZhaoMC(object):

    def __init__(self, config):

        self.mcrelation = availableSims[config.readermodule]()

    #####

    def __call__(self, m, z, overdensity = None):

        # m in Msol / h

        return self.mcrelation(m, z)
        
