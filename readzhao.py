#######################
# Read in m-c relations from zhao
#######################

import numpy as np
import scipy.interpolate as interp

#
#  row 1:  ispec,nbin1,zobs,omegam0,omegaL0,sigma8,tilt,omegab0,
#                omeganu0,T_cmb,N_nu,h/note ('---' will take place of
#                those parameters that needn't be input)
#
#        row >1: ziz,Miz,c,rhohiz,R,V,Ms,rhos,Rs,Vs,Miz_200c,c_200c,
#                rhohiz_200c,Miz_200m,c_200m,rhohiz_200m,uniage_iz



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


######

def readBCC():

    zs = '00 10 20 25 50 75 100'.split()
    zi = [0.0, 0.1, 0.2, 0.25, 0.5, 0.75, 1.0]
    mcs = [readzhao.readZhao('zhaodat/bcc_cosmo_mc_z%s.dat' % x) for x in zs]
    xvals = np.zeros((7,47))
    yvals = np.zeros(7*47)

    for i in range(7):
        xvals[47*i:47*(i+1),0] = zi[i]
        xvals[47*i:47*(i+1),1] = mcs[i]['M200c']/1e15
        yvals[47*i:47*(i+1)] = mcs[i]['c200c']


    interpmc = interp.LinearNDInterpolator(xvals,yvals)

    return interpmc


    
    
        


        
