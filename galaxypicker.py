'''
Select Galaxies from a simulation for analysis
'''

############

import simutils
import numpy as np

###########

class Composite(object):

    def __init__(self, *pickers):
        self.pickers = pickers

    def __call__(self, sim):

        curcat = sim
        for picker in self.pickers:
            curcat = picker(curcat)

        return curcat

##################

class GalaxyPicker(object):

    def __init__(self, prevpicker, config):

        self.config = config

    def __call__(self, sim):

        curmask = self.mask(sim)

        return sim.filter(curmask)

###################

class AllGalaxyPicker(object):

    def __call__(self, sim):

        return sim

###################

class InsufficientGalaxiesException(Exception): pass

class DensityPicker(GalaxyPicker):
    '''assumes that the input catalog is rectalinear'''

    def mask(self, sim):

        x_arcmin = sim.x_arcmin
        y_arcmin = sim.y_arcmin
        zcluster = sim.zlens

        targetdensity = self.config.nperarcmin        

        max_x = np.max(x_arcmin)
        min_x = np.min(x_arcmin)
        delta_x = max_x - min_x

        max_y = np.max(y_arcmin)
        min_y = np.min(y_arcmin)
        delta_y = max_y - min_y

        area = delta_x*delta_y
        
        if targetdensity == -1:
        #take all
            targetnumber = len(x_arcmin)

        elif 'targetz' in self.config:
        #adjust for different redshift
            curr_angdist = nfwutils.global_cosmology.angulardist(zcluster)
            newangdist = nfwutils.global_cosmology.angulardist(self.config.targetz)
            ratio = curr_angdist/newangdist
            newarea = area*ratio**2
            targetnumber = targetdensity*newarea
        else:

            targetnumber = targetdensity*area

        availablenumber = len(x_arcmin)

        if targetnumber > availablenumber:
            raise InsufficientGalaxiesException

        accept = float(targetnumber) / availablenumber

        randomthrow = np.random.random(len(x_arcmin))

        selected = randomthrow < accept

        return selected

################################################
#####################

def squaremask(x_arcmin, y_arcmin, x=0,y=0,theta=0, sidelength = 1.):
    #x,y in arcminutes

    theta_rad = np.pi*theta/180.

    dX = x_arcmin - x   
    dY = y_arcmin - y

    rdX = np.cos(theta_rad)*dX - np.sin(theta_rad)*dY
    rdY = np.sin(theta_rad)*dX + np.cos(theta_rad)*dY

    return (np.abs(rdX) + np.abs(rdY)) <= (np.sqrt(2)*sidelength/2.)

def rectanglemask(x_arcmin, y_arcmin, x=0,y=0,theta=0, xlength=1., ylength=1.):

    #x,y in arcminutes

    theta_rad = np.pi*theta/180.

    dX = x_arcmin - x   
    dY = y_arcmin - y

    rdX = np.cos(theta_rad)*dX - np.sin(theta_rad)*dY
    rdY = np.sin(theta_rad)*dX + np.cos(theta_rad)*dY

    return np.logical_and(np.abs(rdX) < xlength/2., np.abs(rdY) < ylength/2.)


def circlemask(x_arcmin, y_arcmin, x=0, y=0, rad=1.):
    #arcmnutes

    dX = x_arcmin - x
    dY = y_arcmin - y

    return np.sqrt(dX**2 + dY**2) < rad

acsmask = lambda x_arcmin,y_arcmin,x,y: squaremask(x_arcmin,y_arcmin,x,y,sidelength=3.2)
wfc3mask = lambda x_arcmin,y_arcmin,x,y: squaremask(x_arcmin,y_arcmin,x,y,theta=45.,
                                                    sidelength=(3.2*4./5.))

acscentered = lambda x_arcmin,y_arcmin: np.logical_or(acsmask(x_arcmin,y_arcmin,0,0), 
                                                      wfc3mask(x_arcmin,y_arcmin,0, 6.))

offsetpointing = lambda x_arcmin,y_arcmin: np.logical_or(acsmask(x_arcmin,y_arcmin,0.,-3.), 
                                                         wfc3mask(x_arcmin,y_arcmin,0., 3.0))

rotatedoffset = lambda x_arcmin,y_arcmin: np.logical_or(acsmask(x_arcmin,y_arcmin,-3.,0.), 
                                                        wfc3mask(x_arcmin,y_arcmin,3., 0.0))

offsetmosaic = lambda x_arcmin,y_arcmin: np.logical_or(offsetpointing(x_arcmin,y_arcmin), 
                                                       rotatedoffset(x_arcmin,y_arcmin))

offset3 = lambda x_arcmin,y_arcmin: np.logical_or(offsetmosaic(x_arcmin,y_arcmin), 
                                                  acscentered(x_arcmin,y_arcmin))

pisco3 = lambda x_arcmin,y_arcmin: np.logical_or(np.logical_or(rectanglemask(x_arcmin,y_arcmin,0,0,0,8,6), 
                                                               rectanglemask(x_arcmin,y_arcmin,0,7,0,6,8)), 
                                                 rectanglemask(x_arcmin,y_arcmin,0,-7,0,6,8))

pisco4 = lambda x_arcmin,y_arcmin: np.logical_or(np.logical_or(rectanglemask(x_arcmin,y_arcmin,-5.5, 0, 0, 8, 6), 
                                                               rectanglemask(x_arcmin,y_arcmin,5.5, 0, 0, 8, 6)), 
                               np.logical_or(rectanglemask(x_arcmin,y_arcmin,0, 5.5, 0, 6, 8), 
                                             rectanglemask(x_arcmin,y_arcmin,0, -5.5, 0,6, 8)))

def randomoffset(x_arcmin,y_arcmin):
    posangle = np.random.uniform(0, 2*np.pi)
    rotmatrix = np.array([[np.cos(posangle), -np.sin(posangle)],
                          [np.sin(posangle), np.cos(posangle)]])
    acspos = np.dot(rotmatrix, np.array([0., -3]))
    wfc3pos = np.dot(rotmatrix, np.array([0., 3]))

    print acspos, wfc3pos

    rotatedoffset = np.logical_or(acsmask(x_arcmin,y_arcmin,acspos[0], acspos[1]), 
                                  wfc3mask(x_arcmin,y_arcmin,wfc3pos[0], wfc3pos[1]))
    return rotatedoffset

randomoffsetmosaic = lambda x_arcmin,y_arcmin: np.logical_or(randomoffset(x_arcmin,y_arcmin), 
                                                             offsetpointing(x_arcmin,y_arcmin))

centerandoffset = lambda x_arcmin,y_arcmin: np.logical_or(acscentered(x_arcmin,y_arcmin), 
                                                          offsetpointing(x_arcmin,y_arcmin))

centerandrandoffset = lambda x_arcmin,y_arcmin: np.logical_or(acscentered(x_arcmin,y_arcmin), 
                                                              randomoffset(x_arcmin,y_arcmin))


acsdiag = np.sqrt(2.)*3.2/2.
squaremosaic = lambda x_arcmin,y_arcmin: np.logical_or(np.logical_or(acsmask(x_arcmin,y_arcmin,0., acsdiag),
                                                                     acsmask(x_arcmin,y_arcmin,0.,-acsdiag)), 
                                      np.logical_or(acsmask(x_arcmin,y_arcmin,-acsdiag,0.),
                                                    acsmask(x_arcmin,y_arcmin,acsdiag,0.)))

def selectMask(config):

    maskcase = {'squaremask' : lambda x_arcmin, y_arcmin: squaremask(x_arcmin, y_arcmin, 
                                                                     config.maskx, 
                                                                     config.masky, 
                                                                     config.masktheta, 
                                                                     config.masksidelength),
                'circlemask' : lambda x_arcmin, y_arcmin: circlemask(x_arcmin, y_arcmin,
                                                                     config.maskx, 
                                                                     config.masky, 
                                                                     config.maskrad),
                'acsmask' : lambda x_arcmin, y_arcmin: acsmask(x_arcmin, y_arcmin,
                                                               config.maskx, config.masky),
                'wfc3mask' : lambda x_arcmin, y_arcmin: wfc3mask(x_arcmin, y_arcmin,
                                                                 config.maskx, config.masky),
                'acscentered' : acscentered,
                'offsetpointing' : offsetpointing,
                'rotatedoffset' : rotatedoffset,
                'offsetmosaic' : offsetmosaic,
                'offset3' : offset3,
                'pisco3' : pisco3,
                'pisco4' : pisco4,
                'randomoffset' : randomoffset,
                'randomoffsetmosaic' : randomoffsetmosaic,
                'centerandoffset' : centerandoffset,
                'centerandrandoffset' : centerandrandoffset,
                'squaremosaic' : squaremosaic}



class FoVPicker(GalaxyPicker):

    def mask(self, sim):

        x_arcmin = sim.x_arcmin
        y_arcmin = sim.y_arcmin
        zcluster = sim.zcluster

        if 'targetz' in self.config:
            #adjust for different redshift
            curr_angdist = nfwutils.global_cosmology.angulardist(zcluster)
            newangdist = nfwutils.global_cosmology.angulardist(self.config.targetz)
            ratio = curr_angdist/newangdist
            x_arcmin = x_arcmin*ratio
            y_arcmin = y_arcmin*ratio

        selected = selectMask(self.config)(x_arcmin, y_arcmin)
        
        return selected
