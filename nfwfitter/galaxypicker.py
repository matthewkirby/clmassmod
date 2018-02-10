'''
Select Galaxies from a simulation for analysis
'''

############

import numpy as np

import utils.simutils as simutils
import nfwutils

###########


##################

class GalaxyPicker(object):

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

    def configure(self, config):

        self.nperarcmin = config['nperarcmin']

    def mask(self, sim):

        x_arcmin = sim.x_arcmin
        y_arcmin = sim.y_arcmin

        targetdensity = self.nperarcmin        

        max_x = np.max(x_arcmin)
        min_x = np.min(x_arcmin)
        delta_x = max_x - min_x

        max_y = np.max(y_arcmin)
        min_y = np.min(y_arcmin)
        delta_y = max_y - min_y

        area = delta_x*delta_y

        targetnumber = targetdensity*area        
        if targetdensity == -1: # take all
            targetnumber = len(x_arcmin)



        availablenumber = len(x_arcmin)

        if targetnumber > availablenumber:
            raise InsufficientGalaxiesException

        accept = float(targetnumber) / availablenumber

        randomthrow = np.random.random(len(x_arcmin))

        selected = randomthrow < accept

        return selected

################################################
#####################

class FoVPicker(GalaxyPicker):

    def configure(self, config):
        pass

    def mask(self, sim):

        x_arcmin = sim.x_arcmin
        y_arcmin = sim.y_arcmin

        selected = self._mask(x_arcmin, y_arcmin)
        
        return selected


###########

class SquareMask(FoVPicker):

    def configure(self, config):

        self.x=0.
        self.y=0.
        self.theta=0.
        self.sidelength=1.
        if 'maskx' in config:
            #x,y in arcminutes
            self.x = config['maskx']
        if 'masky' in config:
            self.y = config['masky']
        if 'masktheta' in config:
            self.theta = config['masktheta']
        if 'masksidelength' in config:
            self.sidelength = config['masksidelength']

    def _mask(self, x_arcmin, y_arcmin, x = None, y = None, theta = None, sidelength = None):

        if x is None:
            x = self.x
        if y is None:
            y = self.y
        if theta is None:
            theta = self.theta
        if sidelength is None:
            sidelength = self.sidelength

        theta_rad = np.pi*theta/180.

        dX = x_arcmin - x   
        dY = y_arcmin - y

        rdX = np.cos(theta_rad)*dX - np.sin(theta_rad)*dY
        rdY = np.sin(theta_rad)*dX + np.cos(theta_rad)*dY

        return (np.abs(rdX) + np.abs(rdY)) <= (np.sqrt(2)*sidelength/2.)

###

class RectangleMask(FoVPicker):

    def configure(self, config):



        self.x = 0.
        self.y = 0.
        self.theta=0.
        self.xlength=1.
        self.ylength=1.

        if 'maskx' in config:
            self.x = config['maskx']
        if 'masky' in config:
            self.y = config['masky']
        if 'maskxlength' in config:
            self.xlength = config['maskxlength']
        if 'maskylength' in config:
            self.ylength = config['maskylength']

    #x,y in arcminutes

    def _mask(x_arcmin, y_arcmin, x = None, y = None, theta = None, xlength = None, ylength = None):

        if x is None:
            x = self.x
        if y is None:
            y = self.y
        if theta is None:
            theta = self.theta
        if xlength is None:
            xlength = self.xlength
        if ylength is None:
            ylength = self.ylength


        theta_rad = np.pi*theta/180.

        dX = x_arcmin - x   
        dY = y_arcmin - y

        rdX = np.cos(theta_rad)*dX - np.sin(theta_rad)*dY
        rdY = np.sin(theta_rad)*dX + np.cos(theta_rad)*dY

        return np.logical_and(np.abs(rdX) < xlength/2., np.abs(rdY) < ylength/2.)

###

class CircleMask(FoVPicker):

    def configure(self, config):

        self.x = 0.
        self.y = 0.
        self.rad = 1.

        if 'maskx' in config:
            self.x = config['maskx']
        if 'masky' in config:
            self.y = config['masky']
        if 'maskrad' in config:
            self.rad = config['maskrad']



    def _mask(x_arcmin, y_arcmin, x = None, y = None, rad = None):
    #arcmnutes

        if x is None:
            x = self.x
        if y is None:
            y = self.y
        if rad is None:
            rad = self.rad

        dX = x_arcmin - x
        dY = y_arcmin - y

        return np.sqrt(dX**2 + dY**2) < rad

###

class ACSMask(SquareMask):

    def configure(self, config):

        super(ACSMask, self).configure(config)
        self.sidelength = 3.2
        



###
        
class Wfc3Mask(SquareMask):

    def configure(self, config):
        super(Wfc3Mask, self).configure(config)
        self.theta = 45.
        self.sidelength = 3.2*4./5.


###

class ACSCenteredMask(FoVPicker):

    def configure(self, config):
        self.acsmask = ACSMask()
        self.wfc3mask = Wfc3Mask()
        
        self.acsmask.configure({})
        self.wfc3mask.configure({})


    def _mask(self, x_arcmin,y_arcmin):

        return np.logical_or(self.acsmask._mask(x_arcmin,y_arcmin,x=0,y=0), 
                             self.wfc3mask._mask(x_arcmin,y_arcmin,x=0, y=6.))


###

class OffsetPointing(FoVPicker):

    def configure(self, config):
        self.acsmask = ACSMask()
        self.wfc3mask = Wfc3Mask()
        
        self.acsmask.configure({})
        self.wfc3mask.configure({})
        

    def _mask(self, x_arcmin,y_arcmin):
        return np.logical_or(self.acsmask._mask(x_arcmin,y_arcmin,x=0.,y=-3.), 
                             self.wfc3mask._mask(x_arcmin,y_arcmin,x=0., y=3.0))
        

###

class RotatedOffset(FoVPicker):

    def configure(self, config):
        self.acsmask = ACSMask()
        self.wfc3mask = Wfc3Mask()
        
        self.acsmask.configure({})
        self.wfc3mask.configure({})


    def _mask(self, x_arcmin,y_arcmin): 
        
        return np.logical_or(self.acsmask._mask(x_arcmin,y_arcmin,x=-3.,y=0.), 
                             self.wfc3mask._mask(x_arcmin,y_arcmin,x=3., y=0.0))

###

class OffsetMosaic(FoVPicker):
    
    def configure(self, config):

        self.offsetpointing = OffsetPointing()
        self.rotatedoffset = RotatedOffset()

        self.offsetpointing.configure({})
        self.rotatedoffset.configure({})

        

    def _mask(self, x_arcmin,y_arcmin):

        return np.logical_or(self.offsetpointing._mask(x_arcmin,y_arcmin), 
                             self.rotatedoffset._mask(x_arcmin,y_arcmin))

###

class Offset3(FoVPicker):

    def configure(self, config):
        self.offsetmosaic = OffsetMosaic()
        self.acscentered = ACSCenteredMask()

        self.offsetmosaic.configure({})
        self.acscentered.configure({})

    def _mask(self, x_arcmin,y_arcmin): 

        return np.logical_or(self.offsetmosaic._mask(x_arcmin,y_arcmin), 
                             self.acscentered._mask(x_arcmin,y_arcmin))

###

class Pisco3(FoVPicker):

    def configure(self, config):

        self.rectanglemask = RectangleMask()
        self.rectanglemask.configure({})

    def _mask(self, x_arcmin,y_arcmin): 


        return np.logical_or(np.logical_or(rectanglemask._mask(x_arcmin,y_arcmin,x=0,y=0,theta=0,xlength=8,ylength=6), 
                                           rectanglemask._mask(x_arcmin,y_arcmin,x=0,y=7,theta=0,xlength=6,ylength=8)), 
                             rectanglemask._mask(x_arcmin,y_arcmin,x=0,y=-7,theta=0,xlength=6,ylength=8))

###

class Pisco4(FoVPicker):

    def configure(self, config):
        self.rectanglemask = RectangleMask()
        self.rectanglemask.configure({})

    def _mask(self, x_arcmin,y_arcmin): 

        
        return np.logical_or(np.logical_or(self.rectanglemask._mask(x_arcmin,y_arcmin,x=-5.5, y=0, theta=0, xlength=8, ylength=6), 
                                           self.rectanglemask._mask(x_arcmin,y_arcmin,x=5.5, y=0, theta=0, xlength=8, ylength=6)), 
                             np.logical_or(self.rectanglemask._mask(x_arcmin,y_arcmin,x=0, y=5.5, theta=0, xlength=6, ylength=8), 
                                           self.rectanglemask._mask(x_arcmin,y_arcmin,x=0, y=-5.5, theta=0,xlength=6, ylength=8)))

###

class RandomOffset(FoVPicker):
    '''Not sure this does what it is supposed to do'''

    def configure(self, config):
        self.acsmask = ACSMask()
        self.wfc3mask = Wfc3Mask()

        self.acsmask.configure({})
        self.wfc3mask.configure({})

        

    def _mask(self, x_arcmin,y_arcmin):
        posangle = np.random.uniform(0, 2*np.pi)
        rotmatrix = np.array([[np.cos(posangle), -np.sin(posangle)],
                              [np.sin(posangle), np.cos(posangle)]])
        acspos = np.dot(rotmatrix, np.array([0., -3]))
        wfc3pos = np.dot(rotmatrix, np.array([0., 3]))


        rotatedoffset = np.logical_or(self.acsmask._mask(x_arcmin,y_arcmin,x=acspos[0], y=acspos[1]), 
                                      self.wfc3mask._mask(x_arcmin,y_arcmin,x=wfc3pos[0], y=wfc3pos[1]))
        return rotatedoffset

###

class RandomOffsetMosaic(FoVPicker):

    def configure(self, config):
        self.randomoffset = RandomOffset()
        self.randomoffset.configure({})

        self.offsetpointing = OffsetPointing()
        self.offsetpointing.configure({})

    def _mask(self, x_arcmin,y_arcmin): 

        

        return np.logical_or(self.randomoffset._mask(x_arcmin,y_arcmin), 
                             self.offsetpointing._mask(x_arcmin,y_arcmin))

###

class CenterAndOffset(FoVPicker):

    def configure(self, config):
        self.acscentered = ACSCenteredMask()
        self.offsetpointing = OffsetPointing()

        self.acscentered.configure({})
        self.offsetpointing.configure({})


    def _mask(self, x_arcmin,y_arcmin): 

        return np.logical_or(self.acscentered._mask(x_arcmin,y_arcmin), 
                             self.offsetpointing._mask(x_arcmin,y_arcmin))

###

class CenterAndRandOffset(FoVPicker):

    def configure(self, config):

        self.acscentered = ACSCenteredMask()
        self.randomoffset = RandomOffset()

        self.acscentered.configure({}) 
        self.randomoffset.configure({}) 
        

    def _mask(self, x_arcmin,y_arcmin): 
        
        return np.logical_or(self.acscentered._mask(x_arcmin,y_arcmin), 
                             self.randomoffset._mask(x_arcmin,y_arcmin))

###


acsdiag = np.sqrt(2.)*3.2/2.

class SquareMosaic(FoVPicker):

    def configure(self, config):

        self.acsmask = ACSMask()
        self.acsmask.configure({})


    def _mask(self, x_arcmin,y_arcmin): 


        return np.logical_or(np.logical_or(self.acsmask._mask(x_arcmin,y_arcmin,x=0., y=acsdiag),
                                           self.acsmask._mask(x_arcmin,y_arcmin,x=0.,y=-acsdiag)), 
                             np.logical_or(self.acsmask._mask(x_arcmin,y_arcmin,x=-acsdiag,y=0.),
                                           self.acsmask._mask(x_arcmin,y_arcmin,x=acsdiag,y=0.)))


