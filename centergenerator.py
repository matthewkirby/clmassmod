'''
Distributions of where the profile center lands
'''

#############

import numpy as np

import astropy.io.ascii as asciireader
import readtxtfile
import nfwutils
import voigt_tools as vt

#############

def CenterGeneratorFactory(config):

    if 'centergenerator' in config:
        cgmodule, cgclass = config['centergenerator'].split(':')
        centergenerator = simutils.buildObject(cgmodule, cgclass, config = config)
    else:
        centergenerator = NoOffset()

    return centergenerator


#######

class NoOffset(object):

    def __init__(self, config):
        pass
        
    def __call__(self, sim):
        return 0., 0.

########

class CenterGenerator(object):

    def __init__(self, config):
        self.config = config

    def __call__(self, sim):

        return self.offset(sim, self.config)

############


class SZSimOffset(CenterGenerator):

    def __init__(self, *args, **kwds):
        super(SZSimOffset, self).__init__(*args, **kwds)
        self.szsim_offsetcat = asciireader.read('{}/SPT_SN_offset.dat'.format(self.config['offsetdistro_dir']))

    def offset(sim, config):

        offsetingcoresize = self.szsim_offsetcat[self.szsim_offsetcat['coresize[arcmin]'] == config.coresize]

        selectedsim =np.random.uniform(0, len(offsetingcoresize))

        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    
        targetDl = nfwutils.global_cosmology.angulardist(config.targetz)
        anglescale = targetDl/dL  #account for the fact that the fixed angular scatter turns into different effective r_mpc scatter
        offsetx = anglescale*(offsetingcoresize['peak_xpix[arcmin]'] - offsetingcoresize['cluster_xpix'])[selectedsim]  #arcmin
        offsety = anglescale*(offsetingcoresize['peak_ypix'] - offsetingcoresize['cluster_ypix'])[selectedsim]


        offset_phi = np.random.uniform(0, 2*np.pi)

        newoffsetx = offsetx*np.cos(offset_phi) - offsety*np.sin(offset_phi)
        newoffsety = offsetx*np.sin(offset_phi) + offsety*np.cos(offset_phi)

        return newoffsetx, newoffsety


#######


class SZLensingPeakOffset(CenterGenerator):

    def offset(sim, config):

        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    
        targetDl = nfwutils.global_cosmology.angulardist(config.targetz)

        scatter = 0.237*targetDl/dL #arcmin, scaled
        centeroffsetx, centeroffsety = scatter*np.random.standard_normal(2)

        return centeroffsetx, centeroffsety


####


class SZXVPTheoryOffset(CenterGenerator):

    def __init__(self, *args, **kwds):
        super(SZXVPTheoryOffset, self).__init__(*args, **kwds)
        self.xvpoffset = XrayXVPOffset(self.config)

    def offset(sim, config):

        #physical scatter in arcmin, approp for target redshift
        xvp_offsetx, xvp_offsety = self.xvpoffset(sim)


        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    
        targetDl = nfwutils.global_cosmology.angulardist(config.targetz)

        sz_noisescatter = 0.3*targetDl/dL #arcmin, scaled
        centeroffsetx, centeroffsety = sz_noisescatter*np.random.standard_normal(2)

        return (centeroffsetx + xvp_offsetx, 
                centeroffsety + xvp_offsety)


####



class SZXVPBCGOffset(CenterGenerator):

    def __init__(self, *args, **kwds):
        super(SZXVPBCGOffset, self).__init__(*args, **kwds)
        self.sz_xvp_bcg_offsets_deg = readtxtfile.readtxtfile('{}/sptxvp_bcgsz'.format(self.config['offsetdistro_dir']))[:,1]

    def offset(sim, config):


        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    
        targetDl = nfwutils.global_cosmology.angulardist(config.targetz)
        anglescale = targetDl/dL  #account for the fact that the fixed angular scatter turns into different effective r_mpc scatter

        offset = 60*anglescale*sz_xvp_bcg_offsets_deg[np.random.randint(0, len(sz_xvp_bcg_offsets_deg), 1)]

        offset_phi = np.random.uniform(0, 2*np.pi)

        newoffsetx = offset*np.cos(offset_phi)
        newoffsety = offset*np.sin(offset_phi)

        return newoffsetx, newoffsety


####

class SZAnalytic(CenterGenerator):

    def offset(sim, config):

        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    
        targetDl = nfwutils.global_cosmology.angulardist(config.targetz)
        anglescale = targetDl/dL  #account for the fact that the fixed angular scatter turns into different effective r_mpc scatter

        sz_noisescatter = anglescale*np.sqrt(config.szbeam**2 + config.coresize**2)/config.sz_xi

        centeroffsetx, centeroffsety = sz_noisescatter*np.random.standard_normal(2)

        return (centeroffsetx, 
                centeroffsety)
    
####    




class XrayWTGOffset(CenterGenerator):

    def __init__(self, *args, **kwds):
        super(XrayWTGOffset, self).__init__(*args, **kwds)
        self.wtg_offsets_mpc = [x[0] for x in readtxtfile.readtxtfile('{}/wtg_offsets.dat'.format(self.config['offsetdistro_dir']))]

    def offset(sim, config):


        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    

        radial_offset_mpc = self.wtg_offsets_mpc[np.random.randint(0, len(self.wtg_offsets_mpc), 1)]
        radial_offset_arcmin = (radial_offset_mpc/(dL))*(180./np.pi)*60.
        phi_offset = np.random.uniform(0, 2*np.pi)
        centeroffsetx = radial_offset_arcmin*np.cos(phi_offset)
        centeroffsety = radial_offset_arcmin*np.sin(phi_offset)

        return centeroffsetx, centeroffsety


###

class XraySPTHSTOffset(CenterGenerator):

    def offset(sim, config):



        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    

        #offset distribution simple log delta r ~ N(mu, sig) fit to SPT-HST xray bcg offset distro (from Inon)
        centeroffset_mpc = np.exp(-2.625 + 1.413*np.random.standard_normal())
        offset_radial = (centeroffset_mpc/dL)*(180./np.pi)*60.
        offset_phi = np.random.uniform(0, 2*np.pi)

        centeroffsetx = offset_radial*np.cos(offset_phi)
        centeroffsety = offset_radial*np.sin(offset_phi)

        return centeroffsetx, centeroffsety

###



class XrayXVPOffset(CenterGenerator):

    def __init__(self, *args, **kwds):
        super(XrayXVPOffset, self).__init__(*args, **kwds)
        self.xvp_offsets_mpc = readtxtfile.readtxtfile('{}/sptxvp_bcgxray'.format(self.config['offsetdistro_dir']))[:,0]


    def offset(sim, config):


        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    


        radial_offset_mpc = self.xvp_offsets_mpc[np.random.randint(0, len(self.xvp_offsets_mpc), 1)][0]
        radial_offset_arcmin = (radial_offset_mpc/(dL))*(180./np.pi)*60.
        phi_offset = np.random.uniform(0, 2*np.pi)
        centeroffsetx = radial_offset_arcmin*np.cos(phi_offset)
        centeroffsety = radial_offset_arcmin*np.sin(phi_offset)

        return centeroffsetx, centeroffsety


###


class XrayCCCPOffset(CenterGenerator):

    def __init__(self, *args, **kwds):
        super(XrayCCCPOffset, self).__init__(*args, **kwds)
        self.offsets_kpc = [x[0] for x in readtxtfile.readtxtfile('{}/cccp_offsets.dat'.format(self.config['offsetdistro_dir']))]

    def offset(sim, config):

        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    



        radial_offset_kpc = self.offsets_kpc[np.random.randint(0, len(self.offsets_kpc), 1)]
        radial_offset_arcmin = (radial_offset_kpc/(1000.*dL))*(180./np.pi)*60.
        phi_offset = np.random.uniform(0, 2*np.pi)
        centeroffsetx = radial_offset_arcmin*np.cos(phi_offset)
        centeroffsety = radial_offset_arcmin*np.sin(phi_offset)

        return centeroffsetx, centeroffsety


###


class XrayLensingPeakOffset(CenterGenerator):

    def offset(sim, config):

        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    

        delta_mpc = 0.107*np.random.standard_normal(2)

        centeroffsetx, centeroffsety = (delta_mpc/dL)*(180./np.pi)*60 #arcmin

        return centeroffsetx, centeroffsety

###

class XrayLensingPeakVoigtOffset(CenterGenerator):

    def offset(sim, config):

        print 'Xray Lensing Peak Voigt'

        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    

        delta_mpc = vt.voigtSamples(0.048, 0.0565, 2, limits=(-0.3, 0.3))

        centeroffsetx, centeroffsety = (delta_mpc/dL)*(180./np.pi)*60 #arcmin

        return centeroffsetx, centeroffsety

###



class XrayMagneticumOffset(CenterGenerator):

    def __init__(self, *args, **kwds):
        super(XrayMagneticumOffset, self).__init__(*args, **kwds)
        self.xray_magneticum_distro = asciireader.read('{}/magneticum_offsets.dat'.format(self.config['offsetdistro_dir']))


    def offset(sim, config):

        dL = nfwutils.global_cosmology.angulardist(sim.zcluster)    

        delta_kpc = self.xray_magneticum_distro['xrayoffset'][np.random.randint(0, len(self.xray_magneticum_distro), 1)][0]

        radial_offset_arcmin = (delta_kpc/(1000*dL))*(180./np.pi)*60 #arcmin
        phi_offset = np.random.uniform(0, 2*np.pi)

        centeroffsetx = radial_offset_arcmin*np.cos(phi_offset)
        centeroffsety = radial_offset_arcmin*np.sin(phi_offset)



        return centeroffsetx, centeroffsety


    

###


