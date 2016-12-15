'''
Distributions of where the profile center lands
'''

#############

import numpy as np
import pkg_resources

import astropy.io.ascii as asciireader
import readtxtfile
import nfwutils
import voigt_tools as vt
import nfwutils

#############


class NoOffset(object):
        
    def __call__(self, sim):
        return 0., 0.

########

class SZSimOffset(object):

    def __init__(self):

        szsim_offsetcat = asciireader.read(pkg_resources.resource_stream('nfwfitter', 'data/SPT_SN_offset.dat'))

        self.szsim_offsetcat = szsim_offsetcat[szsim_offsetcat['SN'] >= 5]

    def configure(self, config):

        self.coresize = config['coresize']

    def __call__(self, sim):

        offsetingcoresize = self.szsim_offsetcat[self.szsim_offsetcat['coresize[arcmin]'] == self.coresize]

        selectedsim =np.random.randint(0, len(offsetingcoresize))

        offsetx = (offsetingcoresize['peak_xpix[arcmin]'] - offsetingcoresize['cluster_xpix'])[selectedsim]  #arcmin
        offsety = (offsetingcoresize['peak_ypix'] - offsetingcoresize['cluster_ypix'])[selectedsim]


        offset_phi = np.random.uniform(0, 2*np.pi)

        newoffsetx = offsetx*np.cos(offset_phi) - offsety*np.sin(offset_phi)
        newoffsety = offsetx*np.sin(offset_phi) + offsety*np.cos(offset_phi)

        return newoffsetx, newoffsety


class SZSimOffsetCoreIgnored(object):

    def __init__(self):

        szsim_offsetcat = asciireader.read(pkg_resources.resource_stream('nfwfitter', 'data/SPT_SN_offset.dat'))


        self.szsim_offsetcat = szsim_offsetcat[szsim_offsetcat['SN'] >= 5]

    def configure(self, config):

        pass

    def __call__(self, sim):

        selectedsim = np.random.randint(0, len(self.szsim_offsetcat))

        offsetx = (self.szsim_offsetcat['peak_xpix[arcmin]'] - self.szsim_offsetcat['cluster_xpix'])[selectedsim]  #arcmin
        offsety = (self.szsim_offsetcat['peak_ypix'] - self.szsim_offsetcat['cluster_ypix'])[selectedsim]


        offset_phi = np.random.uniform(0, 2*np.pi)

        newoffsetx = offsetx*np.cos(offset_phi) - offsety*np.sin(offset_phi)
        newoffsety = offsetx*np.sin(offset_phi) + offsety*np.cos(offset_phi)

        return newoffsetx, newoffsety


#######


class SZLensingPeakOffset(object):


    def __call__(self, sim):

        scatter = 0.237
        centeroffsetx, centeroffsety = scatter*np.random.standard_normal(2)

        return centeroffsetx, centeroffsety


####


class SZXVPTheoryOffset(object):

    def __init__(self):
        self.xvpoffset = XrayXVPOffset()


    def __call__(self, sim):

        #physical scatter in arcmin, approp for target redshift
        xvp_offsetx, xvp_offsety = self.xvpoffset(sim)

        sz_noisescatter = 0.3
        centeroffsetx, centeroffsety = sz_noisescatter*np.random.standard_normal(2)

        return (centeroffsetx + xvp_offsetx, 
                centeroffsety + xvp_offsety)


####



class SZXVPBCGOffset(object):

    def __init__(self):

        self.sz_xvp_bcg_offsets_deg = readtxtfile.readtxtfile(pkg_resources.resource_filename('nfwfitter', 'data/sptxvp_bcgsz'))[:,1]


    def __call__(self, sim):



        offset = 60*sz_xvp_bcg_offsets_deg[np.random.randint(0, len(sz_xvp_bcg_offsets_deg), 1)]

        offset_phi = np.random.uniform(0, 2*np.pi)

        newoffsetx = offset*np.cos(offset_phi)
        newoffsety = offset*np.sin(offset_phi)

        return newoffsetx, newoffsety


####

class SZAnalytic(object):

    def configure(self, config):

        self.szbeam = config['szbeam']
        self.coresize = config['coresize']
        self.sz_xi = config['sz_xi']

    def __call__(self, sim):


        sz_noisescatter = np.sqrt(self.szbeam**2 + self.coresize**2)/self.sz_xi

        offset = sz_noisescatter*np.random.standard_normal()

        offset_phi = np.random.uniform(0, 2*np.pi)

        newoffsetx = offset*np.cos(offset_phi)
        newoffsety = offset*np.sin(offset_phi)

        return newoffsetx, newoffsety


    
####    




class XrayWTGOffset(object):

    def __init__(self):

        self.wtg_offsets_mpc = [x[0] for x in readtxtfile.readtxtfile(pkg_resources.resource_filename('nfwfitter', 'data/wtg_offsets.dat'))]


    def __call__(self, sim):


        dL = nfwutils.global_cosmology.angulardist(sim.zlens)

        radial_offset_mpc = self.wtg_offsets_mpc[np.random.randint(0, len(self.wtg_offsets_mpc), 1)]
        radial_offset_arcmin = (radial_offset_mpc/(dL))*(180./np.pi)*60.
        phi_offset = np.random.uniform(0, 2*np.pi)
        centeroffsetx = radial_offset_arcmin*np.cos(phi_offset)
        centeroffsety = radial_offset_arcmin*np.sin(phi_offset)

        return centeroffsetx, centeroffsety


###

class XraySPTHSTOffset(object):

    def __call__(self, sim):

        dL = nfwutils.global_cosmology.angulardist(sim.zlens)    

        #offset distribution simple log delta r ~ N(mu, sig) fit to SPT-HST xray bcg offset distro (from Inon)
        centeroffset_mpc = np.exp(-2.625 + 1.413*np.random.standard_normal())
        offset_radial = (centeroffset_mpc/dL)*(180./np.pi)*60.
        offset_phi = np.random.uniform(0, 2*np.pi)

        centeroffsetx = offset_radial*np.cos(offset_phi)
        centeroffsety = offset_radial*np.sin(offset_phi)

        return centeroffsetx, centeroffsety

###


class XrayXVPOffset(object):

    def __init__(self):

        self.xvp_offsets_mpc = readtxtfile.readtxtfile(pkg_resources.resource_filename('nfwfitter', 'data/sptxvp_bcgxray'))[:,0]


    def __call__(self, sim):


        dL = nfwutils.global_cosmology.angulardist(sim.zlens)    


        radial_offset_mpc = self.xvp_offsets_mpc[np.random.randint(0, len(self.xvp_offsets_mpc), 1)][0]
        radial_offset_arcmin = (radial_offset_mpc/(dL))*(180./np.pi)*60.
        phi_offset = np.random.uniform(0, 2*np.pi)
        centeroffsetx = radial_offset_arcmin*np.cos(phi_offset)
        centeroffsety = radial_offset_arcmin*np.sin(phi_offset)

        return centeroffsetx, centeroffsety


###


class XrayCCCPOffset(object):

    def __init__(self):

        self.offsets_kpc = [x[0] for x in readtxtfile.readtxtfile(pkg_resources.resource_filename('nfwfitter', 'data/cccp_offsets.dat'))]

    def __call__(self, sim):

        dL = nfwutils.global_cosmology.angulardist(sim.zlens)    



        radial_offset_kpc = self.offsets_kpc[np.random.randint(0, len(self.offsets_kpc), 1)]
        radial_offset_arcmin = (radial_offset_kpc/(1000.*dL))*(180./np.pi)*60.
        phi_offset = np.random.uniform(0, 2*np.pi)
        centeroffsetx = radial_offset_arcmin*np.cos(phi_offset)
        centeroffsety = radial_offset_arcmin*np.sin(phi_offset)

        return centeroffsetx, centeroffsety


###


class XrayLensingPeakOffset(object):

    def __call__(self, sim):

        dL = nfwutils.global_cosmology.angulardist(sim.zlens)    

        delta_mpc = 0.107*np.random.standard_normal(2)

        centeroffsetx, centeroffsety = (delta_mpc/dL)*(180./np.pi)*60 #arcmin

        return centeroffsetx, centeroffsety

###

class XrayLensingPeakVoigtOffset(object):

    def __call__(self, sim):

        print 'Xray Lensing Peak Voigt'

        dL = nfwutils.global_cosmology.angulardist(sim.zlens)    

        delta_mpc = vt.voigtSamples(0.048, 0.0565, 2, limits=(-0.3, 0.3))

        centeroffsetx, centeroffsety = (delta_mpc/dL)*(180./np.pi)*60 #arcmin

        return centeroffsetx, centeroffsety

###



class XrayMagneticumOffset(object):

    def __init__(self):

        self.xray_magneticum_distro = asciireader.read(pkg_resources.resource_filename('nfwfitter', 'data/xray_offsets_may2016.txt'))


    def __call__(self, sim):

        dL = nfwutils.global_cosmology.angulardist(sim.zlens)    

        delta_kpc = self.xray_magneticum_distro['xray_offset_kpc'][np.random.randint(0, len(self.xray_magneticum_distro), 1)][0]

        radial_offset_arcmin = (delta_kpc/(1000*dL))*(180./np.pi)*60 #arcmin
        phi_offset = np.random.uniform(0, 2*np.pi)

        centeroffsetx = radial_offset_arcmin*np.cos(phi_offset)
        centeroffsety = radial_offset_arcmin*np.sin(phi_offset)



        return centeroffsetx, centeroffsety


    

###


