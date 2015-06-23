############
# Show miscentering distributions for publication
############


import publication_plots as pp
import pylab, matplotlib
import nfwfit
import nfwutils
import astroutils
import numpy as np
import astropy.io.ascii as asciireader

####

def plotcdf(sample, ax, *args, **keywords):
    sorted_sample = np.sort(sample)
    xs = np.hstack([[sorted_sample[0]], sorted_sample])
    ys = np.hstack([[0], (np.arange(len(sample))+1)/float(len(sample))])
    ax.plot(xs, ys , *args, **keywords)

####

sptdat = asciireader.read('sptdat')

###

def makeXrayPlot():

    fig = pylab.figure()
    ax = fig.add_subplot(1,1,1)

    #xrayxvp

    xrayxvpoffsets = 1000*nfwfit.xvp_offsets_mpc
    plotcdf(xrayxvpoffsets, ax, label='XVP X-ray-BCG', linewidth=3, color=pp.colors[1])

    #xraylensing

    dLs = np.array([nfwutils.global_cosmology.angulardist(z) for z in sptdat['z_l']])
    xray_lensing_offsets = 1000*dLs*(np.pi/180.)*astroutils.greatCircleDistance(sptdat['xray_ra'], sptdat['xray_dec'], sptdat['lensing_ra'], sptdat['lensing_dec'])

    plotcdf(xray_lensing_offsets, ax, label='SPT13 X-ray-Lensing', linewidth=1.5, color=pp.colors[0])

    fit_deltas_kpc = 0.107*np.random.standard_normal((2, 2000))*1000
    fit_offsets = np.sqrt(fit_deltas_kpc[0]**2 + fit_deltas_kpc[1]**2)
    
    plotcdf(fit_offsets, ax, c='k', linestyle='--', label='Gaussian Model')

    ax.set_xlim(1, 1000)
    ax.set_xscale('log')
    ax.set_xlabel(r'Offset $\Delta$R [kpc]')
    ax.set_ylabel(r'P($\le\Delta$R)')
    ax.legend(loc='upper left', numpoints=1, frameon=False)

    fig.tight_layout()

    fig.savefig('figures/xrayoffset.png')
    fig.savefig('figures/xrayoffset.eps')
    fig.savefig('figures/xrayoffset.ps')
    fig.savefig('figures/xrayoffset.pdf')

    return fig

####

def makeSZPlot():

    fig = pylab.figure()
    ax = fig.add_subplot(1,1,1)

    #hydro sims

    hydrocat = nfwfit.szsim_offsetcat
    hydro_offset = np.sqrt((hydrocat['peak_xpix[arcmin]'] - hydrocat['cluster_xpix'])**2 + \
                           (hydrocat['peak_ypix'] - hydrocat['cluster_ypix'])**2)

    plotcdf(hydro_offset, ax, label='Hydro Sims', linewidth=3, color=pp.colors[1])

    #lensing sz offsets

    sz_lensing = 60*astroutils.greatCircleDistance(sptdat['sz_ra'], sptdat['sz_dec'], sptdat['lensing_ra'], sptdat['lensing_dec'])

    plotcdf(sz_lensing, ax, label='SPT13 SZ-Lensing', linewidth=1.5, color=pp.colors[0])

    fit_deltas = 0.237*np.random.standard_normal((2, 2000))
    fit_offsets = np.sqrt(fit_deltas[0]**2 + fit_deltas[1]**2)
    
    plotcdf(fit_offsets, ax, c='k', linestyle='--', label='Gaussian Model')

    ax.set_xlim(.01, 10)
    ax.set_xscale('log')
    ax.set_xlabel(r'Offset $\Delta$R [arcmin]')
    ax.set_ylabel(r'P($\le\Delta$R)')
    ax.legend(loc='upper left', numpoints=1, frameon=False)

    fig.tight_layout()

    fig.savefig('figures/szoffset.png')
    fig.savefig('figures/szoffset.eps')
    fig.savefig('figures/szoffset.ps')
    fig.savefig('figures/szoffset.pdf')

    return fig

    

    
