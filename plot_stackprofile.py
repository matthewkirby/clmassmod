##############################
# Overplot the predicted NFW halo over the stack
##############################

import pylab
import nfwmodeltools
import nfwutils
import readtxtfile, glob
import ldac
import numpy as np
import readMXXL

###################

massedges = np.array([0, 3.8e14, 4.2e14, 4.9e14, 5e15])
concenedges = np.array([0, 0.2, 0.26, 0.38, 4.38, 10])


def residual(binbase):

    mass, concen =  readtxtfile.readtxtfile('%s.dat' % binbase)[0]

    cat = ldac.openObjectFile('%s.cat' % binbase)

    zlens = cat.hdu.header['ZLENS']
    
    rscale = nfwutils.rscaleConstM(mass/nfwutils.global_cosmology.h, concen, zlens, 200)
    

    gamma = nfwmodeltools.NFWShear(cat['r_mpc'], concen, rscale, zlens)
    kappa = nfwmodeltools.NFWKappa(cat['r_mpc'], concen, rscale, zlens)

    gpred = cat['beta_s']*gamma / (1 - (cat['beta_s2']*kappa/cat['beta_s']))


    fig = pylab.figure()
    pylab.errorbar(cat['r_mpc'], cat['ghat']/gpred, cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])*gpred), fmt='bo')
    pylab.axhline(1.0, c='k', linewidth=2)
    pylab.xlabel('Radius [Mpc]', fontsize=16)
    pylab.ylabel('g_meas / g_pred', fontsize=16)
    pylab.title('Redshift=%1.1f Mass=%1.2fx10^14 Concen=%1.2f' % (zlens, mass/1e14, concen))

    pylab.savefig('%s.png' % binbase)

    return fig

    
def run(stackdir, simreader):

    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    stackfiles = glob.glob('%s/*.cat' % stackdir)
    stackbases = [x[:-4] for x in stackfiles]
    figs = []
    for stackbase in stackbases:
        try:
            figs.append(residual(stackbase))
        except ValueError:
            pass
                       

    return figs

    
    

    
    

    

    



