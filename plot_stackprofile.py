##############################
# Overplot the predicted NFW halo over the stack
##############################

import pylab, matplotlib
import nfwmodeltools
import nfwutils
import readtxtfile, glob
import ldac
import numpy as np
#import readMXXL, readBCC

###################

massedges = np.array([0, 3.8e14, 4.2e14, 4.9e14, 5e15])
concenedges = np.array([0, 0.2, 0.26, 0.38, 4.38, 10])


def residual(binbase):

    mass, concen, redshift =  readtxtfile.readtxtfile('%s.dat' % binbase)[0]

    cat = ldac.openObjectFile('%s.cat' % binbase)

    zlens = cat.hdu.header['ZLENS']
    
    rscale = nfwutils.rscaleConstM(mass/nfwutils.global_cosmology.h, concen, zlens, 200)

    

    gamma = nfwmodeltools.NFWShear(cat['r_mpc'], concen, rscale, zlens)
    kappa = nfwmodeltools.NFWKappa(cat['r_mpc'], concen, rscale, zlens)

    gpred = cat['beta_s']*gamma / (1 - (cat['beta_s2']*kappa/cat['beta_s']))


    fig = pylab.figure()
#    pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat']/gpred, cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])*gpred), fmt='bo')
#    pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat']/gpred, cat['ghatdistrosigma']/(gpred), fmt='bo')
    pylab.errorbar(cat['r_mpc'], cat['ghat'], cat['ghatdistrosigma']/np.sqrt(cat['ndat']), fmt='bo')

#
    ax = pylab.gca()
    ax.set_xscale('log')
    pylab.axhline(0.0, c='k', linewidth=2)
    pylab.xlabel('Radius [Mpc/h]', fontsize=16)
    pylab.ylabel('<g_meas/g_pred> - 1)', fontsize=16)
#    pylab.ylabel('<g_meas - g_pred>', fontsize=16)

    pylab.axis([0.05, 10, -.55, 0.35])


#    ax2 = ax.twinx()
#    pylab.plot(cat['r_mpc'], gpred, 'k--')
#    ax2.errorbar(cat['r_mpc'], cat['ghat'], cat['ghatdistrosigma']/np.sqrt(cat['ndat']), fmt='rs')
#    ax2.set_ylabel('<g_meas - g_pred>', color='r', fontsize=16)
#    ax2.set_ylim(-0.2, 0.1)
#
#    pylab.title('Redshift=%1.2f Mass=%1.2fx10^14 Concen=%1.2f' % (zlens, mass/1e14, concen))
#
    pylab.tight_layout()

    pylab.savefig('%s.png' % binbase)

    return fig

############

c = 'SlateGray r b m DodgerBlue g DarkSalmon'.split()

############

def plotProfile(datfilebase, c, label, labelZ = True):

    mass, concen, redshift =  readtxtfile.readtxtfile('%s.dat' % datfilebase)[0]

    cat = ldac.openObjectFile('%s.cat' % datfilebase)
    
    zlens = cat.hdu.header['ZLENS']
    rscale = nfwutils.rscaleConstM(mass/nfwutils.global_cosmology.h, concen, zlens, 200)
    gamma = nfwmodeltools.NFWShear(cat['r_mpc'], concen, rscale, zlens)
    kappa = nfwmodeltools.NFWKappa(cat['r_mpc'], concen, rscale, zlens)
    gpred = cat['beta_s']*gamma / (1 - (cat['beta_s2']*kappa/cat['beta_s']))



    signalpred = gpred / (cat['beta_s']*nfwutils.global_cosmology.beta([1e6], zlens)*nfwutils.global_cosmology.angulardist(zlens))
            
    if labelZ:
        flabel='%s z=%1.2f' % (label, zlens)
    else:
        flabel =label

    pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat'], cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])), 
    linestyle='None', marker='o', color=c, label=flabel)
    pylab.plot(cat['r_mpc']*nfwutils.global_cosmology.h, signalpred, marker='None', linestyle=':', color=c, linewidth=2)
    
####

def multibinshear():

    matplotlib.rcParams['figure.figsize'] = [32,8]

    fig = pylab.figure()

    curplot = 1
    for curm in range(4):
        
        pylab.subplot(1,4,curplot)

        colori = -1

        #first plot bk11

        setBK11()

        for snap in [124,141]:

            colori += 1

 #           try:
            plotProfile('bk11stack_massbins_shearprofile_%d/bk11stack_%d' % (snap, curm),
                        c[colori], 'BK11')
#            except:
#                pass


        # then bcc

        setBCC()

        for curz in range(3):

            colori += 1

 #           try:
            plotProfile('bccstack_massbins_shearprofile/bccstack_%d_%d' % (curz, curm),
                        c[colori], 'BCC')
#            except:
#                pass


        # then mxxl

        setMXXL()

        if curm == 3:

            for snap in [41]:

                for mxxlmass, mxxllabel in enumerate('Low High'.split()):

                    colori += 1

#                    try:
                    plotProfile('mxxlstack_massbins_shearprofile_%d/mxxlstack_%d' % (snap, mxxlmass),
                                    c[colori], 'MXXL-%s' % mxxllabel)
#                    except:
#                        pass



        ax = pylab.gca()
        ax.set_xscale('log')
        pylab.axhline(0.0, c='k', linewidth=2)
        pylab.axvline(0.2*0.7, c='g', linestyle='--')
        pylab.axis([0.01, 10, 0., 0.0003])
        pylab.xlabel('Radius [Mpc/h]')
        

        curplot += 1

    pylab.subplot(1,4,1)
    pylab.ylabel('shear')

    pylab.subplot(1,4,4)
    pylab.legend(loc='upper right', fontsize=10)
            
        
    pylab.tight_layout()


    pylab.savefig('multibin_shear.png')

    return fig

######################################################


def multibinresidual(binbase, fig = None):

    matplotlib.rcParams['figure.figsize'] = [16,16]

    if fig is None:
        fig = pylab.figure()
    curplot = 1
    for curm in range(4):
        for curc in range(4):
            pylab.subplot(4,4,curplot)
            colori = 0
#            for curz in range(3):
            for curz in [0]:
                
                mass, concen, redshift =  readtxtfile.readtxtfile('%s_%d_%d_%d.dat' % (binbase, curz, curm, curc))[0]

                cat = ldac.openObjectFile('%s_%d_%d_%d.cat' % (binbase, curz, curm, curc))

                zlens = cat.hdu.header['ZLENS']

                rscale = nfwutils.rscaleConstM(mass/nfwutils.global_cosmology.h, concen, zlens, 200)



                gamma = nfwmodeltools.NFWShear(cat['r_mpc'], concen, rscale, zlens)
                kappa = nfwmodeltools.NFWKappa(cat['r_mpc'], concen, rscale, zlens)

                gpred = cat['beta_s']*gamma / (1 - (cat['beta_s2']*kappa/cat['beta_s']))




                pylab.errorbar(cat['r_mpc']*0.72, cat['ghat'], cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])), 
                               linestyle='None', marker='o', color=c[colori], label='BCC z=%1.1f' % redshift)
                pylab.plot(cat['r_mpc']*0.72, gpred, 'k-', linewidth=2)

                ax = pylab.gca()
                ax.set_xscale('log')
                pylab.axhline(0.0, c='k', linewidth=2)
                pylab.legend(loc='lower center', fontsize=10)



                pylab.axis([0.05, 10, -.55, 0.35])
#                pylab.axis([0.05, 10, -.10, 0.05])



                colori+= 1
            pylab.title('M=%1.1fx10^14 C=%1.1f' % ( mass/1e14, concen))
            curplot += 1

            
    for i in range(4):
        pylab.subplot(4,4,13+i)
        pylab.xlabel('Radius [Mpc/h]')
        pylab.subplot(4,4,4*i+1)
#        pylab.ylabel('<g_m-g_p>/g_p(<M>,<c>)')
#        pylab.ylabel('<g_m/g_p - 1>')



    pylab.tight_layout()

#    pylab.savefig('%s_multibin_fracresid.png' % binbase)
    pylab.savefig('%s_multibin_resid.png' % binbase)

    return fig

######################################################

def BK11multibinresidual(dirbase):

    setBK11()

    matplotlib.rcParams['figure.figsize'] = [16,16]


    fig = pylab.figure()
    curplot = 1
    for curm in range(4):
        for curc in range(4):
            pylab.subplot(4,4,curplot)
            colori = 1

            for snap in [124,141]:
        
                try:

                    mass, concen, redshift =  readtxtfile.readtxtfile('%s_%d/bk11stack_%d_%d.dat' % (dirbase, snap, 
                                                                                                     curm, curc))[0]

                    cat = ldac.openObjectFile('%s_%d/bk11stack_%d_%d.cat' % (dirbase, snap, 
                                                                             curm, curc))


                    zlens = cat.hdu.header['ZLENS']

                    rscale = nfwutils.rscaleConstM(mass/nfwutils.global_cosmology.h, concen, zlens, 200)



                    gamma = nfwmodeltools.NFWShear(cat['r_mpc'], concen, rscale, zlens)
                    kappa = nfwmodeltools.NFWKappa(cat['r_mpc'], concen, rscale, zlens)
                    
                    gpred = cat['beta_s']*gamma / (1 - (cat['beta_s2']*kappa/cat['beta_s']))

                    pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat'], cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])), 
                                   linestyle='None', marker='o', color=c[colori], label='M=%1.1fx10^14' % (mass/1e14))
                    
                    pylab.plot(cat['r_mpc']*nfwutils.global_cosmology.h, gpred, 'k-', linewidth=2)
                    
                    ax = pylab.gca()
                    ax.set_xscale('log')
                    pylab.axhline(0.0, c='k', linewidth=2)
                    
                    pylab.axis([0.05, 10, -0.03, 0.4])


    




                    ax = pylab.gca()
                    ax.set_xscale('log')
                    pylab.axhline(0.0, c='k', linewidth=2)
#                    pylab.legend(loc='lower center', fontsize=10)

#                    pylab.title('M=%1.1fx10^14 C=%1.1f' % ( mass/1e14, concen))

                except:
                    pass

                colori+= 1

            curplot += 1

            
    for i in range(4):
        pylab.subplot(4,4,13+i)
        pylab.xlabel('Radius [Mpc/h]')
        pylab.subplot(4,4,4*i+1)
        pylab.ylabel('<g_m/g_p - 1>')


    pylab.tight_layout()

    pylab.savefig('%s_multibin_shearprofile.png' % dirbase)
#    pylab.savefig('%s_multibin_resid.png' % binbase)

    return fig


#######################################################

def multibinresidualoverplot(binbase, fig = None):


    setBCC()

    if fig is None:

        matplotlib.rcParams['figure.figsize'] = [16,8]
        fig = pylab.figure()


    curplot = 1
    for curc in range(4):
        pylab.subplot(2,4,curplot)
        for curm in [3]:

            colori = 0
            for curz in range(3):
                
                mass, concen, redshift =  readtxtfile.readtxtfile('%s_%d_%d_%d.dat' % (binbase, curz, curm, curc))[0]

                cat = ldac.openObjectFile('%s_%d_%d_%d.cat' % (binbase, curz, curm, curc))


                zlens = cat.hdu.header['ZLENS']

                rscale = nfwutils.rscaleConstM(mass/nfwutils.global_cosmology.h, concen, zlens, 200)



                gamma = nfwmodeltools.NFWShear(cat['r_mpc'], concen, rscale, zlens)
                kappa = nfwmodeltools.NFWKappa(cat['r_mpc'], concen, rscale, zlens)

                gpred = cat['beta_s']*gamma / (1 - (cat['beta_s2']*kappa/cat['beta_s']))

                pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat'], cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])), 
                                   linestyle='None', marker='o', color=c[colori], label='M=%1.1fx10^14' % (mass/1e14))

                pylab.plot(cat['r_mpc']*nfwutils.global_cosmology.h, gpred, 'k-', linewidth=2)

                ax = pylab.gca()
                ax.set_xscale('log')
                pylab.axhline(0.0, c='k', linewidth=2)

                pylab.axis([0.05, 10, -0.03, 0.4])
                zlens = cat.hdu.header['ZLENS']

                rscale = nfwutils.rscaleConstM(mass/nfwutils.global_cosmology.h, concen, zlens, 200)



                gamma = nfwmodeltools.NFWShear(cat['r_mpc'], concen, rscale, zlens)
                kappa = nfwmodeltools.NFWKappa(cat['r_mpc'], concen, rscale, zlens)

                gpred = cat['beta_s']*gamma / (1 - (cat['beta_s2']*kappa/cat['beta_s']))



#                pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat']/gpred, cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])*gpred), 
#                               linestyle='None', marker='o', color=c[colori], label='z=%1.1f' % redshift)
                pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat'], cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])), 
                               linestyle='None', marker='o', color=c[colori], label='z=%1.1f' % redshift)

                pylab.plot(cat['r_mpc']*nfwutils.global_cosmology.h, gpred, 'k-', linewidth=2)

            #    pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat']/gpred, cat['ghatdistrosigma']/(gpred), fmt='bo')

                pylab.plot(cat['r_mpc']*nfwutils.global_cosmology.h, gpred, 'r-')

            #
                ax = pylab.gca()
                ax.set_xscale('log')
                pylab.axhline(0.0, c='k', linewidth=2)



                pylab.axis([0.05, 10, -0.03, 0.4])



                colori+= 1
        pylab.title('C=%1.1f' % (concen))
        curplot += 1

            
    for i in range(4):
        pylab.subplot(2,4,i+1)
        pylab.xlabel('Radius [Mpc/h]')
        pylab.text(0.25, 0.2, 'BCC Mass=%1.1fx10^14' % (mass/1e14))
    pylab.subplot(2,4,1)
    pylab.ylabel('<g_m/g_p-1>')
    pylab.subplot(2,4,4)
    pylab.legend(loc='upper right')


    pylab.tight_layout()

#    pylab.savefig('%s_multibin_resid_overplot.png' % binbase)


    return fig




#######################################################

def setMXXL():

    nfwutils.global_cosmology.set_cosmology(nfwutils.Cosmology(omega_m = 0.25, omega_l = 0.75, h=0.73))

def setBCC():

    nfwutils.global_cosmology.set_cosmology(nfwutils.Cosmology(omega_m = 0.23, omega_l = 0.77, h=.72))

def setBK11():

    nfwutils.global_cosmology.set_cosmology(nfwutils.Cosmology(omega_m = 0.27, omega_l = 0.73, h = 0.7))


#############################################################


def MXXLmultibinresidualoverplot(binbase, fig = None):

    setMXXL()

    if fig is None:
        matplotlib.rcParams['figure.figsize'] = [16,8]
        fig = pylab.figure()


    curplot = 5
    for curc in range(4):
        pylab.subplot(2,4,curplot)

        colori = 0
        for curm in range(2):

            try:
                
                mass, concen =  readtxtfile.readtxtfile('%s_%d_%d.dat' % (binbase, curm, curc))[0]

                cat = ldac.openObjectFile('%s_%d_%d.cat' % (binbase, curm, curc))

                zlens = cat.hdu.header['ZLENS']

                rscale = nfwutils.rscaleConstM(mass/nfwutils.global_cosmology.h, concen, zlens, 200)



                gamma = nfwmodeltools.NFWShear(cat['r_mpc'], concen, rscale, zlens)
                kappa = nfwmodeltools.NFWKappa(cat['r_mpc'], concen, rscale, zlens)

                gpred = cat['beta_s']*gamma / (1 - (cat['beta_s2']*kappa/cat['beta_s']))

                pylab.errorbar(cat['r_mpc']*nfwutils.global_cosmology.h, cat['ghat'], cat['ghatdistrosigma']/(np.sqrt(cat['ndat'])), 
                                   linestyle='None', marker='o', color=c[colori], label='M=%1.1fx10^14' % (mass/1e14))

                pylab.plot(cat['r_mpc']*nfwutils.global_cosmology.h, gpred, 'k-', linewidth=2)

                ax = pylab.gca()
                ax.set_xscale('log')
                pylab.axhline(0.0, c='k', linewidth=2)

                pylab.axis([0.05, 10, -0.03, 0.4])


            except:
                pass
    


            colori+= 1

        pylab.title('C=%1.1f' % (concen))
        curplot += 1

            
    for i in range(4):
        pylab.subplot(2,4,4+i+1)
        pylab.xlabel('Radius [Mpc/h]')
        pylab.text(0.1, 0.03, 'MXXL')
        pylab.minorticks_on()
    pylab.subplot(2,4,5)
    pylab.ylabel('<g_m/g_p-1>')
    pylab.subplot(2,4,8)
    pylab.legend(loc='upper right')


    pylab.tight_layout()

    pylab.savefig('%s_shearstack.png' % binbase)


    return fig




#######################################################

    
def run(stackdir):

#    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    stackfiles = glob.glob('%s/*.cat' % stackdir)
    stackbases = [x[:-4] for x in stackfiles]
    figs = []
    for stackbase in stackbases:
#        try:
        figs.append(residual(stackbase))
#        except ValueError:
#            pass
                       

    return figs

    
    

    
    

    

    



