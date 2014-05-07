#/usr/bin/env python
#######################
# Take a list of catalogs, and compile mean profiles, standard deviations, and any other appropriate info
#######################

import os, sys, cPickle
import numpy as np
import nfwfit, readtxtfile, nfwutils
import astropy.io.fits as pyfits, ldac
import consolidate_fits

###################

def calcOnlineStats(nA, meanA, mom2A, nB, meanB, mom2B):

    nX = nA + nB

    delta = meanB - meanA
    meanX = meanA + delta*nB/nX
    mom2X = mom2A + mom2B + (delta**2)*nA*nB/nX

    return meanX, mom2X
    

###################


class OnlineStatistics(object):

    def __init__(self, binedges):

        self.ncats = 0.
        self.binedges = binedges
        self.nbins = len(binedges) - 1
        self.ndata = np.zeros(self.nbins)
        self.meanradii = np.zeros(self.nbins)
        self.meanshear = np.zeros(self.nbins)
        self.mom2shear = np.zeros(self.nbins)
        self.profileCol = 'r_mpc'
        self.shearCol = 'ghat'

    #########

    @property
    def varianceshear(self):

        return self.mom2shear/(self.ndata - 1)

    #########


    def accumulate(self, catalog, truth):

        zlens = catalog.hdu.header['ZLENS']
        m200 = truth['m200']
        concen = truth['concen']
        rscale = nfwutils.rscaleConstM(m200/nfwutils.global_cosmology.h, concen, zlens, 200)

        gamma = nfwmodeltools.NFWShear(self.profileCol, concen, rscale, zlens)
        kappa = nfwmodeltools.NFWKappa(self.profileCol, concen, rscale, zlens)

        gpred = cat['beta_s']*gamma / (1 - (cat['beta_s']*kappa))

                                       
        g_resid = cat[self.shearCol] - gpred

        
        for curbin in range(self.nbins):
            mintake = self.binedges[curbin]
            maxtake = self.binedges[curbin+1]
            inbin = np.logical_and(catalog[self.profileCol] >= mintake,
                                                     catalog[self.profileCol] < maxtake)
            selected = catalog.filter(inbin)

            catndata = len(selected)

            catmeanradius = np.mean(selected[self.profileCol])
            catmeanshear = np.mean(g_resid[selected])
            catvarianceshear = np.std(g_resid[selected])**2
            

            self.meanradii[curbin], junk = calcOnlineStats(self.ndata[curbin], self.meanradii[curbin], 0.,
                                                       catndata, catmeanradius, 0.)
            
            self.meanshear[curbin], self.mom2shear[curbin] = calcOnlineStats(self.ndata[curbin], 
                                                                             self.meanshear[curbin],
                                                                             self.mom2shear[curbin],
                                                                             catndata, 
                                                                             catmeanshear, 
                                                                             catvarianceshear)


            

            self.ndata[curbin] += catndata

        self.ncats += 1


    ########

    def writeSimCat(self, outfile):

        cols = [pyfits.Column(name = 'r_mpc', format='E', array = self.meanradii),
                pyfits.Column(name = 'ghat', format='E', array = self.meanshear),
                pyfits.Column(name = 'ghatdistrosigma', format='E', array = np.sqrt(self.varianceshear)),
                pyfits.Column(name = 'ndat', format='E', array = self.ndata)]
        catalog = ldac.LDACCat(pyfits.new_table(pyfits.ColDefs(cols)))
        catalog.hdu.header.update('ZLENS', self.meanzlens)

        catalog.saveas(outfile, clobber=True)



############################


def stackCats(stackfile, configname, answerfile, outfile):

    with open(answerfile, 'rb') as input:
        answers = cPickle.load(input)

    tostack = [x[0] for x in readtxtfile.readtxtfile(stackfile)]

    config = nfwfit.readConfiguration(configname)

    simreader = nfwfit.buildSimReader(config)

    nfwutils.global_cosmology.set_cosmology(simreader.getCosmology())

    fitter = nfwfit.buildFitter(config)

    profile = fitter.profileBuilder

    if profile.binspacing == 'linear':
        binedges = np.linspace(profile.minradii, profile.maxradii, profile.nbins+1)
    else:
        binedges = np.logspace(np.log10(profile.minradii), np.log10(profile.maxradii), profile.nbins+1)


    stackedprofile = OnlineStatistics(binedges)

    for catalogname in tostack:

        filebase = os.path.basename(output)
        
        match = consolidate_fits.idpattern.match(filebase)
        
        haloid = int(match.group(1))

        try:
            truth = answers[haloid]
        except KeyError:
            print 'Failure at {0}'.format(output)
            raise


        catalog = nfwfit.readSimCatalog(catalogname, simreader, config)

        stackedprofile.accumulate(catalog, truth)


    stackedprofile.writeSimCat(outfile)
 

    return stackedprofile


############################

def assignMXXLStacks(outdir, massedges = np.array([0, 3.8e14, 4.2e14, 4.9e14, 5e15]),
              concenedges = np.array([0, 0.2, 0.26, 0.38, 4.38, 10])):

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with open('mxxl_answers.pkl', 'rb') as input:
        answers = cPickle.load(input)

    nclusters = len(answers)

    masses = np.zeros(len(answers))
    concens = np.zeros(len(answers))

    halos = np.arange(nclusters)

    for i in range(nclusters):

        masses[i] = answers[i]['m200']
        concens[i] = answers[i]['concen']

    haloassignments = {}

    
    condorfile = open('%s/mxxlstack.submit' % outdir, 'w')
    condorfile.write('executable = /vol/braid1/vol1/dapple/mxxl/mxxlsims/stack_condorwrapper.sh\n')
    condorfile.write('universe = vanilla\n')
    
    for curmass_i in range(len(massedges)-1):
        for curconcen_i in range(len(concenedges)-1):

            massselect = np.logical_and(masses >= massedges[curmass_i],
                                        masses < massedges[curmass_i+1])
            concenselect = np.logical_and(concens >= concenedges[curconcen_i],
                                          concens < concenedges[curconcen_i + 1])
            binselect = np.logical_and(massselect, concenselect)
            selected = halos[binselect]

            with open('%s/mxxlstack_%d_%d.list' % (outdir, curmass_i, curconcen_i), 'w') as output:

                for curhalo in selected:
                    output.write('/vol/braid1/vol1/dapple/mxxl/snap41/halo_cid%d\n' % curhalo)

            with open('%s/mxxlstack_%d_%d.dat' % (outdir, curmass_i, curconcen_i), 'w') as output:

                output.write('# <Mass> <Concen>\n')
                output.write('%f %f\n' % (np.mean(masses[binselect]), np.mean(concens[binselect])))

            condorfile.write('Error = %s/consolidate.mxxlstack_%d_%d.stderr\n' % (outdir, curmass_i, curconcen_i))
            condorfile.write('Output = %s/consolidate.mxxlstack_%d_%d.stdout\n' % (outdir, curmass_i, curconcen_i))
            condorfile.write('Log = %s/consolidate.mxxlstack_%d_%d.batch.log\n' % (outdir, curmass_i, curconcen_i))
            condorfile.write('Arguments = %s/mxxlstack_%d_%d.list /vol/braid1/vol1/dapple/mxxl/mxxlsims/stackconfig.sh %s/mxxlstack_%d_%d.cat\n' % (outdir, curmass_i, curconcen_i,
                                                                                                                                                  outdir, curmass_i, curconcen_i))
            condorfile.write('queue\n')



    condorfile.close()
                                                                                                                                                  
            

#############################

def assignBCCStacks(outdir, redshiftedges = np.array([0., 0.4, 0.6, 1.0]),
                    massedges = np.array([0., 2.2e14, 2.6e14, 3.2e14, 1e16]),
                    concenedges = np.array([0, 2.9, 3.7, 4.7, 10])):

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with open('bcc_answers.pkl', 'rb') as input:
        answers = cPickle.load(input)

    nclusters = len(answers)

    masses = np.zeros(len(answers))
    concens = np.zeros(len(answers))
    redshifts = np.zeros(len(answers))

    halos = np.array([x for x in answers.keys()])

    for i in range(nclusters):

        haloid = halos[i]

        masses[i] = answers[haloid]['m200']
        concens[i] = 1000*answers[haloid]['concen']
        redshifts[i] = answers[haloid]['redshift']

    haloassignments = {}

    for curz_i in range(len(redshiftedges) - 1):

        redshiftselect = np.logical_and(redshifts >= redshiftedges[curz_i], redshifts < redshiftedges[curz_i+1])

        for curmass_i in range(len(massedges) - 1):

            massselect = np.logical_and(masses >= massedges[curmass_i], masses < massedges[curmass_i +1])

            for curconcen_i in range(len(concenedges) - 1):

                concenselect = np.logical_and(concens >= concenedges[curconcen_i], concens < concenedges[curconcen_i+1])

                inbin = np.logical_and(np.logical_and(redshiftselect, massselect), concenselect)
                


                with open('%s/bccstack_%d_%s_%s.list' % (outdir, curz_i, curmass_i, curconcen_i), 'w') as output:

                    for curhalo in halos[inbin]:
                        output.write('/u/ki/dapple/nfs/bcc_clusters/recentered/cluster_%d.hdf5\n' % curhalo)

                with open('%s/bccstack_%d_%s_%s.dat' % (outdir, curz_i, curmass_i, curconcen_i), 'w') as output:

                    output.write('# <Mass> <Concen> <redshift>\n')
                    output.write('%f %f %f\n' % (np.mean(masses[inbin]), np.mean(concens[inbin]), np.mean(redshifts[inbin])))




                                                                                                                                                  
            
                             





##############################

if __name__ == '__main__':

    stackfile = sys.argv[1]
    config = sys.argv[2]
    answerkey = sys.argv[3]
    outfile = sys.argv[4]

    stackCats(stackfile, config, answerkey, outfile)
