############
# This script compares bias and scatter calculations for different source densities
#############

import os
import bcc_becker_baseline as bbb
import nfwutils


mxxldir='/users/dapple/astro/mxxlsims/mxxl_imperial/snap41'
bccdir='/users/dapple/braid1/bcc'

mcrs = 'c4 cfree'.split()
radialranges = ['r{0}'.format(i) for i in '1 5 10'.split()]
#samplings = ['n0_0'] + ['n{0}_2'.format(i) for i in range(6)]
samplings = ['n{0}_0'.format(i) for i in range(6)]


#mxxldata = {}
#bccdata = {}
##
#for mcr in mcrs:
#
#    if mcr not in mxxldata:
#        mxxldata[mcr] = {}
#    if mcr not in bccdata:
#        bccdata[mcr] = {}
#
#    for radrange in radialranges:
#
#        if radrange not in mxxldata[mcr]:
#            mxxldata[mcr][radrange] = {}
#        if radrange not in bccdata[mcr]:
#            bccdata[mcr][radrange] = {}
#
#        for sample in samplings:
#
#            simname = '{0}-{1}-{2}'.format(mcr, radrange, sample)
#            mxxlfilename='{0}/{1}/consolidated.pkl'.format(mxxldir, simname)
#
#            mxxldata[mcr][radrange][sample] = bbb.loadData(mxxlfilename)
#
#            bccfilename='{0}/{1}/consolidated.pkl'.format(bccdir, simname)
#            if not os.path.exists(bccfilename):
#                print 'Skipping BCC {0}'.format(simname)
#                continue
#            bccdata[mcr][radrange][sample] = bbb.loadData(bccfilename)
#   
#


c = 'SlateGray r b m DodgerBlue g DarkSalmon'.split()

linestyle='-- - :'.split()
#
radrange_trans = {
'r1' : '0.25 - 0.5 Mpc',
'r2' : '0.25 - 1.5 Mpc',
'r3' : '0.25 - 2.5 Mpc',
'r4' : '0.25 - 3.0 Mpc',
'r5' : '0.50 - 1.5 Mpc',
'r6' : '0.50 - 2.5 Mpc',
'r7' : '0.50 - 3.0 Mpc',
'r8' : '0.75 - 1.5 Mpc',
'r9' : '0.75 - 2.5 Mpc',
'r10' :'0.75 - 3.0 Mpc',
'r11' :'0.25 - 0.75 Mpc',
'r12' :'0.25 - 1.0 Mpc',
'r13' :'0.50 - 0.75 Mpc',
'r14' :'0.50 - 1.0 Mpc',
'r15' :'0.75 - 1.0 Mpc'
}
#
#
sampling_trans = {
'n0_0' : '-1',
'n1_0' : '100',
'n2_0' : '20',
'n3_0' : '7',
'n4_0' : '4',
'n5_0' : '1',
'n0_1' : '0.16',
'n0_2' : '0.33',
'n0_3' : '0.50',
'n1_0' : '100 0.00',
'n1_1' : '100 0.16',
'n1_2' : '100 0.33',
'n1_3' : '100 0.50',
'n2_0' : '20 0.00',
'n2_1' : '20 0.16',
'n2_2' : '20 0.33',
'n2_3' : '20 0.50',

'n3_0' : '7 0.00',
'n3_1' : '7 0.16',
'n3_2' : '7 0.33',
'n3_3' : '7 0.50',

'n4_2' : '4 0.33',
'n5_2' : '1 0.33',

}
#
#
#figsize(8,8)
#
#####
## Bias
#
#for mcr in mcrs:
#
#    for radrange in radialranges:
#
#        figure()
#
#
#        ccount = 0        
#        for curdens, sampling in enumerate(samplings):
#
#            axhline(1.0, c='DimGray', linestyle='--')
#
##            linecount = 0
##            if sampling in bccdata[mcr][radrange]:
##                simdat = bccdata[mcr][radrange][sampling]
##                massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
##                plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], 
##                     label='BCC {0}'.format(sampling_trans[sampling]))
##
#
#
#            linecount = 1
#            simdat = mxxldata[mcr][radrange][sampling]
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
#            plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], 
#                 label='MXXL {0}'.format(sampling_trans[sampling]))
#
#            ccount += 1
#
#
#
#
#
#
#
#
#        xlabel(r'Mass $M_{200}$', fontsize=16)
#        ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=16)
#        axis([14.1, 15.2, 0.7, 1.3])
#        legend(loc='upper left')
#        title('MC Relation: {0} RadialRange: {1}'.format(mcr, radrange_trans[radrange]))
#        tight_layout()
#        savefig('density_noise2_{0}-{1}.png'.format(mcr, radrange))
#
####
## Scatter


#for mcr in mcrs:
#
#    for radrange in radialranges:
#
#        figure()
#
#
#        ccount = 0        
#        for curdens, sampling in enumerate(samplings):
#
#            axhline(1.0, c='DimGray', linestyle='--')
#
#            linecount = 0
#            if sampling in bccdata[mcr][radrange]:
#                simdat = bccdata[mcr][radrange][sampling]
#                massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
#                plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount],
#                     label='BCC {0}'.format(sampling_trans[sampling]))
#
#
#
#            linecount = 1
#            simdat = mxxldata[mcr][radrange][sampling]
#            massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
#            plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount],
#                 label='MXXL {0}'.format(sampling_trans[sampling]))
#
#            ccount += 1
#
#
#    
#        xlabel(r'Mass $M_{200}$', fontsize=16)
#        ylabel(r'Distribution Scatter', fontsize=16)
#        axis([14.1, 15., 0.0, 0.6])
#        legend(loc='upper left')
#        title('Scatter - MC Relation: {0} Radial Range: {1}'.format(mcr, radrange_trans[radrange]))
#        tight_layout()
#        savefig('density_{0}-{1}_scatter.png'.format(mcr, radrange))
#
#

################
# Investigation of BCC Scatter
#
#
#
#
simdat = bccdata['cfree']['r10']['n0_0']
#
##bias
#figure()
#ccount = 0
#linecount=1
#zbins = [(0., 0.3), (0.35, 0.5), (0.5, 0.89), (0.9, 1.5)]
#
#for zbin in zbins:
#    
#    massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat, np.logical_and(simdat['redshifts'] >= zbin[0],
#                                                                           simdat['redshifts'] < zbin[1]))
#    plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], 
#         label='{0} {1}'.format(*zbin))
#
#    ccount += 1
#
#
#xlabel(r'Mass $M_{200}$', fontsize=16)
#ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=16)
#axis([14.1, 15.2, 0.7, 1.3])
#legend(loc='upper left')
#title('BCC Bias Redshift Evolution')
#tight_layout()
#savefig('bcc_bias_redshift_evo_r10_n0_0.png')
#
#

####scatter
#
figure()
ccount = 0
linecount=1
zbins = [(0., 0.3), (0.35, 0.5), (0.5, 0.89), (0.9, 1.5)]

for zbin in zbins:

    massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat, np.logical_and(simdat['redshifts'] >= zbin[0],simdat['redshifts'] < zbin[1]))
    plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount],
         label='{0} {1}'.format(*zbin))

    ccount += 1


xlabel(r'Mass $M_{200}$', fontsize=16)
ylabel(r'Distribution Scatter', fontsize=16)
axis([14.1, 15., 0.0, 0.6])
legend(loc='upper left')
title('BCC Scatter Redshift Evolution')
tight_layout()
savefig('bcc_scatter_redshift_evo_r10_n0_0.png')




############
# Evolution with noise at fixed sampling

noisesamplings = ['n2_{0}'.format(i) for i in range(4)]

###load data
#for mcr in mcrs:
#
#    if mcr not in mxxldata:
#        mxxldata[mcr] = {}
#    if mcr not in bccdata:
#        bccdata[mcr] = {}
#
#    for radrange in radialranges:
#
#        if radrange not in mxxldata[mcr]:
#            mxxldata[mcr][radrange] = {}
#        if radrange not in bccdata[mcr]:
#            bccdata[mcr][radrange] = {}
#
#        for sample in noisesamplings:
#
#            simname = '{0}-{1}-{2}'.format(mcr, radrange, sample)
#            mxxlfilename='{0}/{1}/consolidated.pkl'.format(mxxldir, simname)
#
#            mxxldata[mcr][radrange][sample] = bbb.loadData(mxxlfilename)
#
#            bccfilename='{0}/{1}/consolidated.pkl'.format(bccdir, simname)
#            if not os.path.exists(bccfilename):
#                print 'Skipping BCC {0}'.format(simname)
#                continue
#            bccdata[mcr][radrange][sample] = bbb.loadData(bccfilename)
#
#
#
figsize(8,8)

####
# Bias

#for mcr in mcrs:
#
#    for radrange in radialranges:
#
#        figure()
#
#
#        ccount = 0        
#        for curdens, sampling in enumerate(noisesamplings):
#
#            axhline(1.0, c='DimGray', linestyle='--')
#
#            linecount = 0
#            if sampling in bccdata[mcr][radrange]:
#                simdat = bccdata[mcr][radrange][sampling]
#                massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
#                plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], 
#                     label='BCC {0}'.format(sampling_trans[sampling]))
#
#
#
#            linecount = 1
#            simdat = mxxldata[mcr][radrange][sampling]
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
#            plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], 
#                 label='MXXL {0}'.format(sampling_trans[sampling]))
#
#            ccount += 1
#
#
#
#
#
#
#
#
#        xlabel(r'Mass $M_{200}$', fontsize=16)
#        ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=16)
#        axis([14.1, 15.2, 0.7, 1.3])
#        legend(loc='upper left')
#        title('MC Relation: {0} RadialRange: {1}'.format(mcr, radrange_trans[radrange]))
#        tight_layout()
##        savefig('noise_{0}-{1}.png'.format(mcr, radrange))
#
####
## Scatter
#
#
#for mcr in mcrs:
#
#    for radrange in radialranges:
#
#        figure()
#
#
#        ccount = 0        
#        for curdens, sampling in enumerate(noisesamplings):
#
#            axhline(1.0, c='DimGray', linestyle='--')
#
#            linecount = 0
#            if sampling in bccdata[mcr][radrange]:
#                simdat = bccdata[mcr][radrange][sampling]
#                massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
#                plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount],
#                     label='BCC {0}'.format(sampling_trans[sampling]))
#
#
#
#            linecount = 1
#            simdat = mxxldata[mcr][radrange][sampling]
#            massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
#            plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount],
#                 label='BCC {0}'.format(sampling_trans[sampling]))
#
#            ccount += 1
#
#
#    
#        xlabel(r'Mass $M_{200}$', fontsize=16)
#        ylabel(r'Distribution Scatter', fontsize=16)
#        axis([14.1, 15., 0.0, 0.6])
#        legend(loc='upper left')
#        title('Scatter - MC Relation: {0} Radial Range: {1}'.format(mcr, radrange_trans[radrange]))
#        tight_layout()
#        savefig('noise_{0}-{1}_scatter.png'.format(mcr, radrange))



########################


