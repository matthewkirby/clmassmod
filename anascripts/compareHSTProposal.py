############
# This script compares bias and scatter calculations for different HST proposal set ups
#  with source densities and betas approximately scaled for different cluster redshifts
#############

import os
import bcc_becker_baseline as bbb
import nfwutils


mxxldir='/users/dapple/astro/mxxlsims/mxxl_imperial/snap41'

mcrs = 'c4 cfree'.split()
redshifts = range(0, 5, 2)
strategies = 'mosaic deepmosaic cvzsquare'.split()

mxxldata = {}

for mcr in mcrs:

    if mcr not in mxxldata:
        mxxldata[mcr] = {}

    for redshift in redshifts:

        if redshift not in mxxldata[mcr]:
            mxxldata[mcr][redshift] = {}

        for strategy in strategies:

            simname = '{0}-hst_{1}_{2}'.format(mcr, strategy, redshift)
            mxxlfilename='{0}/{1}/consolidated.pkl'.format(mxxldir, simname)

            mxxldata[mcr][redshift][strategy] = bbb.loadData(mxxlfilename)

   


c = 'r b m DodgerBlue g DarkSalmon'.split()

linestyle='- -- :'.split()
#
redshift_trans = {
    0 : 'z = 0.7',
    1 : 'z = 0.8',
    2 : 'z = 0.9',
    3 : 'z = 1.0',
    4 : 'z = 1.1'
}



figsize(8,6)

#####
## Bias
#
for mcr in mcrs:

    figure()

    linecount = 0

    
    for strategy in strategies:

        ccount = 0        

        for redshift in redshifts:

            axhline(1.0, c='DimGray', linestyle='--')

            simdat = mxxldata[mcr][redshift][strategy]
            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
            plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], 
                 label='{0} {1}'.format(strategy, redshift_trans[redshift]))

            ccount += 1

        linecount += 1



    xlabel(r'Mass $M_{200}$', fontsize=16)
    ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=16)
    axis([14.1, 15.2, 0.7, 1.3])
    legend(loc='upper left')
    title('MC Relation: {0}'.format(mcr))
    tight_layout()
    savefig('hstproposal_{0}.png'.format(mcr))
    
####
## Scatter

for mcr in mcrs:
        
    figure()

    linecount = 0

    
    for strategy in strategies:

        ccount = 0        

        for redshift in redshifts:

            axhline(1.0, c='DimGray', linestyle='--')

            simdat = mxxldata[mcr][redshift][strategy]
            
            massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
            plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount],
                 label='{0} {1}'.format(strategy, redshift_trans[redshift]))
            
            ccount += 1
          
        linecount += 1

    
    xlabel(r'Mass $M_{200}$', fontsize=16)
    ylabel(r'Distribution Scatter', fontsize=16)
    axis([14.1, 15., 0.0, 0.6])
    legend(loc='upper left')
    title('Scatter - MC Relation: {0}'.format(mcr))
    tight_layout()
    savefig('hstproposal_{0}_scatter.png'.format(mcr))

