############
# This script compares bias and scatter calculations for different HST proposal set ups
#  with source densities and betas approximately scaled for different cluster redshifts
#############

import os
import bcc_becker_baseline as bbb
import nfwutils


mxxldir='/users/dapple/astro/mxxlsims/mxxl_imperial/snap41'

mcrs = 'c4 cfree'.split()
strategies = 'subaru_z4 pisco_z4_3 pisco_z4_4 pisco_z4_4square'.split()

mxxldata = {}

for mcr in mcrs:

    if mcr not in mxxldata:
        mxxldata[mcr] = {}


    for strategy in strategies:
        
        simname = '{0}-sour_{1}'.format(mcr, strategy)
        mxxlfilename='{0}/{1}/consolidated.pkl'.format(mxxldir, simname)

        mxxldata[mcr][strategy] = bbb.loadData(mxxlfilename)

   


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
    ccount = 0
    
    for strategy in strategies:





        axhline(1.0, c='DimGray', linestyle='--')

        simdat = mxxldata[mcr][strategy]
        massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
        plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], 
             label=strategy)

        ccount += 1





    xlabel(r'Mass $M_{200}$', fontsize=16)
    ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=16)
    axis([14.1, 15.2, 0.7, 1.3])
    legend(loc='upper left')
    title('MC Relation: {0}'.format(mcr))
    tight_layout()
    savefig('piscoproposal_{0}.png'.format(mcr))
    
####
## Scatter

for mcr in mcrs:
        
    figure()

    linecount = 0
    ccount = 0
    
    for strategy in strategies:




        axhline(1.0, c='DimGray', linestyle='--')

        simdat = mxxldata[mcr][strategy]
            
        massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
        plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount],
             label=strategy)
            
        ccount += 1
          


    
    xlabel(r'Mass $M_{200}$', fontsize=16)
    ylabel(r'Distribution Scatter', fontsize=16)
    axis([14.1, 15., 0.0, 0.6])
    legend(loc='upper left')
    title('Scatter - MC Relation: {0}'.format(mcr))
    tight_layout()
    savefig('piscoproposal_{0}_scatter.png'.format(mcr))

