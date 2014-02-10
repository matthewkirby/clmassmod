import bcc_becker_baseline as bbb
import nfwutils


mxxldir='/users/dapple/astro/mxxlsims/mxxl_imperial/snap41'

mcrs = 'c4 duffy'.split()
radialranges = ['r{0}'.format(i) for i in range(1,11)]


#data = {}
#
#for mcr in mcrs:
#
#    data[mcr] = {}
#
#    for radrange in radialranges:
#
#        simname = '{0}_{1}'.format(mcr, radrange)
#        filename='{0}/{1}/consolidated.pkl'.format(mxxldir, simname)
#
#        data[mcr][radrange] = bbb.loadData(filename)
#
#
#


c = 'r b m c g'.split()
linestyle=': -- -'.split()

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
'r10' :'0.75 - 3.0 Mpc'
}

####
# Bias

for mcr in mcrs:

    figure()


    linecount = 0
    ccount = 0

    for currad, radrange in enumerate(radialranges):

        simdat = data[mcr][radrange]

        if currad == 0:
            linecount = 0
            ccount = 3
        if currad == 1:
            linecount = 0
            ccount = 0

        if currad == 4:
            linecount = 1
            ccount = 0
        if currad == 7:
            linecount = 2
            ccount = 0

            

        
        massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
        plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], label=radrange_trans[radrange])
        ccount += 1





    

    
    xlabel(r'Mass $M_{200}$', fontsize=16)
    ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=16)
    axis([14.1, 15.2, 0.7, 1.3])
    legend(loc='upper left')
    title('MC Relation: {0}'.format(mcr))
    savefig('{0}.png'.format(mcr))

###
# Scatter

for mcr in mcrs:

    figure()


    linecount = 0
    ccount = 0

    for currad, radrange in enumerate(radialranges):

        simdat = data[mcr][radrange]

        if currad == 0:
            linecount = 0
            ccount = 3
        if currad == 1:
            linecount = 0
            ccount = 0

        if currad == 4:
            linecount = 1
            ccount = 0
        if currad == 7:
            linecount = 2
            ccount = 0

            

        massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
        plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount], label=radrange_trans[radrange])
        ccount += 1





    

    
    xlabel(r'Mass $M_{200}$', fontsize=16)
    ylabel(r'Distribution Scatter', fontsize=16)
    axis([14.1, 15., 0.0, 0.6])
    legend(loc='upper left')
    title('Scatter - MC Relation: {0}'.format(mcr))
    savefig('{0}_scatter.png'.format(mcr))
