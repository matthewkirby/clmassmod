import bcc_becker_baseline as bbb
import nfwutils


mxxldir='/users/dapple/astro/mxxlsims/mxxl_imperial/snap41'

mcrs = 'c4 cfree'.split()
radialranges = ['r{0}'.format(i) for i in '11 5 6 10'.split()]
densities = 'n0_0 n3_0'.split()
binnings = 'equalbins50 equalbins200 logbins6 logbins12 logbins18 linearbins6 linearbins12 linearbins18'.split()


#data = {}
#
#for mcr in mcrs:
#
#    data[mcr] = {}
#
#    for radrange in radialranges:
#
#        data[mcr][radrange] = {}
#
#        for density in densities:
#
#            data[mcr][radrange][density] = {}
#
#            for binning in binnings:
#
#                simname = '{0}-{1}-{2}-{3}'.format(mcr, radrange, density, binning)
#                filename='{0}/{1}/consolidated.pkl'.format(mxxldir, simname)
#
#                data[mcr][radrange][density][binning] = bbb.loadData(filename)
#
#        
#




c = 'r b m DodgerBlue g DarkSalmon BlueViolet Teal'.split()
#linestyle=': -- -'.split()
linestyle='-- - :'.split()

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


figsize(8,8)

####
# Bias

for mcr in mcrs:

    for radrange in radialranges:

        figure()

        ccount = 0
        for binning in binnings:

            
            linecount = 0
            for density in densities:



                simdat = data[mcr][radrange][density][binning]            


        
                massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
                axhline(1.0, c='DimGray', linestyle='--')
                plot(massbin, ratiobin, color=c[ccount], linestyle=linestyle[linecount], label='{0} {1}'.format(binning, density))

    
                linecount += 1

            ccount += 1

    
        xlabel(r'Mass $M_{200}$', fontsize=16)
        ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=16)
        axis([14.1, 15.2, 0.7, 1.3])
        legend(loc='upper left')
        figlabel = '{0} {1}'.format(mcr, radrange_trans[radrange])
        title(figlabel)
        tight_layout()
        savefig('{0}.png'.format(figlabel))

###
# Scatter

for mcr in mcrs:

    for radrange in radialranges:

        figure()

        ccount = 0
        for binning in binnings:

            
            linecount = 0
            for density in densities:



                simdat = data[mcr][radrange][density][binning]            


                massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
                plot(massbin, sym_distro, color=c[ccount], linestyle=linestyle[linecount], label='{0} {1}'.format(binning, density))

                linecount += 1

            ccount += 1

        xlabel(r'Mass $M_{200}$', fontsize=16)
        ylabel(r'Distribution Scatter', fontsize=16)
        axis([14.1, 15., 0.0, 0.6])
        legend(loc='upper left')
        figlabel = '{0} {1}'.format(mcr, radrange_trans[radrange])
        title('Scatter - {0}'.format(figlabel))
        tight_layout()
        savefig('{0}_scatter.png'.format(figlabel))


    


    

    
