import bcc_becker_baseline as bbb
import nfwutils


mxxldir='/users/dapple/astro/mxxlsims/mxxl_imperial/snap41'

mcrs = 'cfree c4'.split()
radialranges = ['r{0}'.format(i) for i in '1 11 12 2 3 4'.split()]



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



outerrange = [0.5, 0.75, 1.0, 1.5, 2.5, 3.0]

c = 'r b m DodgerBlue g DarkSalmon'.split()

linestyle='- -- :'.split()

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





###################################
# Mass Instrinsc Scatter


#massbins = []
#scatters = []
#
#for radrange in radialranges:
#
#
#    simdat = data['c4'][radrange]
#
#    massbin, massbin_min, massbin_max, intrinsic_scatter, median_err, sym_distro, logratiodists, sigmadists = bbb.scatterSummary(simdat)
#
#    massbins.append(massbin)
#    scatters.append(sym_distro)
#
#

figsize(4.5,4)
#
#figure()
#
#mlow500 = nfwutils.Mdelta(nfwutils.rscaleConstM(10**massbins[0][0], 4.0, 1.0, 200), 4.0, 1.0, 500)/ 1e14
#mid500 = nfwutils.Mdelta(nfwutils.rscaleConstM(10**massbins[0][8], 4.0, 1.0, 200), 4.0, 1.0, 500)/ 1e14
#mhigh500 = nfwutils.Mdelta(nfwutils.rscaleConstM(10**massbins[0][-1], 4.0, 1.0, 200), 4.0, 1.0, 500)/ 1e14
#
#
#
##plot(outerrange[:-1], [scatter[0] for scatter in scatters[:-1]], c='SlateGray', linewidth=1.5, label='%1.1f x10^14' % mlow500)
##plot(outerrange[:-1], [scatter[-1] for scatter in scatters[:-1]], c='k', linewidth = 3.0, label='%1.1f x10^14' % mhigh500)
#plot(outerrange[:-1], [scatter[8] for scatter in scatters[:-1]], c='k', linewidth = 3.0, label='%1.1f x10^14' % mid500)
#
#axvline((1.9/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(0.75), c='r', linestyle=':', linewidth=2.)
#axvline((1.9/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(1.15), c='r', linestyle='--', linewidth=2.)
#text(0.76, 0.23, '1xACS', fontsize=14, fontweight='semibold')
#
#axvline((2.65/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(0.75), c='g', linestyle=':', linewidth=2.)
#axvline((2.65/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(1.15), c='g', linestyle='--', linewidth=2.)
#text(1.105, 0.196, '3xACS', fontsize=14, fontweight='semibold')
#
#axvline((3.3/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(0.75), c='b', linestyle=':', linewidth=2.)
#axvline((3.3/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(1.15), c='b', linestyle='--', linewidth=2.)
#text(1.405, 0.16, '4xACS', fontsize=14, fontweight='semibold')
#
#axis([0.5, 2.0, 0.15, 0.55])
#
#ax = gca()
#ax.set_xticks(map(float, '0.5 1.0 1.5 2.0'.split()))
#ax.set_xticklabels('0.5 1.0 1.5 2.0'.split(), fontsize=16)
#
#ax.set_yticks(np.arange(0.15, 0.60, 0.1))
#ax.set_yticklabels(['%1.2f' % x for x in np.arange(0.15, 0.60, 0.1)], fontsize=16)
#
#
#
#xlabel('Outer Fit Radius [Mpc]', fontsize=18)
#ylabel('Intrinsic Scatter', fontsize=18)
#legend()
#
#
#tight_layout()


#############################
# Bias Uncertainty


#massbins2 = []
#offset = []
#
#for radrange in radialranges:
#
#
#    simdat1 = data['c4'][radrange]
#    simdat2 = data['cfree'][radrange]
#
#    massbin1, ratiobin1, ratioerr1 = bbb.summary2DMass(simdat1)
#    massbin2, ratiobin2, ratioerr2 = bbb.summary2DMass(simdat2)
#
#    massbins2.append(massbin1)
#    offset.append((np.array(ratiobin1)/ np.array(ratiobin2)) - 1.)
#
##




figure()

mlow500 = nfwutils.Mdelta(nfwutils.rscaleConstM(10**massbins[0][0], 4.0, 1.0, 200), 4.0, 1.0, 500)/ 1e14
mid500 = nfwutils.Mdelta(nfwutils.rscaleConstM(10**massbins[0][8], 4.0, 1.0, 200), 4.0, 1.0, 500)/ 1e14
mhigh500 = nfwutils.Mdelta(nfwutils.rscaleConstM(10**massbins[0][-1], 4.0, 1.0, 200), 4.0, 1.0, 500)/ 1e14



plot(outerrange[:-1], 100*np.array([scatter[8] for scatter in scatters[:-1]]), c='k', linewidth = 1.8, label='Intrinsic Scatter')

twincolor = '#B366B3'
#twincolor = 'SlateGray'
plot(outerrange[:-1], 100*np.array([np.abs(x[8]) for x in offset[:-1]]), c=twincolor, linewidth = 3.0, label='M-C Error Bound')


axvline((1.9/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(0.75), c='r', linestyle=':', linewidth=2.)
axvline((1.9/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(1.15), c='r', linestyle='--', linewidth=2.)
text(0.78, 100*0.17, '1xACS', fontsize=14, fontweight='semibold')

axvline((2.65/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(0.75), c='g', linestyle=':', linewidth=2.)
axvline((2.65/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(1.15), c='g', linestyle='--', linewidth=2.)
text(1.12, 100*0.17, '3xACS', fontsize=14, fontweight='semibold')

axvline((3.3/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(0.75), c='b', linestyle=':', linewidth=2.)
axvline((3.3/60)*(np.pi/180.)*nfwutils.global_cosmology.angulardist(1.15), c='b', linestyle='--', linewidth=2.)
text(1.425, 100*0.17, '4xACS', fontsize=14, fontweight='semibold')

axis([0.72, 1.8, 0.0, 100*0.45])

ax = gca()
minorticks_on()

ax.set_xticks(map(float, '0.75 1.25 1.75'.split()))
ax.set_xticks(map(float, '1.0 1.5'.split()), minor=True)
ax.set_xticklabels('0.75 1.25 1.75'.split(), fontsize=16)

ax.set_yticks(np.arange(0.0, 100*0.50, 100*0.1))
ax.set_yticks(np.arange(100*0.05, 100*0.50, 100*0.1), minor=True)
ax.set_yticklabels(['%2.0f' % x for x in np.arange(0.0, 100*0.50, 100*0.1)], fontsize=16)

ax.tick_params('both', length=5, width=1, which = 'major')
ax.tick_params('both', length=2.5, width=0.5, which = 'minor')

xlabel('Outer Fit Radius [Mpc]', fontsize=18)
ylabel('% Uncertainty', fontsize=18)
legend()

ax.set_xlim(0.72, 1.8)



#ax2.set_ylabel('M-C Relation Uncertainty', fontsize=18, color=twincolor)
#
#ax2.set_xlim(0.72, 1.8)
#ax2.set_ylim(0.0, 0.13)
#
#ax2.set_yticks(np.arange(0.0, 0.13, 0.03))
#ax2.set_yticklabels(['%1.2f' % x for x in np.arange(0.0, 0.13, 0.03)], fontsize=16)
#
#for tl in ax2.get_yticklabels():
#    tl.set_color(twincolor)
#
#
#


tight_layout()



