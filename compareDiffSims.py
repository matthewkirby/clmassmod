import bcc_becker_baseline as bbb
import nfwutils


bccsimdir='/u/ki/dapple/nfs/bcc_clusters/recentered'

bk11snap141dir='/u/ki/dapple/nfs/beckersims/snap141/intlength400'
bk11snap124dir='/u/ki/dapple/nfs/beckersims/snap124/intlength400'

mxxldir='/u/ki/dapple/nfs22/mxxl/snap41'

#simdirs = [bccsimdir, bk11snap141dir, bk11snap124dir, mxxldir]
simdirs = [bccsimdir, bk11snap141dir, mxxldir]
#simnames = 'BCC BK11_z=0.25 BK11_z=0.5 MXXL_z=1.0'.split()
simnames = 'BCC BK11_z=0.5 MXXL_z=1.0'.split()

subdirs='c4_r5 cfree_r5'.split()
subdirnames='c4_r5 cfree_r5'.split()

data = {}

for simname, simdir in zip(simnames, simdirs):

    data[simname] = {}
    
    for subdir in subdirs:
            
            data[simname][subdir] = bbb.loadData('{0}/{1}/consolidated.pkl'.format(simdir, subdir))





c = 'y r b m c g'.split()


for subdirname, subdir in zip(subdirnames, subdirs):

    figure()
    ccount = 0
    patches = []
    labels = []

    for simname in simnames:

        simdat = data[simname]

        if simname == 'BCC':

#            zlow = simdat[subdir]['redshifts'] < 0.3
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat[subdir], selection = zlow)
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} z < 0.3'.format(simname))
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
#            labels.append('{0} z < 0.3'.format(simname))
#            ccount += 1
#

            zlow = np.logical_and(simdat[subdir]['redshifts'] > 0.36, simdat[subdir]['redshifts'] < 0.5)
            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat[subdir], selection = zlow)
            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} 0.36 < z < 0.5'.format(simname))
            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
            labels.append('{0} 0.36 < z < 0.5'.format(simname))
            ccount += 1

            
            zhigh = np.logical_and(simdat[subdir]['redshifts'] > 0.5, simdat[subdir]['redshifts'] < 0.89)
            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat[subdir], selection = zhigh)
            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} 0.5 < z < 0.89'.format(simname))
            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
            labels.append('{0} 0.5 < z < 0.89'.format(simname))
            ccount += 1


        else:
        
            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat[subdir])
            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label=simname)
            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
            labels.append(simname)
            ccount += 1



    

    
    xlabel(r'Mass $M_{200}$', fontsize=16)
    ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=16)
    axis([14.1, 15.2, 0.83, 1.12])
    legend(patches, labels, loc='upper left')
    title('MC Relation: {0}'.format(subdirname))
    savefig('{0}.png'.format(subdirname))
