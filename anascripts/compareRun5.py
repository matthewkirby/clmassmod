import bcc_becker_baseline as bbb
import nfwutils
import re


rundir = '/users/dapple/astro/mxxlsims/mxxl_imperial/rundirs/run5consolidated'

sims = 'bcc bk11snap141 bk11snap124 mxxlsnap54 mxxlsnap41'.split()
simnames = 'BCC BK11_z=0.25 BK11_z=0.5 MXXL_z=0.25 MXXL_z=1.0'.split()


#mcs = 'cfree c4 zhaomc'.split()
mcs = ['cfree']

#radranges = 'r8 r9 r10'.split()
radranges = ['r10']
#samplings = 'n0 n3 n4'.split()
samplings = ['n0']


##########################
#Load Data

#data = {}
#
#for sim in sims:
#
#    for mc in mcs:
#
#        for radrange in radranges:
#
#            for sampling in samplings:
#
#                key = (sim, mc, radrange, sampling)
#
#                name ='%s.%s-%s-%s_0.consolidated.pkl' % key
#                filename='%s/%s' % (rundir, name)
#
#                try:
#                    data[key] = bbb.loadData(filename)
#                except:
#                    print 'Cannot load %s' % filename
#    
#

################################

#c = 'y r b m c g DodgerBlue DarkSalmon'.split()
c = [(.9,.6,0), (.35, .7, .9), (0,.6,.5)]
hatch = [None, '-']

########
## Q1: What is the spread between sims using different mc relations
#
#    
#for mc in mcs:
#
#    figure()
#    ccount = 0
#    patches = []
#    labels = []
#
#
#    for sim, simname in zip(sims, simnames):
#
#
#
#
#        key = (sim, mc, 'r10', 'n0') 
#
#        if key not in data:
#            print 'Missing %s' % ' '.join(key)
#            continue
#
#        simdat = data[key]
#
#        print
#        print key
#
#
#        if simname == 'BCC':
#
#
#            zlow = simdat['redshifts'] < 0.36
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat, selection = zlow)
#
#            print len(simdat['redshifts'][zlow]), len(massbin)
#
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} z < 0.36'.format(simname))
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
#            labels.append('{0} z < 0.36'.format(simname))
#            ccount += 1
#            
#
#            zlow = np.logical_and(simdat['redshifts'] > 0.36, simdat['redshifts'] < 0.5)
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat, selection = zlow)
#
#            print len(simdat['redshifts'][zlow]), len(massbin)
#
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} 0.36 < z < 0.5'.format(simname))
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
#            labels.append('{0} 0.36 < z < 0.5'.format(simname))
#            ccount += 1
#
#            
#            zhigh = np.logical_and(simdat['redshifts'] > 0.5, simdat['redshifts'] < 0.89)
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat, selection = zhigh)
#
#            print len(simdat['redshifts'][zhigh]), len(massbin)
#
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} 0.5 < z < 0.89'.format(simname))
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
#            labels.append('{0} 0.5 < z < 0.89'.format(simname))
#            ccount += 1
#
#
#        else:
#            
#            if ccount == 2:
#                ccount = 0
#            
#
#            curhatch = hatch[0]
#            if re.match('MXXL', simname) is not None:
#                curhatch = hatch[1]
#
#            if curhatch is not None and ccount == 1:
#                ccount += 1
#        
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
#
#            print len(simdat['redshifts']), len(massbin)
#
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label=simname, hatch = curhatch)
#
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3, hatch = curhatch))
#            labels.append(simname)
#            ccount += 1
#
#
#
#    
#
#    
#    xlabel(r'Mass $M_{200}$', fontsize=18)
#    ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=18)
#    axis([14.2, 15.6, 0.83, 1.25])
#    legend(patches, labels, loc='upper left')
#    title('MC Relation: {0}'.format(mc), fontsize=18)
#    savefig('plots_run5/{0}.png'.format(mc))
#
#
#######
## Q2: Redshift Slices
#
#
#### z=0.25    


def fillerrs(m, merr, y, yerr):
    
    xpoints = []
    y_low = []
    y_high = []

    for i in range(len(m)):

        xpoints.append(m[i]-merr[0][i])
        xpoints.append(m[i]+merr[1][i])
        y_low.append(y[i]-yerr[i])
        y_low.append(y[i]-yerr[i])
        y_high.append(y[i]+yerr[i])
        y_high.append(y[i]+yerr[i])
        
    return xpoints, y_low, y_high

######
#
#
sims = 'bcc bk11snap141 mxxlsnap54'.split()
simnames = 'BCC BK11_z=0.25 MXXL_z=0.25'.split()
curhatch = None


for mc in mcs:

    figure()
    ccount = 0
    patches = []
    labels = []


    axhline(1.0, c='k', linewidth=2.5, linestyle='--')


    for sim, simname in zip(sims, simnames):




        key = (sim, mc, 'r10', 'n0') 

        if key not in data:
            print 'Missing %s' % ' '.join(key)
            continue

        simdat = data[key]

        print
        print key


        if simname == 'BCC':


            zlow = simdat['redshifts'] < 0.36
            massbin, massedge, ratiobin, ratioerr = bbb.summary2DMass(simdat, selection = zlow)

            

            print len(simdat['redshifts'][zlow]), len(massbin)

#            errorbar(massbin, ratiobin, ratioerr, massedge, marker='None', color=c[ccount], linestyle='None')
            xpoints, y_low, y_high = fillerrs(massbin, massedge, ratiobin, ratioerr)
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} z < 0.36'.format(simname))
            fill_between(xpoints, y_low, y_high, alpha=0.8, color = c[ccount], label='{0} z < 0.36'.format(simname))
            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.8))
            labels.append('{0} z < 0.36'.format(simname))
            ccount += 1


            

        else:
                    
            massbin, massedge, ratiobin, ratioerr = bbb.summary2DMass(simdat)
            xpoints, y_low, y_high = fillerrs(massbin, massedge, ratiobin, ratioerr)
            print len(simdat['redshifts']), len(massbin)
#            errorbar(massbin, ratiobin, ratioerr, massedge, marker='None', color=c[ccount], linestyle='None')
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label=simname, hatch = curhatch)
            fill_between(xpoints, y_low, y_high, alpha=0.8, color = c[ccount], label=simname, hatch = curhatch)
            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.8, hatch = curhatch))
            labels.append(simname)
            ccount += 1

        print simname, ratioerr/ratiobin

    

    

    
    xlabel(r'Mass $M_{200}$', fontsize=18)
    ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=18)
    axis([14.2, 15.4, 0.85, 1.15])
    legend(patches, labels, loc='upper left')
    title('Z = 0.25 MC Relation: {0}'.format(mc), fontsize=18)
    savefig('{0}_z25.png'.format(mc))

##
#### z=0.5    
#
#sims = 'bcc bk11snap124'.split()
#simnames = 'BCC BK11_z=0.5'.split()
#
#curhatch = None
#
#
#for mc in mcs:
#
#    figure()
#    ccount = 0
#    patches = []
#    labels = []
#
#
#    for sim, simname in zip(sims, simnames):
#
#
#
#
#        key = (sim, mc, 'r10', 'n0') 
#
#        if key not in data:
#            print 'Missing %s' % ' '.join(key)
#            continue
#
#        simdat = data[key]
#
#        print
#        print key
#
#
#        if simname == 'BCC':
#
#
#            zlow = np.logical_and(simdat['redshifts'] > 0.36, simdat['redshifts'] < 0.5)
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat, selection = zlow)
#
#            print len(simdat['redshifts'][zlow]), len(massbin)
#
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} z < 0.36'.format(simname))
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
#            labels.append('{0} 0.36 < z < 0.5'.format(simname))
#            ccount += 1
#            
#
#        else:
#                    
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
#
#            print len(simdat['redshifts']), len(massbin)
#
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label=simname, hatch = curhatch)
#
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3, hatch = curhatch))
#            labels.append(simname)
#            ccount += 1
#
#
#
#    
#
#    
#    xlabel(r'Mass $M_{200}$', fontsize=18)
#    ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=18)
#    axis([14.2, 15.6, 0.83, 1.25])
#    legend(patches, labels, loc='upper left')
#    title('Z = 0.5 MC Relation: {0}'.format(mc), fontsize=18)
#    savefig('{0}.z5.png'.format(mc))
#
#
## z=1.0    
#
#sims = 'bcc mxxlsnap41'.split()
#simnames = 'BCC MXXL_z=1.0'.split()
#
#
#curhatch = None
#
#
#for mc in mcs:
#
#    figure()
#    ccount = 0
#    patches = []
#    labels = []
#
#
#    axhline(1.0, c='k', linewidth=2.5, linestyle='--')
#
#
#    for sim, simname in zip(sims, simnames):
#
#
#
#
#        key = (sim, mc, 'r10', 'n0') 
#
#        if key not in data:
#            print 'Missing %s' % ' '.join(key)
#            continue
#
#        simdat = data[key]
#
#        print
#        print key
#
#
#        if simname == 'BCC':
#
#
#            zhigh = np.logical_and(simdat['redshifts'] > 0.5, simdat['redshifts'] < 0.89)
#            massbin, massedge, ratiobin, ratioerr = bbb.summary2DMass(simdat, selection = zhigh)
#
#            
#
#            print len(simdat['redshifts'][zlow]), len(massbin)
#
##            errorbar(massbin, ratiobin, ratioerr, massedge, marker='None', color=c[ccount], linestyle='None')
#            xpoints, y_low, y_high = fillerrs(massbin, massedge, ratiobin, ratioerr)
##            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} z < 0.36'.format(simname))
#            fill_between(xpoints, y_low, y_high, alpha=0.8, color = c[ccount])
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.8))
#            labels.append('{0} 0.5 < z < 0.89'.format(simname))
#            ccount += 1
#
#
#            
#
#        else:
#
#            ccount += 1
#                    
#            massbin, massedge, ratiobin, ratioerr = bbb.summary2DMass(simdat)
#            xpoints, y_low, y_high = fillerrs(massbin, massedge, ratiobin, ratioerr)
#            print len(simdat['redshifts']), len(massbin)
##            errorbar(massbin, ratiobin, ratioerr, massedge, marker='None', color=c[ccount], linestyle='None')
##            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label=simname, hatch = curhatch)
#            fill_between(xpoints, y_low, y_high, alpha=0.8, color = c[ccount], label=simname, hatch = curhatch)
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.8, hatch = curhatch))
#            labels.append(simname)
#            ccount += 1
#
#        print simname, ratioerr/ratiobin
#
#    
#
#
#    xlabel(r'Mass $M_{200}$', fontsize=18)
#    ylabel(r'Ratio $M_{\mathrm{lensing}} / M_{\mathrm{true}}$', fontsize=18)
#    axis([14.2, 15.4, 0.85, 1.15])
#    legend(patches, labels, loc='upper left')
#    title('Z ~ 1.0 MC Relation: {0}'.format(mc), fontsize=18)
#    savefig('{0}_z10.png'.format(mc))
#    
#

#
#
#for mc in mcs:
#
#    figure()
#    ccount = 0
#    patches = []
#    labels = []
#
#
#    for sim, simname in zip(sims, simnames):
#
#
#
#
#        key = (sim, mc, 'r10', 'n0') 
#
#        if key not in data:
#            print 'Missing %s' % ' '.join(key)
#            continue
#
#        simdat = data[key]
#
#        print
#        print key
#
#
#        if simname == 'BCC':
#

#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat, selection = zhigh)
#
#            print len(simdat['redshifts'][zhigh]), len(massbin)
#
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label='{0} 0.5 < z < 0.89'.format(simname))
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3))
#            labels.append('{0} 0.5 < z < 0.89'.format(simname))
#            ccount += 1
#            
#
#        else:
#
#            ccount += 1
#                    
#            massbin, ratiobin, ratioerr = bbb.summary2DMass(simdat)
#
#            print len(simdat['redshifts']), len(massbin)
#
#            fill_between(massbin, ratiobin + ratioerr, ratiobin - ratioerr, alpha=0.3, color = c[ccount], label=simname, hatch = curhatch)
#
#            patches.append(Rectangle((0, 0), 1, 1, fc=c[ccount], alpha=0.3, hatch = curhatch))
#            labels.append(simname)
#            ccount += 1
#
#
#
#    
#
#    
