import pylab
import matplotlib.patches as mpatches
import astropy.io.ascii as asciireader
import numpy as np





#data = asciireader.read('hstbiassummary_2040')

def plotOne(data, center, mc, rs, delta, fig = None):

    if fig is None:
        fig = pylab.figure()

    ax = pylab.gca()

    selected = reduce(np.logical_and, (data['center'] == center,
                                       data['mc'] == mc,
                                       data['rad'] == rs,
                                       data['delta'] == delta))

    sims = 'BK11124', 'MXXL41'
    simredshifts = np.array([0.49925070, 0.9887081240357551])
    offsets = -0.02, 0.02
    markers = 'o ^'.split()
    colors = 'b r'.split()

    clusters = []
    redshifts = []
    sortedbiases = []
    sortedbiaserrs = []

    for sim, offset, marker, color in zip(sims, offsets, markers, colors):
        
        records = data[np.logical_and(selected, data['sim'] == sim)]

        if len(records) == 0:
            sortedbiases.append([])
            continue
        
        redshiftorder = np.argsort(records['zcluster'])

        redshifts = records['zcluster'][redshiftorder]
        clusters = records['cluster'][redshiftorder]

        sortedbiases.append(records['b'][redshiftorder])
        sortedbiaserrs.append(records['b_err'][redshiftorder])

        ax.errorbar(records['zcluster'][redshiftorder] + offset, 
                    records['b'][redshiftorder],
                    records['b_err'][redshiftorder],
                    linestyle='none',
                    marker=marker,
                    c = color,
                    label = sim)


    interpted_bias = np.zeros_like(redshifts)
    interpted_biaserr = np.zeros_like(redshifts)
    delta_bias = np.zeros_like(redshifts)

    if len(sortedbiases[0]) == len(sortedbiases[1]) and len(redshifts) > 0:

        for i in range(len(redshifts)):
            xs = np.array([redshifts[i] + offset for offset in offsets])
            ys = np.array([biasvals[i] for biasvals in sortedbiases])
            yerrs = np.array([biaserrs[i] for biaserrs in sortedbiaserrs])
            ax.plot(xs, ys, 'k:')

            delta_bias[i] = ys[1] - ys[0]

            #interp bias
            if redshifts[i] <= simredshifts[0]:
                interpted_bias[i] = ys[0]
                interpted_biaserr[i] = yerrs[0]
            elif redshifts[i] >= simredshifts[1]:
                interpted_bias[i] = ys[1]
                interpted_biaserr[i] = yerrs[1]
            else:
                m = (ys[1] - ys[0])/(simredshifts[1] - simredshifts[0])
                interpted_bias[i] = m*(redshifts[i] - simredshifts[0]) + ys[0]
                interpted_biaserr[i] = np.sqrt(((xs[1] - redshifts[i])*yerrs[0])**2 + \
                                               ((redshifts[i] - xs[0])*yerrs[1])**2)/(xs[1] - xs[0])
                

            
                                          



    ax.axhline(1.0, c='k')
    ax.set_xlabel('Redshift', fontsize=16)
    ax.set_ylabel('Bias', fontsize=16)
    ax.set_title('{} {} {} {}'.format(center, mc, rs, delta), fontsize=16)
    fig.tight_layout()




    

    

    return fig, clusters, redshifts, interpted_bias, interpted_biaserr, delta_bias
        
####


centers = 'xrayNONE xrayXVP szxvptcenter core%d xraylensingpeak szlensingpeak'.split()
mcs = 'c4 duffy'.split()
rss = 'r5 r16'.split()
deltas = [500, 200]


def doAll():

    data = asciireader.read('hstbiassummary_nocomments')

    figs = []

    with open('hstbiassummary_reduced', 'w') as output:

        output.write('cluster zcluster rad mc delta center b b_err b_delta\n')

        for delta in deltas:

            for rs in rss:
                for center in centers:
                    for mc in mcs:

                        results = plotOne(data, center, mc, rs, delta = delta)
                        fig, clusters, redshifts, interpted_bias, interpted_biaserr, delta_bias = results
                        ax = pylab.gca()
                        if rs == 'r5':
                            ax.set_ylim(0.65, 1.15)
                        elif rs == 'r16':
                            ax.set_ylim(0.5, 0.9)

                        for cluster, redshift, bias, biaserr, deltabias in zip(clusters, 
                                                                               redshifts, 
                                                                               interpted_bias, 
                                                                               interpted_biaserr,
                                                                               delta_bias):
                            output.write('%s %1.2f %s %s %d %s %1.4f %1.4f %1.4f\n' % (cluster,
                                                                                       redshift,
                                                                                       rs,
                                                                                       mc,
                                                                                       delta,
                                                                                       center,
                                                                                       bias,
                                                                                       biaserr,
                                                                                       deltabias))
                                                                             
                                                       


                        fig.savefig('hstbiassummary_plots/summary.%s.%s.%s.%d.png' % (rs, center, mc, delta))

                        figs.append(fig)

    return figs



##################

#pubplotdata = asciireader.read('hstbiassummary_reduced')

def pubplots():

    
    #xray

    xrayfig = pylab.figure()
    ax = xrayfig.add_subplot(1,1,1)

    centers = 'xrayNONE xrayXVP'.split()
    center_labels = ('Perfect', 'Empirical')

    mcs = 'c4 duffy'.split()
    mc_labels = ('c=4', 'Duffy et al. 2008')

    offsets = (-0.02, 0.02)
    colors = ('k', 'r')

    ax.axhline(1.0, c='k', linestyle='--')    

    for centeri, center in enumerate(centers):

        xs = []
        ys = []
        yerrs = []


        for mci, mc in enumerate(mcs):

            selection = reduce(np.logical_and, (pubplotdata['rad'] == 'r5',
                                                pubplotdata['mc'] == mc,
                                                pubplotdata['delta'] == 500,
                                                pubplotdata['center'] == center))
            
            selected = pubplotdata[selection]

            x = selected['zcluster'] + offsets[mci]
            y = selected['b']
            yerr = selected['b_err']

            labelargs = {}
            if centeri == 0:
                labelargs = dict(label=mc_labels[mci])
            
            if mci == 0:
                ax.errorbar(x,y,yerr, 
                            color=colors[centeri], 
                            marker='o', 
                            linestyle='none', 
                            **labelargs)
            elif mci == 1:
                ax.errorbar(x,y,yerr, 
                            color = colors[centeri],
                            markerfacecolor='none', 
                            markeredgecolor=colors[centeri],
                            marker = 'o', 
                            linestyle='none',
                            **labelargs)

            xs.append(x)
            ys.append(y)
            yerrs.append(yerr)


        if len(xs[0]) == len(xs[1]):
            
            for i in range(len(xs[0])):
                ax.plot([x[i] for x in xs], 
                        [y[i] for y in ys],
                        marker=None,
                        linestyle=':',
                        color = colors[centeri])

    ax.set_xlabel('Cluster Redshift', fontsize=18)
    ax.set_ylabel('Bias', fontsize=18)
    ax.set_xlim(0.5, 1.2)
    ax.set_ylim(0.5, 1.2)

    mclegend = ax.legend(loc='lower right', 
                       numpoints=1, 
                       frameon=False)


    patches = [mpatches.Patch(color=colors[i], label=center_labels[i]) for i in range(len(center_labels))]
    centerlegend = ax.legend(patches, center_labels, 
                             loc='lower left', 
                             numpoints=1, 
                             frameon=False, 
                             title='X-ray Center')
    pylab.setp(centerlegend.get_title(), fontsize=16)

    ax.add_artist(mclegend) #add back first legend

    xrayfig.tight_layout()
    xrayfig.savefig('hstbiassummary_plots/xray_bias_summary.png')
    xrayfig.savefig('hstbiassummary_plots/xray_bias_summary.eps')
    xrayfig.savefig('hstbiassummary_plots/xray_bias_summary.ps')
    xrayfig.savefig('hstbiassummary_plots/xray_bias_summary.pdf')


    #SZ

    szfig = pylab.figure()
    ax = szfig.add_subplot(1,1,1)

    centers = 'xrayNONE szxvptcenter core%d'.split()
    center_labels = ('Perfect', 'Empirical', 'Hydro Sim')

    mcs = 'c4 duffy'.split()
    mc_labels = ('c=4', 'Duffy et al. 2008')

    offsets = (-0.02, 0.0, 0.02)
    colors = ('k', 'r', 'b')

    
    ax.axhline(1.0, c='k', linestyle='--')
    

    for centeri, center in enumerate(centers):

        xs = []
        ys = []
        yerrs = []


        for mci, mc in enumerate(mcs):

            selection = reduce(np.logical_and, (pubplotdata['rad'] == 'r5',
                                                pubplotdata['mc'] == mc,
                                                pubplotdata['delta'] == 500,
                                                pubplotdata['center'] == center))
            
            selected = pubplotdata[selection]

            x = selected['zcluster'] + offsets[mci]
            y = selected['b']
            yerr = selected['b_err']

            labelargs = {}
            if centeri == 0:
                labelargs = dict(label=mc_labels[mci])
            
            if mci == 0:
                ax.errorbar(x,y,yerr, 
                            color=colors[centeri], 
                            marker='o', 
                            linestyle='none', 
                            **labelargs)
            elif mci == 1:
                ax.errorbar(x,y,yerr, 
                            color = colors[centeri],
                            markerfacecolor='none', 
                            markeredgecolor=colors[centeri],
                            marker = 'o', 
                            linestyle='none',
                            **labelargs)

            xs.append(x)
            ys.append(y)
            yerrs.append(yerr)


        if len(xs[0]) == len(xs[1]):
            
            for i in range(len(xs[0])):
                ax.plot([x[i] for x in xs], 
                        [y[i] for y in ys],
                        marker=None,
                        linestyle=':',
                        color = colors[centeri])

    ax.set_xlabel('Cluster Redshift', fontsize=18)
    ax.set_ylabel('Bias', fontsize=18)
    ax.set_xlim(0.5, 1.2)
    ax.set_ylim(0.5, 1.2)

    mclegend = ax.legend(loc='lower right', 
                       numpoints=1, 
                       frameon=False)


    patches = [mpatches.Patch(color=colors[i], label=center_labels[i]) for i in range(len(center_labels))]
    centerlegend = ax.legend(patches, center_labels, 
                             loc='lower left', 
                             numpoints=1, 
                             frameon=False, 
                             title='SZ Center')
    pylab.setp(centerlegend.get_title(), fontsize=16)

    ax.add_artist(mclegend) #add back first legend


    szfig.tight_layout()
    szfig.savefig('hstbiassummary_plots/sz_bias_summary.png')
    szfig.savefig('hstbiassummary_plots/sz_bias_summary.eps')
    szfig.savefig('hstbiassummary_plots/sz_bias_summary.ps')
    szfig.savefig('hstbiassummary_plots/sz_bias_summary.pdf')



    return xrayfig, szfig
        

            
        


                


    
              
