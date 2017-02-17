import os
import numpy as np
import cPickle
import yaml
import colossus.halo.concentration as chc
import colossusMassCon as cmc
import nfwutils

def createSims(bins, nperbin, zcluster, basedir, simgroup, idstart, cosmology):

    outputdir = '{}/{}'.format(basedir, simgroup)
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    nfwutils.global_cosmology.set_cosmology(cosmology)
    cmc.matchCosmo()
    
    #bins in delta=200
    #masses are M_sol*h

    m200s = []
    c200s = []

    answers = {}

    for mlow, mhigh in bins:

        curmasses = np.exp(np.random.uniform(np.log(mlow), np.log(mhigh), size=nperbin))
        m200s.append(curmasses)

        curconcens = chc.concentration(curmasses, '200c', zcluster, model='diemer15')
        c200s.append(curconcens)

    m200s = np.hstack(m200s)
    c200s = np.hstack(c200s)

    for i, id in enumerate(range(idstart, idstart+len(m200s))):

        config = dict(max_dist = 3.,
                      gridlength = 1024.,
                      zcluster = zcluster,
                      m200 = float(m200s[i]),
                      c200 = float(c200s[i]))

        with open('{}/analytic_{}.yaml'.format(outputdir, id), 'w') as output:
            yaml.dump(config, output)

        
        answers[id] = dict(m200 = m200s[i],
                             concen = c200s[i],
                             redshift = zcluster)


    with open('analytic_{}_answers.pkl'.format(simgroup), 'wb') as output:
        cPickle.dump(answers, output)
