#!/usr/bin/env python

import glob

##HST
#sims = 'bk11snap124 mxxlsnap41'.split()
#rss = 'r5 r16'.split()
#mcs = 'c4 duffy diemer15'.split()
#centers = 'xrayNONE core szxvptcenter szlensingpeak xrayXVP xraylensingpeak'.split()
#deltas = (200, 500, 2500)
#

##Megacam
sims = 'bk11snap124 bk11snap141 mxxlsnap41 mxxlsnap54'.split()
rss = ['r9']
mcs = 'c4 duffy diemer15'.split()
centers = 'corenone sztcenter szxvptcenter core'.split()
deltas = (200, 500)

halosprocessed = {}
for line in open('haloprocessed').readlines():
    tokens = line.split()
    sim, config,loc = tokens
    if sim not in halosprocessed:
        halosprocessed[sim] = {}
    halosprocessed[sim][config] = loc

        

finished = {}
for line in open('/vol/euclid1/euclid1_2/dapple/rundlns/finished').readlines():
    tokens = line.split()
    sim = tokens[0]
    config = tokens[1]
    deltas = map(int, tokens[2:])
    if sim not in finished:
        finished[sim] = {}
    cursimlist = finished[sim]
    cursimlist[config] = deltas
    

needdln = {}
missinghalos = {}
for sim in finished.keys():
    missinghalos[sim] = []
    needdln[sim] = {}
    for delta in deltas:
        needdln[sim][delta] = []



                    
for sim in sims:
    cursimlist = finished[sim]
    for rs in rss:
        for mc in mcs:
            for center in centers:
#                for line in open('shearprofiles/coresizeindex.list').readlines():
#                    cluster, coreindex, coresize = line.split()
                for line in open('configfiles/megacam_siminput.reduced.list').readlines():
                    cluster, zcluster, ndensity, beta, core, coreindex = line.split()

                    curcenter = center
                    if center == 'core':
                        curcenter = 'core{}'.format(coreindex)

#                    config = 'hstnoisebins-{mc}-{rs}-{curcenter}-{cluster}'.format(mc = mc,          
#                                                                                   rs = rs,          
#                                                                                   curcenter = curcenter,  
#                                                                                   cluster = cluster)

                    config = 'mega-{mc}-{rs}-sigma0.25-{curcenter}-{cluster}'.format(mc = mc,
                                                                                   rs = rs,
                                                                                   curcenter = curcenter,
                                                                                   cluster = cluster)

                    if config not in cursimlist:
                        if config in halosprocessed[sim]:
                            for delta in deltas:
                                needdln[sim][delta].append(config)
                        else:
                            missinghalos[sim].append(config)
                    else:
                        for delta in deltas:
                            if delta not in cursimlist[config]:
                                needdln[sim][delta].append(config)


for sim in needdln.keys():
    for delta in deltas:
        with open('dlntorun.{}.{}'.format(sim, delta), 'w') as output:
            for config in needdln[sim][delta]:
                output.write('{}\n'.format(config))

for sim in missinghalos.keys():
    with open('missinghalos.{}'.format(sim),'w') as output:
        for config in missinghalos[sim]:
            output.write('{}\n'.format(config))



            

                    

                
