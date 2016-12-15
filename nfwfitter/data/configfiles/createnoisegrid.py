#!/usr/bin/env python

for i, dens in enumerate([-1, 100, 20, 7, 4, 1]):
    for j, scat in enumerate([0., .16, .33, .5]):
        
        with open('n{0}_{1}.sh'.format(i,j), 'w') as output:
            output.write('nperarcmin={0}\n'.format(dens))
            output.write('shapenoise={0}\n\n'.format(scat))

