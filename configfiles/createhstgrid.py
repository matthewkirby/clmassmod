#!/usr/bin/env python
#######################

redshifts = [0.7, 0.8, 0.9, 1.0, 1.1]

#for i, z in enumerate(redshifts):
#
#    with open('hst_single_{0}.sh'.format(i), 'w') as output:
#
#        output.write('readermodule=readMXXL_HSTBeta\n')
#        output.write('readerclass=MXXLHSTSimReader\n')
#        output.write('maskname=acsmask\n')
#        output.write('maskx=0.\n')
#        output.write('masky=0.\n')
#        output.write('targetz={0}\n'.format(z))
#        output.write('nperarcmin={0}\n'.format(21))
#        output.write('shapenoise={0}\n\n'.format(0.4))
#
#
for i, z in enumerate(redshifts):

    with open('hst_deepmosaic_{0}.sh'.format(i), 'w') as output:

        output.write('readermodule=readMXXL_HSTBeta\n')
        output.write('readerclass=MXXLHSTSimReader\n')
        output.write('maskname=circlemask\n')
        output.write('maskx=0.\n')
        output.write('masky=0.\n')
        output.write('maskrad=2.75\n')
        output.write('targetz={0}\n'.format(z))
        output.write('nperarcmin={0}\n'.format(21))
        output.write('shapenoise={0}\n\n'.format(0.4))


for i, z in enumerate(redshifts):

    with open('hst_cvzsquare_{0}.sh'.format(i), 'w') as output:

        output.write('readermodule=readMXXL_HSTBeta\n')
        output.write('readerclass=MXXLHSTSimReader\n')
        output.write('maskname=squaremosaic\n')
        output.write('targetz={0}\n'.format(z))
        output.write('nperarcmin={0}\n'.format(16))
        output.write('shapenoise={0}\n\n'.format(0.4))

