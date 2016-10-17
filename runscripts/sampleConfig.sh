##########################
# This defines the baseline configuration file to use with nfwfit.py
##########################

readermodule=readMXXL
readerclass=MXXLSimReader

nbootstraps=1


profileModule=basicBinning
profileBuilder=bootstrapfixedbins
binspacing='linear'
nbins=12

profileCol=r_mpc
profileMin=.5
profileMax=3.0 

ngals=200

massconModule=basicMassCon
massconRelation=constant
concentration=4.0


maskname=circlemask
maskx=0.
masky=0.
maskrad=20.
