import profilebuilder
profilebuilder=profilebuilder.ProfileBuilder()


massprior=linear
import rescalecluster
rescalecluster=rescalecluster.RedshiftRescaler()

import galaxypicker
galaxypicker=galaxypicker.DensityPicker()
import betacalcer
betacalcer=betacalcer.FixedBeta()
import shearnoiser
shearnoiser=shearnosier.GaussianShapeNoise()

import basicBinning
binner=GaussianFixedBins()
profileCol='r_mpc'
binspacing='linear'
nbins=12
