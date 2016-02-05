import profilebuilder
profilebuilder=profilebuilder.ProfileBuilder()


massprior='linear'
import rescalecluster
rescalecluster=rescalecluster.RedshiftRescaler()

import galaxypicker
galaxypicker=galaxypicker.DensityPicker()
import betacalcer
betacalcer=betacalcer.FixedBeta()
import shearnoiser
shearnoiser=shearnoiser.GaussianShapeNoise()

import basicBinning
binner=basicBinning.GaussianFixedBins()
profilecol='r_mpc'
binspacing='linear'
nbins=12

import binnoiser
binnoiser=binnoiser.NoBinNoise()

