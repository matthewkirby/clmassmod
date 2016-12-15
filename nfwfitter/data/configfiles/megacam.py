import nfwfitter.profilebuilder as profilebuilder
profilebuilder=profilebuilder.ProfileBuilder()


massprior='linear'
import nfwfitter.rescalecluster as rescalecluster
rescalecluster=rescalecluster.RedshiftRescaler()

import nfwfitter.galaxypicker as galaxypicker
galaxypicker=galaxypicker.DensityPicker()
import nfwfitter.betacalcer as betacalcer
betacalcer=betacalcer.FixedBeta()
import nfwfitter.shearnoiser as shearnoiser
shearnoiser=shearnoiser.GaussianShapeNoise()

import nfwfitter.basicBinning as basicBinning
binner=basicBinning.GaussianFixedBins()
profilecol='r_mpc'
binspacing='linear'
nbins=12

import nfwfitter.binnoiser as binnoiser
binnoiser=binnoiser.NoBinNoise()

