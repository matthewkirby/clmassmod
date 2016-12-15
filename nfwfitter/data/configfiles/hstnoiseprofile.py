import nfwfitter.rescalecluster as rescalecluster
rescalecluster=rescalecluster.RedshiftRescaler()

import nfwfitter.galaxypicker as galaxypicker
galaxypicker=galaxypicker.AllGalaxyPicker()

import nfwfitter.betacalcer as betacalcer
betacalcer=betacalcer.FixedBeta()

import nfwfitter.shearnoiser as shearnoiser
shearnoiser=shearnoiser.NoNoise()

import nfwfitter.hstnoiseBinning as hstnoiseBinning
hstbinning = hstnoiseBinning.HSTBinning()
binwidth=0.1
profilecol='r_mpc'
binner = hstnoiseBinning.HSTBinnerWrapper()
binnoiser = hstnoiseBinning.HSTBinNoiserWrapper()

