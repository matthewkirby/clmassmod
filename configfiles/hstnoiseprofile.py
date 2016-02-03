import galaxypicker
galaxypicker=galaxypicker.AllGalaxyPicker()

import betacalcer
betacalcer=betacalcer.InfiniteRedshift()

import shearnoiser
shearnoiser=shearnoiser.NoNoise()

import hstnoiseBinning
hstbinning = hstnoiseBinning.HSTBinning()
binwidth=0.1
profilecol='r_mpc'
binner = hstnoiseBinning.HSTBinnerWrapper()
binnoiser = hstnoiseBinning.HSTBinNoiserWrapper()

