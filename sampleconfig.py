








import readMXXL
simreader = readMXXL.MXXLSimReader()

import nfwfit
model = nfwfit.NFW_MC_Model()
massprior='linear'

import colossusMassCon
massconRelation = colossusMassCon.ColossusMC()
colossusMCname='diemer15'

import nfwfit
fitter = nfwfit.NFWFitter()
fitmethod='pdf'


import galaxypicker
densitypicker = galaxypicker.DensityPicker()
nperarcmin=10.
targetz=0.85
maskname='squaremosaic'

import galaxypicker
fovpicker = galaxypicker.FoVPicker()
maskname='acsmask'

galaxypicker = galaxypicker.Composite(galaxypicker,fovpicker)

import betacalcer
betacalcer = betacalcer.FixedBeta()
beta_s = 0.5

import shearnoiser
shearnoiser = shearnoiser.GaussianShapeNoise()
shapenoise = 0.25

import centergenerator
centergenerator = centergenerator.XrayMagneticumOffset()

import basicBinning
binner = basicBinning.gaussianfixedbins()
profilemax = 1.5
profilemin = 0.5
binspacing = 'linear'
nbins = 12
profilecol = 'r_mpc'


import binnoiser
binnoiser = binnoiser.NoBinNoise()

import profilebuilder
profilebuilder = profilebuilder.ProfileBuilder()











