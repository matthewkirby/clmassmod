# clmassmod
A galaxy cluster weak-lensing mass modeling and verification code

This repository documents the work done within the SPT-WL working group and the LSST DESC Clusters working group to measure how weak lensing mass estimates relate to the true 3D mass of a cluster. The work focuses on emulating actual observations and mimicing procedures used by weak lensing observers. Mock observations are based on simulated lensing maps from n-body simulations.


Projects for the DESC Hack Day, 17 Feb 2017!

1.) Experiment with jointly fitting individual cluster shear profiles (Mass, concentration, miscentering, 2-halo amplitude) jointly as a population, ie a hierarchical model. Similar to [Lieu et al 2017](https://arxiv.org/abs/1701.00478), but we want to fit more parameters per cluster and fit an order of magnitude more clusters simultaneously. From a statistics & computation standpoint, can we sample O(2k) parameters with an MCMC engine, perhaps by taking advantage of the hierarchical structure of the problem?

2.) Daisuke Nagai & Erwin Lau, Nick Battaglia, and Stefan Hilbert & the Magneticum team have offered lensing maps produced from hydro simulations for analysis. Each simulation needs a reader that will bring the lensing information to a common format for CLMASSMOD use. We need to write that wrapper. We might also need to write an algorithm that coverts 2D mass maps to shear maps.


python setup.py develop
