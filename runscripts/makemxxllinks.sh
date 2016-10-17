#!/bin/bash


cat idconversion | awk '($1 !~ "#"){print}' | {

    while read cid sid proj; do 
	
	ln -s /vol/marvin2/raid6/shilbert/mxxl/lensing/halos/halos_054/halo_54_${sid}_${proj}_a_gauss_0.02_infinite_source_distance.convergence_map halo_cid${cid}.convergence_map
	
	ln -s /vol/marvin2/raid6/shilbert/mxxl/lensing/halos/halos_054/halo_54_${sid}_${proj}_a_gauss_0.02_infinite_source_distance.shear_1_map halo_cid${cid}.shear_1_map

	ln -s /vol/marvin2/raid6/shilbert/mxxl/lensing/halos/halos_054/halo_54_${sid}_${proj}_a_gauss_0.02_infinite_source_distance.shear_2_map halo_cid${cid}.shear_2_map

    done 

}