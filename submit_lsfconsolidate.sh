#!/bin/bash
####################


torun=$1


for config in `cat $torun`; do 
    bsub -q long -oo fullrun/consolidate.$config.log ./consolidate_fits.py . bcc ~/nfs/bcc_clusters/recentered/$config
    bsub -q long -oo fullrun/consolidate.$config.log ./consolidate_fits.py . bk11snap141 ~/nfs/beckersims/snap141/intlength400/$config
    bsub -q long -oo fullrun/consolidate.$config.log ./consolidate_fits.py . bk11snap124 ~/nfs/beckersims/snap124/intlength400/$config

done