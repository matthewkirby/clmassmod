#!/bin/bash

export HOME=/home/dapple
export PATH=/vol/aibn218/data1/dapple/anaconda/bin:/home/dapple/bin:/vol/software/software/astro/theli/THELI//theli/gui/:/vol/software/software/astro/theli/THELI//theli/bin/Linux_64/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
export LD_LIBRARY_PATH=/home/dapple/lib
export SRCLOC=/vol/euclid1/euclid1_raid1/dapple/mxxlsims
#export PYTHONPATH=/home/dapple/braid1/mxxl/mxxlsims:/home/dapple/lib/python2.7/site-packages
python /vol/euclid1/euclid1_raid1/dapple/mxxlsims/rundln.py $@

