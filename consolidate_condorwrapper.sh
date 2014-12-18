#!/bin/bash

export HOME=/home/dapple
export PATH=/home/dapple/anaconda/bin:/home/dapple/bin:/vol/software/software/astro/theli/THELI//theli/gui/:/vol/software/software/astro/theli/THELI//theli/bin/Linux_64/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
export LD_LIBRARY_PATH=/home/dapple/lib
export SRCLOC=/vol/braid1/vol1/dapple/mxxl/mxxlsims
#export PYTHONPATH=/home/dapple/braid1/mxxl/mxxlsims:/home/dapple/lib/python2.7/site-packages
python /home/dapple/braid1/mxxl/mxxlsims/consolidate_fits.py $@


