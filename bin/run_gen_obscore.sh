#!/bin/bash
# This script runs the ObsCore generator for LSST generated images through Butler Gen 3.
#
cd "$(dirname "$0")"

if [ ! -z "$LSST_STACK" ]
then
  source $LSST_STACK/loadLSST.bash
  setup lsst_distrib;
fi

pip install -r ../requirements.txt
pip install ..

# generate the ObScore data for each dataset type
../python/lsst/dax/obscore_generator/gen_img_obscore.py --ds_types raw,calexp,deepCoadd --out_dir ./out

# combine the outputs
../python/lsst/dax/obscore_generator/combine_obscore_csv.py --out_dir ./out

