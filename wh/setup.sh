#!/bin/bash

# Setup the cython proxies, find input ntuple files, and compute luminosity.

source jobid.sh
export datasrc=/scratch/efriis/data/
export jobid=$jobid7


rake "meta:getinputs[$jobid, $datasrc]"
rake "meta:getmeta[inputs/$jobid, mm/metaInfo, 7]"

export jobid=$jobid8
./make_proxies.sh
rake "meta:getinputs[$jobid, $datasrc]"
# Use the 7TeV WH samples for 8TeV
pushd inputs/$jobid/
# Symlink the list of input files and the counts of the number of events.
# For the effectively lumis, we have to recompute using the 8 TeV x-section.
ls ../../inputs/$jobid7/WH_*HWW* | grep -v lumicalc | xargs -n 1 ln -s 
popd
rake "meta:getmeta[inputs/$jobid, mm/metaInfo, 8]"
