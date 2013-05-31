#!/bin/bash

# Get the data. To use my original ntuples, use login06 !
export jobid=SubmitMuTauSingleMU
export datasrc=/scratch/mcepeda/data

export jobid8TeV=$jobid
export afile=`find $datasrc/$jobid | grep root | head -n 1`

## Build the cython wrappers
rake "make_wrapper[$afile, mt/final/Ntuple, MuTauTree]"

ls *pyx | sed "s|pyx|so|" | xargs rake 

#Collect all the info about the samples
rake "meta:getinputs[$jobid, $datasrc,mt/metaInfo]"

#Compute the luminosity
rake "meta:getmeta[inputs/$jobid, mt/metaInfo, 8]"


