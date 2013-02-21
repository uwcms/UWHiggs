#!/bin/bash

# Generate the cython proxies used in the analyses

source jobid.sh

export jobid=$jobid8
export datasrc=`ls -d /scratch/*/data/$jobid | head -n 1`

if [ -z $1 ]; then
    export afile=`find $datasrc/ | grep root | head -n 1`
else
    export afile=$1
fi

echo "Building cython wrappers from file: $afile"

rake "make_wrapper[$afile, eet/final/Ntuple, EETauTree]"
rake "make_wrapper[$afile, emt/final/Ntuple, EMuTauTree]"
rake "make_wrapper[$afile, mmt/final/Ntuple, MuMuTauTree]"
rake "make_wrapper[$afile, mmm/final/Ntuple, MuMuMuTree]"
rake "make_wrapper[$afile, mm/final/Ntuple, MuMuTree]"
rake "make_wrapper[$afile, em/final/Ntuple, EMuTree]"
rake "make_wrapper[$afile, ee/final/Ntuple, EETree]"

ls *pyx | sed "s|pyx|so|" | xargs -n 1 -P 10 rake 
