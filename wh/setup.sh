#!/bin/bash

# Setup the cython proxies, find input ntuple files, and compute luminosity.

source jobid.sh
export jobid=$jobid7
export datasrc=$(ls -d /scratch/*/data/$jobid | awk -F$jobid '{print $1}')

for dir in $datasrc; do
    echo $dir
    rake "meta:getinputs[$jobid, $dir, mm/metaInfo]"
    rake "meta:getmeta[inputs/$jobid, mm/metaInfo, 7]"
done

rm inputs/$jobid/WZJetsTo3LNu_ZToTauTau*
for file in $(ls inputs/$jobid/WZJetsTo3LNu*); do
    newname=`echo $file | sed 's|WZJetsTo3LNu|WZJetsTo3LNu_ZToTauTau|'`
    cp -uv $file $newname
done

export jobid=$jobid8
export datasrc=$(ls -d /scratch/*/data/$jobid | awk -F$jobid '{print $1}')
#./make_proxies.sh
for dir in $datasrc; do
    echo $dir
    rake "meta:getinputs[$jobid, $dir, mm/metaInfo]"
    rake "meta:getmeta[inputs/$jobid, mm/metaInfo, 8]"
done

rm inputs/$jobid/WZJetsTo3LNu_ZToTauTau*
for file in $(ls inputs/$jobid/WZJetsTo3LNu*); do
    newname=`echo $file | sed 's|WZJetsTo3LNu|WZJetsTo3LNu_ZToTauTau|'`
    cp -uv $file $newname
done

# Use the 7TeV WH samples for 8TeV
#pushd inputs/$jobid/
# Symlink the list of input files and the counts of the number of events.
# For the effectively lumis, we have to recompute using the 8 TeV x-section.
#ls ../../inputs/$jobid7/WH_*HWW* | grep -v lumicalc | xargs -n 1 ln -s 
#popd
