#!/bin/bash
export OVERRIDE_META_TREE_data_ET='et/metaInfo'

export IGNORE_LUMI_ERRORS=1

source jobid.sh
#export jobid=$jobid8

export jobid='newNtuple_9Oct'
echo $jobid
export datasrc=/hdfs/store/user/$USER/  #$(ls -d /scratch/*/data/$jobid | awk -F$jobid '{print $1}')
#export datasrc=/nfs_scratch/taroni/data
export MEGAPATH=/hdfs/store/user/$USER
#export MEGAPATH=/nfs_scratch/taroni/data
./make_proxies.sh
rake "meta:getinputs[$jobid, $datasrc,et/metaInfo]"
rake "meta:getmeta[inputs/$jobid, et/metaInfo, 8]"
#export jobid=$jobidmt
#./make_proxies.sh
#rake "meta:getinputs[$jobid, $datasrc,em/metaInfo]"
#rake "meta:getmeta[inputs/$jobid, em/metaInfo, 8]"
#rake "meta:getinputs[$jobid, $datasrc,mt/metaInfo]"
#rake "meta:getmeta[inputs/$jobid, mt/metaInfo, 8]"


unset OVERRIDE_META_TREE_data_ET
unset IGNORE_LUMI_ERRORS
