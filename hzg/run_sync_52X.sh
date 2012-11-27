#!/bin/bash

# A script to get a sync revision and run the ntuplization locally

srmPrefix="srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN"
remoteDataDir="/hdfs/store/user/lgray/HZG_sync/"
localDataDir="/scratch/$USER/data/HZG_sync/"
syncPostfix=${1}

if [ -z "$syncPostfix" ]
then
    echo "Please specify a postfix representing the sync revision!"
    pstfxs=`lcg-ls ${srmPrefix}=${remoteDataDir}`
    echo "Available postfixes: " 
    for fix in $pstfxs
    do
      echo " ${fix##*/}"
    done    
    exit
fi

mcSignalEle=`lcg-ls -b -Dsrmv2 ${srmPrefix}=${remoteDataDir}${syncPostfix}/52X/ | grep MCSignalEle | sed 's/\/hdfs\/store/root\:\/\/cmsxrootd\.hep\.wisc\.edu\/\/store/' | tr '\n' '\,' | head -c -1`
mcBkgEle=`lcg-ls -b -Dsrmv2 ${srmPrefix}=${remoteDataDir}${syncPostfix}/52X/ | grep MCBkgEle | sed 's/\/hdfs\/store/root\:\/\/cmsxrootd\.hep\.wisc\.edu\/\/store/' | tr '\n' '\,' | head -c -1`
dataElectron=`lcg-ls -b -Dsrmv2 ${srmPrefix}=${remoteDataDir}${syncPostfix}/52X/ | grep DataElectron | sed 's/\/hdfs\/store/root\:\/\/cmsxrootd\.hep\.wisc\.edu\/\/store/' | tr '\n' '\,' | head -c -1`
dataMuon=`lcg-ls -b -Dsrmv2 ${srmPrefix}=${remoteDataDir}${syncPostfix}/52X/ | grep DataMuon | sed 's/\/hdfs\/store/root\:\/\/cmsxrootd\.hep\.wisc\.edu\/\/store/' | tr '\n' '\,' | head -c -1`

theCfg=$CMSSW_BASE/src/FinalStateAnalysis/NtupleTools/test/make_ntuples_cfg.py

cmsRun  $theCfg inputFiles="$mcSignalEle" outputFile="MCSignalEle_52X_ntuples.root" makeHZG=1 makeDiObject=1 passThru=1 eventView=1 reportEvery=100 >& MCSignalEle_ntuples.log &

cmsRun  $theCfg inputFiles="$mcBkgEle" outputFile="MCBkgEle_52X_ntuples.root" makeHZG=1 makeDiObject=1 passThru=1 eventView=1 reportEvery=100 >& MCBkgEle_ntuples.log &

cmsRun  $theCfg inputFiles="$dataElectron" outputFile="DataElectron_52X_ntuples.root" makeHZG=1 makeDiObject=1 passThru=1 eventView=1 reportEvery=100 >& DataElectron_ntuples.log &

cmsRun  $theCfg inputFiles="$dataMuon" outputFile="DataMuon_52X_ntuples.root" makeHZG=1 makeDiObject=1 passThru=1 eventView=1 reportEvery=100 >& DataMuon_ntuples.log &

echo "Waiting for ntuplizing jobs to finish..."
wait