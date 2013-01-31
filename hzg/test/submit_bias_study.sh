#!/bin/bash

syncPostfix=${1}
ntoys=${2}



theCfg=$CMSSW_BASE/src/UWHiggs/hzg/test/do_bias_study.py
hdfsOutDir=srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN=/hdfs/store/user/lgray/HZG_bias_study/${syncPostfix}/

for channel in electron muon
  do
  hdfsInDir=root://cmsxrootd.hep.wisc.edu//store/user/lgray/HZG_bias_study/${channel}/

  for cat in 1 2 3 4
    do
    for mass in 120 125 130 135 140 145 150
      do
      for order in 3 4 5 6
      do
	for turnon in erf sigm
	  do
	  for truth in exp pow
	  do
	    fajOpts="--input-files-per-job=2 --infer-cmssw-path --express-queue --output-dir=${hdfsOutDir} --input-dir=${hdfsInDir} --job-generates-output-name --vsize-limit=30000 --fwklite"
	    jobName=hZg_bias_study.${channel}.c${cat}.m${mass}.o${order}.${turnon}.${truth}
	    farmoutAnalysisJobs $fajOpts $jobName $theCfg ${ntoys} ${cat} ${mass} ${channel} ${order} ${turnon} ${truth} inputFiles='$inputFileNames'
	  done
	done 
      done
    done
  done
done
