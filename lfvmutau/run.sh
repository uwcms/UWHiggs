#!/bin/bash

# Analyze the H->mutau channel

set -o nounset
set -o errexit

#export jobid=SubmitMuTauSingleMU
#export jobid=MuTauSingleMUV6
#export jobid=MuTauSingleMuJetReReco
#export jobid=MuTauSingleMuJetReReco
export jobid=MuTauSingleMuJetReRecoSkim
export PU=true
rake mttight
export PU=false
rake pu1mttight


