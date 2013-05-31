#!/bin/bash 

echo "Setting up CMSSW runtime environment"
eval `scramv1 ru -sh`

echo "Sourcing FSA environment"
source $CMSSW_BASE/src/FinalStateAnalysis/environment.sh

echo "Sourcing HiggsAnalysis/HiggsToTauTau environment"
source $CMSSW_BASE/src/HiggsAnalysis/HiggsToTauTau/environment.sh

#Is the analysis blinded?
export blind='YES'

#check if dev area is up to date
check_git_updates.sh
