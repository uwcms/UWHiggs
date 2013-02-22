#!/bin/bash 

echo "Setting up CMSSW runtime environment"
eval `scramv1 ru -sh`
source $CMSSW_BASE/src/FinalStateAnalysis/environment.sh

#Is the analysis blinded?
export blind='YES'

#check if dev area is up to date
check_git_updates.sh
