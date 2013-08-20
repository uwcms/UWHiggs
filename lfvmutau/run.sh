#!/bin/bash

set -o nounset
set -o errexit

export jobid=MuTauSingleMuJetReReco
export jobidPU=MuTauSingleMuJetReReco

rake analyzeMuTauData
rake analyzeMuTauMC
rake analyzeMuTauNOPU
