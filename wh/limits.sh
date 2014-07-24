#!/bin/bash

# Make and copy all datacards

set -o nounset
set -o errexit

source jobid.sh
export jobid=$jobid7
rake f3_postfit_plots
rake compare_limits

export jobid=$jobid8 #'2013-Jun-30-8TeV' #
rake f3_postfit_plots
rake compare_limits
