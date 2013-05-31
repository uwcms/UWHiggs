#!/bin/bash

# Make and copy all datacards

set -o nounset
set -o errexit

source jobid.sh
#export jobid=$jobid7
#rake cards
#rake limits
#rake harvest_limits
#rake plot_limits
#rake compare_limits

export jobid=$jobid8
rake cards
rake limits
rake harvest_limits
rake plot_limits
rake compare_limits