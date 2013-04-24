#!/bin/bash

# Make and copy all datacards

set -o nounset
set -o errexit

source jobid.sh
#export jobid=$jobid7
#rake cards

export jobid=$jobid8
rake cards
