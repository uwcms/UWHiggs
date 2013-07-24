#!/bin/bash
# Run all of the analysis

set -o nounset
set -o errexit

source jobid.sh

export jobid=$jobid8
rake fakerates
rake fits

#rake mmcontrol
#rake emcontrol
rake eecontrol

rake mmt
rake emt
rake eet

#export jobid=$jobid7
#rake fakerates
#rake fits

#rake mmcontrol
#rake emcontrol
#rake eecontrol

#rake mmt
#rake emt
#rake eet

