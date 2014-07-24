#!/bin/bash
# Run all of the analysis

set -o nounset
set -o errexit

source jobid.sh

export jobid=$jobid8
rake fakerates
rake charge_fakes
rake kNN

rake mmt
rake emt
rake eet

rake mmcontrol
rake emcontrol
rake eecontrol

export jobid=$jobid7
rake fakerates
rake charge_fakes
rake kNN
#rake fits

rake mmt
rake emt
rake eet

rake mmcontrol
rake emcontrol
rake eecontrol
