#!/bin/bash
# Run all of the analysis

set -o nounset
set -o errexit

source jobid.sh

export jobid=$jobid8
rake mmcontrol
rake emcontrol

export jobid=$jobid7
rake mmcontrol
rake emcontrol

export jobid=$jobid8
rake fakerates
rake mmt
rake emt

export jobid=$jobid7
rake fakerates
rake mmt
rake emt
