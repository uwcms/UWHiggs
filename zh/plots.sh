#!/bin/bash
# Run all of the analysis

set -o nounset
set -o errexit

source jobid.sh

export jobid=$jobid7
rake plots
rake cards
rake limits
rake limitplots

export jobid=$jobid8
rake plots
rake cards
rake limits
rake limitplots