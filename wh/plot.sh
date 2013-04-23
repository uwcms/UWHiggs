#!/bin/bash
# Plot all of the analysis

set -o nounset
set -o errexit

source jobid.sh
export jobid=$jobid7
rake control_zee 
rake eet_shapes

rake control_em
rake emt_shapes

rake control_zmm 
rake mmt_shapes
python plots_for_prepp.py 

export jobid=$jobid8
rake control_zee 
rake eet_shapes

rake control_em
rake emt_shapes

rake control_zmm 
rake mmt_shapes
python plots_for_prepp.py 
