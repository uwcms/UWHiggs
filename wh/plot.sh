#!/bin/bash
# Plot all of the analysis

set -o nounset
set -o errexit

source jobid.sh
#export jobid=$jobid7
#rake plot_zee 
#rake plot_eet
#
#rake plot_em
#rake plot_emt
#
#rake plot_zmm 
#rake plot_mmt
#python plots_for_prepp.py 
#
export jobid=$jobid8
#rake plot_zee 
rake plot_eet
#
#rake plot_em
rake plot_emt
#
#rake plot_zmm 
rake plot_mmt
python plots_for_prepp.py 
