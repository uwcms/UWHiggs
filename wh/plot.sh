#!/bin/bash
# Plot all of the analysis

set -o nounset
set -o errexit

source jobid.sh
export jobid=$jobid7
echo $jobid
python WHPlotterMMT.py
python WHPlotterEMT.py
python WHPlotterEET.py

python PlotControlZMM.py 
python PlotControlEM.py 
python PlotControlZEE.py 

python plots_for_prepp_kNN.py 

export jobid=$jobid8
echo $jobid
python WHPlotterMMT.py
python WHPlotterEMT.py
python WHPlotterEET.py

python PlotControlZMM.py 
python PlotControlEM.py 
python PlotControlZEE.py 

python plots_for_prepp_kNN.py 

