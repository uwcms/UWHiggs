University of Wisconsin CMS Higgs Analysis
==========================================

Contains final plotting/analysis code for UW H→ττ analyses.

Installation
------------

Start from a CMSSW 5_2_5 or greater working area.
Clone the code from the remote repository::

   git clone --recursive https://github.com/ekfriis/UWHiggs.git

The --recursive option is necessary for the dependencies, which are stored as
git submodules, to be downloaded correctly.  Now check out/install all the
extra dependencies::

   cmsenv
   kinit YOUR_CERN_USERNAME@CERN.CH 
   cd UWHiggs
   ./install.sh
   scram b -j 8


Quick Start
-----------

To run the WH analysis::
   
   source FinalStateAnalysis/environment.sh
   cd UWHiggs/wh
   ./setup.sh # run once
   ./run.sh  # run anytime analyzers are changed
   ./plot.sh # plot results

