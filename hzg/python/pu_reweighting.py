import ROOT
from ROOT import gDirectory,TFile,TH1D

import os
from math import floor

from FinalStateAnalysis.Utilities.version import cmssw_major_version,\
     cmssw_minor_version

pwd = gDirectory.GetPath()

puS7  = None
puS10 = None
puS6  = None

from SimGeneral.MixingModule.mix_E7TeV_Fall2011_Reprocess_50ns_PoissonOOTPU_cfi \
     import mix as S6
puS6 = S6.input.nbPileupEvents

if cmssw_major_version() == 5 and cmssw_minor_version()==2:
    from SimGeneral.MixingModule.mix_2012_Startup_50ns_PoissonOOTPU_cfi \
         import mix as S7
    puS7 = S7.input.nbPileupEvents
elif cmssw_major_version() == 5 and cmssw_minor_version()==3:
    from SimGeneral.MixingModule.mix_2012_Startup_50ns_PoissonOOTPU_cfi \
         import mix as S7
    puS7 = S7.input.nbPileupEvents
    from SimGeneral.MixingModule.mix_2012_Summer_50ns_PoissonOOTPU_cfi \
         import mix as S10
    puS10 = S10.input.nbPileupEvents

def do_truth_reweight(mcPUTruth,datahisto,mcprob):
    data_histo_bin  = datahisto.FindBin(mcPUTruth)
    data_histo_prob = ( datahisto.GetBinContent(data_histo_bin)/
                        datahisto.Integral() )
    return data_histo_prob/mcprob
    
CD_file = TFile.Open(os.environ['CMSSW_BASE']+
                      '/src/UWHiggs/hzg/data/pu_truth_2012_CD.root')
gDirectory.cd(pwd)
CD_truth_histo = CD_file.Get('pileup').Clone()
CD_file.Close()
def pu_S10_CD_reweight(mcPUTruth):
    bin_width = CD_truth_histo.GetBinWidth(1)
    mc_histo_bin = int(floor(mcPUTruth))
    
    return do_truth_reweight(mcPUTruth,
                             CD_truth_histo,
                             ( puS10.probValue[mc_histo_bin]/
                               sum(puS10.probValue) )*bin_width )


def clean_up_pu():
    del CD_truth_histo
