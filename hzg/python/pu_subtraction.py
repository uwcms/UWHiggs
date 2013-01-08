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

def do_truth_reweight(dataPUTruth,datahisto,mcprob):
    data_histo_bin   = datahisto.FindBin(dataPUTruth)
    data_histo_prob = datahisto.GetBinContent(data_histo_bin)

    return data_histo_prob/mcprob/datahisto.Integral()
    
CD_file = TFile.Open(os.environ['CMSSW_BASE']+
                      '/src/UWHiggs/hzg/data/pu_truth_2012_CD.root')
gDirectory.cd(pwd)
CD_truth_histo = CD_file.Get('pileup').Clone()
CD_file.Close()
def pu_S10_CD_reweight(dataPUTruth):
    mc_histo_bin = min(int(floor(dataPUTruth)),S10.probFunctionVariable[-1])
    return do_truth_reweight(dataPUTruth,
                             CD_truth_histo,
                             S10.probValue[mc_histo_bin]/sum(S10.probValue))

    
