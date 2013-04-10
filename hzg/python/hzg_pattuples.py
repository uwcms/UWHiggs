#! /usr/bin/env python

import os

patTupleRoot = os.environ['hzgpattupleroot']
patTupleDate = os.environ['hzgpattupledate']

sevenTeVTuples1 = {
    'tupleName':'7TeV-53X-PatTuple_Master',
    'tupleDate':patTupleDate,
    'tupleRoot':patTupleRoot,
    'ggHToZG':"""
    /GluGluToHToZG_M-120_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /GluGluToHToZG_M-125_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /GluGluToHToZG_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /GluGluToHToZG_M-135_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /GluGluToHToZG_M-140_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /GluGluToHToZG_M-145_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /GluGluToHToZG_M-150_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    """.split(),
    'VBFHToZG':"""
    /VBF_HToZG_M-120_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /VBF_HToZG_M-125_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /VBF_HToZG_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /VBF_HToZG_M-135_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /VBF_HToZG_M-140_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /VBF_HToZG_M-145_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /VBF_HToZG_M-150_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    """.split(),
    'VHToZG':"""
    /WH_ZH_HToZG_M-120_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /WH_ZH_HToZG_M-125_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /WH_ZH_HToZG_M-130_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /WH_ZH_HToZG_M-135_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /WH_ZH_HToZG_M-140_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /WH_ZH_HToZG_M-145_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /WH_ZH_HToZG_M-150_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    """.split(),
    'ttHToZG':"""
    /TTH_HToZG_M-120_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /TTH_HToZG_M-125_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /TTH_HToZG_M-130_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /TTH_HToZG_M-135_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /TTH_HToZG_M-140_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /TTH_HToZG_M-145_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    /TTH_HToZG_M-150_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM
    """.split(),
    'ZGToLLG':"""
    /ZGToEEG_TuneZ2_7TeV-madgraph-upto2jets/lgray-Fall11-REDIGI-AODSIM-START42_V14B-v1_ZGToEEG-RECO-01715716b3165466edf30580d661ec8b/USER
    /ZGToMuMuG_TuneZ2_7TeV-madgraph-upto2jets/lgray-Fall11-REDIGI-AODSIM-START42_V14B-v1_ZGToMuMuG-RECO-01715716b3165466edf30580d661ec8b/USER
    """.split()
    }

sevenTeVTuples2 = {
    'tupleName':'7TeV-53X-PatTuple',
    'tupleDate':patTupleDate,
    'tupleRoot':patTupleRoot,
    'Zjets_M50':"""
    /DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM
    """.split(),
    'data_DoubleElectron_2011AB_16Jan2012_v1':"""
    /DoubleElectron/Run2011A-16Jan2012-v1/AOD
    /DoubleElectron/Run2011B-16Jan2012-v1/AOD
    """.split(),
    'data_DoubleMu_2011AB_16Jan2012_v1':"""
    /DoubleMu/Run2011A-16Jan2012-v1/AOD
    /DoubleMu/Run2011B-16Jan2012-v1/AOD
    """.split()
    }

eightTeVTuples = {
    'tupleName':'8TeV-53X-PatTuple_Master',
    'tupleDate':patTupleDate,
    'tupleRoot':patTupleRoot,
    'ggHToZG':"""
    /GluGluToHToZG_M-120_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /GluGluToHToZG_M-125_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /GluGluToHToZG_M-130_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /GluGluToHToZG_M-135_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /GluGluToHToZG_M-140_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /GluGluToHToZG_M-145_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /GluGluToHToZG_M-150_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    """.split(),
    'VBFHToZG':"""
    /VBF_HToZG_M-120_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /VBF_HToZG_M-125_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /VBF_HToZG_M-130_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /VBF_HToZG_M-135_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /VBF_HToZG_M-140_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /VBF_HToZG_M-145_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /VBF_HToZG_M-150_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    """.split(),
    'VHToZG':"""
    /WH_ZH_HToZG_M-120_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /WH_ZH_HToZG_M-125_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /WH_ZH_HToZG_M-130_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /WH_ZH_HToZG_M-140_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /WH_ZH_HToZG_M-145_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /WH_ZH_HToZG_M-150_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    """.split(),
    'ttHToZG':"""
    /TTH_HToZG_M-120_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /TTH_HToZG_M-125_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /TTH_HToZG_M-130_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /TTH_HToZG_M-135_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /TTH_HToZG_M-140_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /TTH_HToZG_M-145_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    /TTH_HToZG_M-150_8TeV-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    """.split(),
    'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball':"""
    /DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    """.split(),
    'ZGToLLG':"""
    /ZGToLLG_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM
    """.split(),
    'data_DoubleElectron_2012AB_13Jul2012':"""
    /DoubleElectron/Run2012A-13Jul2012-v1/AOD
    /DoubleElectron/Run2012B-13Jul2012-v1/AOD
    """.split(),
    'data_DoubleElectron_2012C':"""
    /DoubleElectron/Run2012C-24Aug2012-v1/AOD
    /DoubleElectron/Run2012C-PromptReco-v2/AOD
    /DoubleElectron/Run2012C-EcalRecover_11Dec2012-v1/AOD
    """.split(),
    'data_DoubleElectron_2012D':"""
    /DoubleElectron/Run2012D-PromptReco-v1/AOD
    """.split(),
    'data_DoubleMu_2012AB':"""
    /DoubleMu/Run2012A-13Jul2012-v1/AOD
    /DoubleMu/Run2012B-13Jul2012-v4/AOD
    """.split(),
    'data_DoubleMu_2012C':"""
    /DoubleMu/Run2012C-24Aug2012-v1/AOD
    /DoubleMu/Run2012C-PromptReco-v2/AOD
    /DoubleMu/Run2012C-EcalRecover_11Dec2012-v1/AOD    
    """.split(),
    'data_DoubleMu_2012D':"""
    /DoubleMu/Run2012D-PromptReco-v1/AOD
    """.split()
    }

allTuples = {'7TeVpart1':sevenTeVTuples1,
             '7TeVpart2':sevenTeVTuples2,
             '8TeV':eightTeVTuples}

#for tuple in sevenTeVTuples.split():
#    print '%s%s/%s-%s'%(patTupleRoot,tuple,patTupleDate,sevenTeV)
#
#for tuple in eightTeVTuples.split():
#    print '%s%s/%s-%s'%(patTupleRoot,tuple,patTupleDate,eightTeV)
