import ROOT
import os
from sys import argv, stdout, stderr
import sys
import glob
import logging

#jobid = os.environ['jobid']
jobid = 'newNtuple_2Sept_old'
from os import listdir
from os.path import isfile, join

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadGridX(True)

from ETauTree import ETauTree
import math

##oldntuple
#f=ROOT.TFile("/hdfs/store/user/taroni/newNtuple_3June/data_SingleElectron_Run2012B_22Jan2013_v1/make_ntuples_cfg-patTuple_cfg-00FB7F10-C77D-E211-B334-003048CF9B4E.root")
##newNtuple
#f = ROOT.TFile("/hdfs/store/user/taroni/newNtuple_2Sept/data_SingleElectron_Run2012B_22Jan2013_v1/make_ntuples_cfg-patTuple_cfg-00FB7F10-C77D-E211-B334-003048CF9B4E.root")
##mc
#f= ROOT.TFile("/hdfs/store/user/taroni/newNtuple_2Sept/Zjets_M50_skimmedLL/make_ntuples_cfg-patTuple_cfg-00037C53-AAD1-E111-B1BE-003048D45F38.root")
#more trees in ntuple
#f= ROOT.TFile("/hdfs/store/user/taroni/testneNtupleMoreTrees/Zjets_M50/make_ntuples_cfg-patTuple_cfg-00037C53-AAD1-E111-B1BE-003048D45F38.root")
#one tree only 
#f= ROOT.TFile("/hdfs/store/user/taroni/testneNtuple/Zjets_M50/make_ntuples_cfg-patTuple_cfg-00037C53-AAD1-E111-B1BE-003048D45F38.root")
f= ROOT.TFile("/hdfs/store/user/taroni/newNtuple_3Dec/data_SingleElectron_Run2012A_22Jan2013_v1/make_ntuples_cfg-patTuple_cfg-08172766-C973-E211-B14C-0030487D8541.root")

f.ls()
mytree = f.Get("et/final/Ntuple")

for n, event in enumerate(mytree) :
    if n > 100 : break
    print event.run, event.lumi, event.evt, event.tPhi, event.ePhi, event.tEta, event.eEta, event.tPt, event.ePt
