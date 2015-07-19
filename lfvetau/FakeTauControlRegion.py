##Correction Factor still to add
from ETauTree import ETauTree
import os
import ROOT
import math
import optimizer
import glob
import array
import mcCorrections
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, sin, cos, acos, sinh
from cutflowtracker import cut_flow_tracker
#Makes the cut flow histogram
cut_flow_step = ['allEvents', 'e1sel', 'e1IDiso', 'e2sel', 'e2IDiso', 'ZMass', 'tsel', 'tAntiMuon', 'tAntiEle',  'MtToMet', 'tRawIso10','tRawIso5',  'tLooseIso', 'tTightIso' 
]

from inspect import currentframe

def get_linenumber():
    cf = currentframe()
    return cf.f_back.f_lineno

def deltaPhi(phi1, phi2):
    PHI = abs(phi1-phi2)
    if PHI<=pi:
        return PHI
    else:
        return 2*pi-PHI
def deltaR(phi1, phi2, eta1, eta2):
    deta = eta1 - eta2
    dphi = abs(phi1-phi2)
    if (dphi>pi) : dphi = 2*pi-dphi
    return sqrt(deta*deta + dphi*dphi);
    
def etDR(row):
    return row.e1_t_DR if row.e1_t_DR < row.e2_t_DR else row.e2_t_DR

def etDPhi(row):
    e1tDPhi=deltaPhi(row.e1Phi, row.tPhi)
    e2tDPhi=deltaPhi(row.e2Phi, row.tPhi)
    return e1tDPhi if e1tDPhi < e2tDPhi else e2tDPhi

def Z(row):
    e1p=ROOT.TVector3(row.e1Pt*cos(row.e1Phi),row.e1Pt*sin(row.e1Phi),row.e1Pt*sinh(row.e1Eta))
    e2p=ROOT.TVector3(row.e2Pt*cos(row.e2Phi),row.e2Pt*sin(row.e2Phi),row.e2Pt*sinh(row.e2Eta))
    e1FourVector= ROOT.TLorentzVector(e1p, sqrt(e1p.Mag2()+row.e1Mass*row.e1Mass))
    e2FourVector= ROOT.TLorentzVector(e2p, sqrt(e2p.Mag2()+row.e2Mass*row.e2Mass))
    zFourVector = e1FourVector+e2FourVector
    return zFourVector


class FakeTauControlRegion(MegaBase):
    tree = 'et/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='ET'
        super(FakeTauControlRegion, self).__init__(tree, outfile, **kwargs)
        self.tree = ETauTree(tree)
        self.out=outfile
        self.histograms = {}
        self.pucorrector = mcCorrections.make_puCorrector('singlee')

        #optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith(self.channel) ]
        #self.grid_search = {}
        #if len(optimizer_keys) > 1:
        #    for key in optimizer_keys:
        #        self.grid_search[key] = optimizer.grid_search[key]
        #else:
        #    self.grid_search[''] = optimizer.grid_search[optimizer_keys[0]]


    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return 1.

            
            
            
        return self.pucorrector(row.nTruePU) * \
            mcCorrections.eid_correction( row, 'e') * \
            mcCorrections.eiso_correction(row, 'e') * \
            mcCorrections.trig_correction(row, 'e')
            
    def begin(self):
        
        tauiso = [ 'tLoose', 'tTigh']
        folder = []
        sign = ['ss','os']
        for iso in tauiso:
            for s in sign:
                folder.append(s+'/'+iso)
                folder.append(s+'/'+iso+'/tptregion')
                j=0
                while j < 4 :
                    folder.append(s+'/'+iso+'/'+str(j))
                    folder.append(s+'/'+iso+'/'+str(j)+'/tptregion')
                    j+=1
                    
        for f in folder: 
            
            self.book(f,"ePt", "e p_{T}", 200, 0, 200)
            self.book(f,"ePhi", "e phi",  100, -3.2, 3.2)
            self.book(f,"eEta", "e eta", 50, -2.5, 2.5)
            self.book(f,"tPt", "t p_{T}", 200, 0, 200)
            self.book(f,"tPhi", "t phi",  100, -3.2, 3.2)
            self.book(f,"tEta", "t eta", 50, -2.5, 2.5)
            self.book(f,"tAbsEta", "t abs eta", 50, -2.5, 2.5)
 
                     
    def fill_histos(self, row, folder='os/Loose', fakeRate = False):
        weight = self.event_weight(row)
        histos = self.histograms
 
        histos[folder+'/ePt'].Fill(row.ePt, weight)
        histos[folder+'/eEta'].Fill(row.eEta, weight)
        histos[folder+'/ePhi'].Fill(row.ePhi, weight) 

        histos[folder+'/tPt'].Fill(row.tPt, weight)
        histos[folder+'/tEta'].Fill(row.tEta, weight)
        histos[folder+'/tAbsEta'].Fill(abs(row.tEta), weight)
        histos[folder+'/tPhi'].Fill(row.tPhi, weight) 
 
    def process(self):
        
        #print self.tree.inputfilename
        for row in self.tree:
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            
            if not bool(row.singleE27WP80Pass) : continue
           
            if not bool(row.eMatchesSingleE27WP80) : continue
            
            if row.bjetCSVVeto30!=0 : continue 
            if not selections.eLowPtSelection(row, 'e'): continue
            if not selections.lepton_id_iso(row, 'e', 'eid13Loose_idantiso'): continue
            if abs(row.eEta) > 1.4442 and abs(row.eEta < 1.566) : continue
            
            
            if not selections.tauSelection(row, 't'): continue
            if row.tPt < 30 : continue
        
            if not row.tAntiMuon2Loose: continue
            if not row.tAntiElectronMVA5Tight: continue #was 3
            if row.tauVetoPt20EleTight3MuLoose : continue 
            #if row.tauHpsVetoPt20 : continue
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoCicLooseIso : continue # change it with Loose

            sign = 'ss' if row.e_t_SS else 'os'
            if  row.tLooseIso3Hits : 
                tauiso = 'tLoose'
                folder = sign+'/'+tauiso
                self.fill_histos(row,  folder)
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
                
            if row.tTightIso3Hits :
                tauiso = 'tTigh' 
                folder = sign+'/'+tauiso
                self.fill_histos(row,  folder)
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
                                                
             
            
    def finish(self):
        self.write_histos()


