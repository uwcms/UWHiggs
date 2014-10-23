##Correction Factor still to add
from EETauTree import EETauTree
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


class TauFakeRateAnalyzerMVA(MegaBase):
    tree = 'eet/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EET'
        super(TauFakeRateAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        self.tree = EETauTree(tree)
        self.out=outfile
        self.histograms = {}
        self.pucorrector = mcCorrections.make_puCorrector('singlee')

        optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith(self.channel) ]
        self.grid_search = {}
        if len(optimizer_keys) > 1:
            for key in optimizer_keys:
                self.grid_search[key] = optimizer.grid_search[key]
        else:
            self.grid_search[''] = optimizer.grid_search[optimizer_keys[0]]


    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return 1.

        etrig = 'e1'
        if row.e2Pt > row.e1Pt : etrig = 'e2'
        if bool(row.e1MatchesSingleE27WP80) and  not bool(row.e2MatchesSingleE27WP80) : etrig = 'e1'
        if not bool(row.e1MatchesSingleE27WP80) and  bool(row.e2MatchesSingleE27WP80) :  etrig = 'e2'

        #if bool(row.e1MatchesEle27WP80) and  not bool(row.e2MatchesEle27WP80) : etrig = 'e1'
        #if not bool(row.e1MatchesEle27WP80) and  bool(row.e2MatchesEle27WP80) :  etrig = 'e2'
            
            
            
        return self.pucorrector(row.nTruePU) * \
            mcCorrections.eid_correction( row, 'e1', 'e2') * \
            mcCorrections.eiso_correction(row, 'e1', 'e2') * \
            mcCorrections.trig_correction(row, etrig     )
            
    def begin(self):
        
        tauiso = ['tNoCuts', 'tSuperSuperLoose', 'tSuperLoose', 'tLoose', 'tTigh']
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
            
            self.book(f,"e1Pt", "e1 p_{T}", 200, 0, 200)
            self.book(f,"e1Phi", "e1 phi",  100, -3.2, 3.2)
            self.book(f,"e1Eta", "e1 eta", 50, -2.5, 2.5)

            self.book(f,"e2Pt", "e2 p_{T}", 200, 0, 200)
            self.book(f,"e2Phi", "e2 phi",  100, -3.2, 3.2)
            self.book(f,"e2Eta", "e2 eta", 50, -2.5, 2.5)

            self.book(f, "e1e2Mass",  "e1e2 Inv Mass",  32, 0, 320)

            self.book(f, "tMtToPFMET", "#tau Met MT", 100, 0, 100)
            
            self.book(f,"tRawIso3Hits", "tRawIso3Hits", 500, 0, 100) 
            self.book(f,"tPt", "t p_{T}", 200, 0, 200)
            self.book(f,"tPtbarrel", "t p_{T} barrel", 200, 0, 200)
            self.book(f,"tPtendcap", "t p_{T} endcap", 200, 0, 200)
            self.book(f,"tPhi", "t phi",  100, -3.2, 3.2)
            self.book(f,"tEta", "t eta", 50, -2.5, 2.5)
            self.book(f,"tAbsEta", "t abs eta", 50, -2.5, 2.5)
 
            self.book(f,"etDR", "e t DR", 50, 0, 10)
            self.book(f,"etDPhi", "e t DPhi", 32, 0, 3.2)

            self.book(f,"ztDR", "Z #tau DR", 50, 0, 10)
            self.book(f,"ztDPhi", "Z #tau DPhi", 32, 0, 3.2)
            self.book(f,"Zpt", "Z p_{T}", 200, 0, 200)

            

            self.book(f, "type1_pfMetEt", "type1_pfMetEt",200, 0, 200) 
            self.book(f, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5)
            self.book(f, "bjetCSVVeto30", "number of bjets", 10, -0.5, 9.5)

        for s in sign:
            self.book(s+'/tNoCuts', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
            
            xaxis = self.histograms[s+'/tNoCuts/CUT_FLOW'].GetXaxis()
            self.cut_flow_histo = self.histograms[s+'/tNoCuts/CUT_FLOW']
            self.cut_flow_map   = {}
            for i, name in enumerate(cut_flow_step):
                xaxis.SetBinLabel(i+1, name)
                self.cut_flow_map[name] = i+0.5
                    
    def fill_histos(self, row, folder='os/tSuperLoose', fakeRate = False):
        weight = self.event_weight(row)
        histos = self.histograms
 
        histos[folder+'/e1Pt'].Fill(row.e1Pt, weight)
        histos[folder+'/e1Eta'].Fill(row.e1Eta, weight)
        histos[folder+'/e1Phi'].Fill(row.e1Phi, weight) 

        histos[folder+'/e2Pt'].Fill(row.e2Pt, weight)
        histos[folder+'/e2Eta'].Fill(row.e2Eta, weight)
        histos[folder+'/e2Phi'].Fill(row.e2Phi, weight)

        histos[folder+'/e1e2Mass'].Fill(row.e1_e2_Mass, weight)
        histos[folder+'/tMtToPFMET'].Fill(row.tMtToPFMET,weight)
    
        #histos[folder+'/tRawIso3Hits'].Fill(row.tRawIso3Hits, weight)
        histos[folder+'/tPt'].Fill(row.tPt, weight)
        if abs(row.tEta) < 1.5 :  histos[folder+'/tPtbarrel'].Fill(row.tPt, weight)
        if abs(row.tEta) > 1.5 :  histos[folder+'/tPtendcap'].Fill(row.tPt, weight)
        histos[folder+'/tEta'].Fill(row.tEta, weight)
        histos[folder+'/tAbsEta'].Fill(abs(row.tEta), weight)
        histos[folder+'/tPhi'].Fill(row.tPhi, weight) 
 
        histos[folder+'/type1_pfMetEt'].Fill(row.type1_pfMetEt)
        histos[folder+'/etDR'].Fill(etDR(row)) 
        histos[folder+'/etDPhi'].Fill(etDPhi(row)) 
        histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight) 
        histos[folder+'/bjetCSVVeto30'].Fill(row.bjetCSVVeto30, weight) 
       
        histos[folder+'/ztDR'].Fill(deltaR(Z(row).Phi(), row.tPhi, Z(row).Eta(), row.tEta))
        histos[folder+'/ztDPhi'].Fill(deltaPhi(Z(row).Phi(), row.tPhi))
        histos[folder+'/Zpt'].Fill(Z(row).Pt())
            

    def process(self):
        
        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)
        myevent =()
        #print self.tree.inputfilename
        for row in self.tree:
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            
            #if row.run > 2: 
            if not bool(row.singleE27WP80Pass) : continue
            #            if hasattr(self.tree, 'row.e1MatchesEle27WP80') and hasattr(self.tree, 'row.e2MatchesEle27WP80') :
            #if not bool(row.e1MatchesEle27WP80) and not bool(row.e2MatchesEle27WP80) : continue
            
            #else :
            if not bool(row.e1MatchesSingleE27WP80) and  not bool(row.e2MatchesSingleE27WP80) : continue
                #if not bool(row.singleEPass) : continue
                #if not bool(row.e1MatchesSingleE) and  not bool(row.e2MatchesSingleE) : continue
                
            if row.bjetCSVVeto30!=0 : continue 
            if row.e1Pt < 30 : continue
            if row.e2Pt < 30 : continue
                       
#        for i, row in enumerate(self.tree):
#            if  i >= 100:
#                return
           # print bool(cut_flow_trk.disabled)
            cut_flow_trk.new_row(row.run,row.lumi,row.evt)
            #print row.run,row.lumi,row.evt
            cut_flow_trk.Fill('allEvents')
            if not selections.eSelection(row, 'e1'): continue
            cut_flow_trk.Fill('e1sel')
            if not selections.lepton_id_iso(row, 'e1', 'eid13Tight_etauiso01'): continue
            if abs(row.e1Eta) > 1.4442 and abs(row.e1Eta < 1.566) : continue
            
            
            cut_flow_trk.Fill('e1IDiso')
            if not selections.eSelection(row, 'e2'): continue
            cut_flow_trk.Fill('e2sel')
            if not selections.lepton_id_iso(row, 'e2', 'eid13Tight_etauiso01'): continue
            if abs(row.e2Eta) > 1.4442 and abs(row.e2Eta) < 1.566 : continue
            
            cut_flow_trk.Fill('e2IDiso')
            if not abs(row.e1_e2_Mass-91.2) < 20: continue
            cut_flow_trk.Fill('ZMass')
            if not selections.tauSelection(row, 't'): continue
            cut_flow_trk.Fill('tsel')

            if not row.tAntiMuon2Loose: continue
            cut_flow_trk.Fill('tAntiMuon')
            if not row.tAntiElectronMVA5Tight: continue #was 3
            cut_flow_trk.Fill('tAntiEle')

            if row.tauVetoPt20EleTight3MuLoose : continue 
            #if row.tauHpsVetoPt20 : continue
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoCicLooseIso : continue # change it with Loose
            
            # if not row.tMtToMET < 50:  continue
            cut_flow_trk.Fill('MtToMet')
            
            #            if  etDR(row) < 1. : continue 
            if (row.run, row.lumi, row.evt, row.e1Pt, row.e2Pt)==myevent: continue
            myevent=(row.run, row.lumi, row.evt, row.e1Pt, row.e2Pt)

            tauiso = 'tNoCuts'
            sign = 'ss' if row.e1_e2_SS else 'os'
            folder = sign+'/'+tauiso
          
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)
            
            
            if  row.tPt < 90 and row.tPt>65 : 
                folder = folder+'/tptregion'
                self.fill_histos(row, folder)
                folder = sign+'/'+tauiso+'/tptregion'
                self.fill_histos(row, folder)
                
            if not row.tRawIso3Hits < 10 : continue
            cut_flow_trk.Fill('tRawIso10')
            tauiso = 'tSuperSuperLoose'
            folder = sign+'/'+tauiso
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)
            if  row.tPt < 90 and row.tPt>65 : 
                folder = folder+'/tptregion'
                self.fill_histos(row, folder)
                folder = sign+'/'+tauiso+'/tptregion'
                self.fill_histos(row, folder)
                
            if not row.tRawIso3Hits < 5 : continue
            cut_flow_trk.Fill('tRawIso5')
            tauiso = 'tSuperLoose'
            folder = sign+'/'+tauiso
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)
            if  row.tPt < 90 and row.tPt>65 : 
                folder = folder+'/tptregion'
                self.fill_histos(row, folder)
                folder = sign+'/'+tauiso+'/tptregion'
                self.fill_histos(row, folder)

                        
            if  row.tLooseIso3Hits : 
                tauiso = 'tLoose'
                folder = sign+'/'+tauiso
                self.fill_histos(row,  folder)
                cut_flow_trk.Fill('tLooseIso')
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
                if  row.tPt < 90 and row.tPt>65 : 
                    folder = folder+'/tptregion'
                    self.fill_histos(row, folder)
                    folder = sign+'/'+tauiso+'/tptregion'
                    self.fill_histos(row, folder)

            if row.tTightIso3Hits :
                tauiso = 'tTigh' 
                folder = sign+'/'+tauiso
                self.fill_histos(row,  folder)
                cut_flow_trk.Fill('tTightIso')
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
                if  row.tPt < 90 and row.tPt>65 : 
                    folder = folder+'/tptregion'
                    self.fill_histos(row, folder)
                    folder = sign+'/'+tauiso+'/tptregion'
                    self.fill_histos(row, folder)
                    
 
             
        cut_flow_trk.flush()
                                
             
            
    def finish(self):
        self.write_histos()


