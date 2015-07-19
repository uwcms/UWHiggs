##Correction Factor still to add
from EETree import EETree
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
    

class NewEleFakeRateFromW(MegaBase):
    tree = 'ee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EE'
        super(NewEleFakeRateFromW, self).__init__(tree, outfile, **kwargs)
        self.tree = EETree(tree)
        self.out=outfile
        self.histograms = {}
        self.pucorrector = mcCorrections.make_puCorrector('singlee')
        self.mye1 = 'e1'
        self.mye2 = 'e2'
        #optimizer_keys   = [ i for i in optimizer.grid_search.keys() if i.startswith(self.channel) ]
        self.grid_search = {}
        #if len(optimizer_keys) > 1:
        #    for key in optimizer_keys:
        #        self.grid_search[key] = optimizer.grid_search[key]
        #else:
        #    self.grid_search[''] = optimizer.grid_search[optimizer_keys[0]]


    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return 1.

 
        return self.pucorrector(row.nTruePU) * \
            mcCorrections.eid_correction( row, 'e1', 'e2') * \
            mcCorrections.eiso_correction(row, 'e1', 'e2') * \
            mcCorrections.trig_correction(row, 'e2'  )
         


##add the trigger correction 

    def begin(self):
        
        eiso = ['eLoose', 'eTigh']
        folder = []
        sign = ['ss','os']
        for iso in eiso:
            for s in sign:
                folder.append(s+'/'+iso)
                j=0
                while j < 4 :
                    folder.append(s+'/'+iso+'/'+str(j))
                    j+=1
                    
        for f in folder: 
            
            self.book(f,"e1Pt", "e1 p_{T}", 200, 0, 200)
            ##self.book(f,"e1Phi", "e1 phi",  100, -3.2, 3.2)
            ##self.book(f,"e1Eta", "e1 eta", 46, -2.3, 2.3)
            ##
            self.book(f,"e2Pt", "e2 p_{T}", 200, 0, 200)
            self.book(f,"e2Phi", "e2 phi",  100, -3.2, 3.2)
            self.book(f,"e2Eta", "e2 eta", 46, -2.3, 2.3)
            self.book(f,"e2AbsEta", "e2 abs eta", 23, 0, 2.3)
            self.book(f,"e2Pt_vs_e2AbsEta", "e2 pt vs e2 abs eta", 23, 0, 2.3,  20, 0, 200.,  type=ROOT.TH2F)


            self.book(f, "e1e2Mass",  "e1e2 Inv Mass",  32, 0, 320)
            self.book(f, "e1MtToPFMET",  "e1MtToPFMET", 20, 0, 200 )
            self.book(f, "e2MtToPFMET",  "e2MtToPFMET", 20, 0, 200 )
            self.book(f, "e1e2_DeltaPhi", "e1-e2 DeltaPhi" , 32, 0, 3.2)
            self.book(f, "e1e2_DeltaR", "e1-e2 DeltaR" , 32, 0, 3.2)


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
 

        

        ##histos[folder+'/e1Pt'].Fill( getattr(row, self.mye1+'Pt'), weight)
        ##histos[folder+'/e1Eta'].Fill(getattr(row, self.mye1+'Eta'), weight)
        ##histos[folder+'/e1Phi'].Fill(getattr(row, self.mye1+'Phi'), weight) 
        ##
        histos[folder+'/e2Pt'].Fill( row.e2Pt , weight)
        histos[folder+'/e2Eta'].Fill(row.e2Eta, weight)
        histos[folder+'/e2Phi'].Fill(row.e2Phi, weight)
        histos[folder+'/e2AbsEta'].Fill(abs(row.e2Eta), weight)
        histos[folder+'/e2Pt_vs_e2AbsEta'].Fill(abs(row.e2Eta), row.e2Pt,  weight)

        histos[folder+'/e1e2Mass'].Fill(getattr(row, 'e1_e2_Mass'), weight)
        histos[folder+'/e1MtToPFMET'].Fill(row.e1MtToPfMet,weight)
        histos[folder+'/e2MtToPFMET'].Fill(row.e2MtToPfMet,weight)
        
        histos[folder+'/e1e2_DeltaPhi'].Fill(deltaPhi(row.e1Phi, row.e2Phi), weight)
        histos[folder+'/e1e2_DeltaR'].Fill(row.e1_e2_DR, weight)


        ###histos[folder+'/type1_pfMetEt'].Fill(row.type1_pfMet_Et)
        ##histos[folder+'/ee3DR'].Fill(self.ee3DR(row)) 
        ##histos[folder+'/ee3DPhi'].Fill(self.ee3DPhi(row)) 
        ##histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight) 
        ##histos[folder+'/bjetCSVVeto30'].Fill(row.bjetCSVVeto30, weight) 
        ##
        ##histos[folder+'/ze3DR'].Fill(deltaR(self.Z(row).Phi(), getattr(row, self.mye3+'Phi'), self.Z(row).Eta(), getattr(row, self.mye3+'Eta')))
        ##histos[folder+'/ze3DPhi'].Fill(deltaPhi(self.Z(row).Phi(), getattr(row, self.mye3+'Phi')))
        ##histos[folder+'/Zpt'].Fill(self.Z(row).Pt())
            

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
            if not bool(row.e2MatchesSingleE27WP80) : continue
            if not bool(row.e1MatchesSingleE27WP80) : continue
                #if not bool(row.singleEPass) : continue
                #if not bool(row.e1MatchesSingleE) and  not bool(row.e2MatchesSingleE) : continue
                
            if row.bjetCSVVeto30!=0 : continue 
            #if jn !=0: continue
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
            if not selections.lepton_id_iso(row, 'e1', 'eid13Tight_etauiso01'): continue ##good iso point 0.5
            if abs(row.e1Eta) > 1.4442 and abs(row.e1Eta < 1.566) : continue
            
            
            cut_flow_trk.Fill('e1IDiso')
            if not selections.eSelection(row, 'e2'): continue
            cut_flow_trk.Fill('e2sel')
            if not selections.lepton_id_iso(row, 'e2', 'eid13Loose_idiso05'): continue # goodisopoint 0.5
            if abs(row.e2Eta) > 1.4442 and abs(row.e2Eta) < 1.566 : continue
            

            ##if not row.e1_e2_SS: continue
            ##            cut_flow_trk.Fill('e2IDiso')
            if  abs(row.e1_e2_Mass-91.2) < 20: continue
            ##            cut_flow_trk.Fill('ZMass')'

            if row.e1MtToPfMet < 20 : continue
            

            cut_flow_trk.Fill('tsel')


            if row.tauVetoPt20EleTight3MuLoose : continue 
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoCicLooseIso : continue # change it with Loose
            
            #if not row.e3MtToMET < 50:  continue
            cut_flow_trk.Fill('MtToMet')
            
            
            #if (row.run, row.lumi, row.evt, row.e1Pt, row.e2Pt)==myevent: continue
            #myevent=(row.run, row.lumi, row.evt, row.e1Pt, row.e2Pt)

            eleiso = 'eLoose'
            sign = 'ss' if row.e1_e2_SS else 'os'
            folder = sign+'/'+eleiso
          
            self.fill_histos(row, folder)
            folder=folder+'/'+str(int(jn))
            self.fill_histos(row, folder)
            
            if selections.lepton_id_iso(row, 'e2', 'eid13Tight_etauiso01'):
                eleiso = 'eTigh' 
                folder = sign+'/'+eleiso
                self.fill_histos(row,  folder)
                cut_flow_trk.Fill('tTightIso')
                folder=folder+'/'+str(int(jn))
                self.fill_histos(row, folder)
                
                    
 
             
        cut_flow_trk.flush()
                                
             
            
    def finish(self):
        self.write_histos()
