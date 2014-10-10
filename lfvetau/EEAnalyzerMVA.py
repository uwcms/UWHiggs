##Correction Factor still to add
from EETree import EETree
import os
import ROOT
import math
from pdb import set_trace
import glob
import array
import baseSelections as selections
import mcCorrections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, sin, cos, acos, sinh

def Z(row):
    e1p=ROOT.TVector3(row.e1Pt*cos(row.e1Phi),row.e1Pt*sin(row.e1Phi),row.e1Pt*sinh(row.e1Eta))
    e2p=ROOT.TVector3(row.e2Pt*cos(row.e2Phi),row.e2Pt*sin(row.e2Phi),row.e2Pt*sinh(row.e2Eta))
    e1FourVector= ROOT.TLorentzVector(e1p, sqrt(e1p.Mag2()+row.e1Mass*row.e1Mass))
    e2FourVector= ROOT.TLorentzVector(e2p, sqrt(e2p.Mag2()+row.e2Mass*row.e2Mass))
    zFourVector = e1FourVector+e2FourVector
    return zFourVector

def zPhi(pt1, eta1, phi1, pt2, eta2, phi2):
    px1 = pt1*cos(phi1)
    py1 = pt1*sin(phi1)
    pz1 = pt1*sinh(eta1)
    px2 = pt2*cos(phi2)
    py2 = pt2*sin(phi2)
    pz2 = pt2*sinh(eta2)
    
    px = px1+px2
    py = py1+py2
    pt = sqrt(px*px+py*py)
    phi = acos(px/pt)
    return phi

def deltaPhi(phi1, phi2):
    PHI = abs(phi1-phi2)
    if PHI<=pi:
        return PHI
    else:
        return 2*pi-PHI
def deltaR(phi1, ph2, eta1, eta2):
    deta = eta1 - eta2
    dphi = abs(phi1-phi2)
    if (dphi>pi) : dphi = 2*pi-dphi
    return sqrt(deta*deta + dphi*dphi);

class EEAnalyzerMVA(MegaBase):
    tree = 'ee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='EE'
        super(EEAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        self.tree = EETree(tree)
        self.out=outfile
        self.histograms = {}
        self.pucorrector = mcCorrections.make_puCorrector('singlee')
        self.pucorrectorUp = mcCorrections.make_puCorrectorUp('singlee')
        self.pucorrectorDown = mcCorrections.make_puCorrectorDown('singlee')

        
    def event_weight(self, row):
        if row.run > 2: #FIXME! add tight ID correction
            return [1.]
        
        etrig = 'e1'
        if row.e2Pt > row.e1Pt : etrig = 'e2'
        if bool(row.e1MatchesSingleE27WP80) and  not bool(row.e2MatchesSingleE27WP80) : etrig = 'e1'
        if not bool(row.e1MatchesSingleE27WP80) and  bool(row.e2MatchesSingleE27WP80) :  etrig = 'e2'

        allmcCorrections=    mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
                          mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                          mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                          mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
                       

        trUp_mcCorrections =   mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
                          mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                          mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                          mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_p1s_MVA(row, etrig) 
        trDown_mcCorrections = mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
                               mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                               mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                               mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_m1s_MVA(row, etrig) 

        e1idUp_mcCorrections=  mcCorrections.get_electronId_corrections13_p1s_MVA(row, 'e1') *\
                              mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                              mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                              mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
        e1idDown_mcCorrections= mcCorrections.get_electronId_corrections13_m1s_MVA(row, 'e1') * \
                                mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                                mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                                mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
        e2idUp_mcCorrections=  mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
                              mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                              mcCorrections.get_electronId_corrections13_p1s_MVA(row, 'e2') * \
                              mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
        e2idDown_mcCorrections=   mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
                                mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                                mcCorrections.get_electronId_corrections13_m1s_MVA(row, 'e2') * \
                                mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
        
        e1isoUp_mcCorrections=    mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
                               mcCorrections.get_electronIso_corrections13_p1s_MVA(row, 'e1') * \
                               mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                               mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
        e1isoDown_mcCorrections= mcCorrections.get_electronId_corrections13_m1s_MVA(row, 'e1') * \
                                 mcCorrections.get_electronIso_corrections13_p1s_MVA(row, 'e1') * \
                                 mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                                 mcCorrections.get_electronIso_corrections13_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
        e2isoUp_mcCorrections=   mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
                               mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                               mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                               mcCorrections.get_electronIso_corrections13_p1s_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
        e2isoDown_mcCorrections= mcCorrections.get_electronId_corrections13_MVA(row, 'e1') * \
                                 mcCorrections.get_electronIso_corrections13_MVA(row, 'e1') * \
                                 mcCorrections.get_electronId_corrections13_MVA(row, 'e2') * \
                                 mcCorrections.get_electronIso_corrections13_m1s_MVA(row, 'e2') * mcCorrections.get_trigger_corrections_MVA(row, etrig) 
        
    
        #pucorrlist = self.pucorrector(row.nTruePU)
        
        weight =  self.pucorrector(row.nTruePU) *\
                 allmcCorrections
        weight_up =  self.pucorrectorUp(row.nTruePU) *\
                    allmcCorrections
        weight_down =  self.pucorrectorDown(row.nTruePU) *\
                      allmcCorrections
        
        weight_tr_up = self.pucorrector(row.nTruePU) *\
                       trUp_mcCorrections
        weight_tr_down = self.pucorrector(row.nTruePU) *\
                         trDown_mcCorrections

        
        weight_e1id_up =  self.pucorrector(row.nTruePU) *\
                 e1idUp_mcCorrections
        weight_e2id_up =  self.pucorrector(row.nTruePU) *\
                 e2idUp_mcCorrections
        weight_e1id_down =  self.pucorrector(row.nTruePU) *\
                 e1idDown_mcCorrections
        weight_e2id_down =  self.pucorrector(row.nTruePU) *\
                 e2idDown_mcCorrections
        weight_e1iso_up =  self.pucorrector(row.nTruePU) *\
                 e1isoUp_mcCorrections
        weight_e2iso_up =  self.pucorrector(row.nTruePU) *\
                 e2isoUp_mcCorrections
        weight_e1iso_down =  self.pucorrector(row.nTruePU) *\
                 e1isoDown_mcCorrections
        weight_e2iso_down =  self.pucorrector(row.nTruePU) *\
                 e2isoDown_mcCorrections
        
        return  [weight, weight_up, weight_down, weight_tr_up,  weight_tr_down, weight_e1id_up, weight_e1id_down, weight_e1iso_up, weight_e1iso_down, weight_e2id_up, weight_e2id_down,  weight_e2iso_up, weight_e2iso_down]


 
## 
    def begin(self):
        sign=['os', 'ss']
        jets = [0,1,2,3]
        folder=[]
        pudir = ['','p1s/', 'm1s/','trp1s/', 'trm1s/', 'e1idp1s/','e1idm1s/',  'e1isop1s/','e1isom1s/','e2idp1s/','e2idm1s/',  'e2isop1s/','e2isom1s/']

            
        for d  in pudir :
            for i in sign:
                folder.append(d+i)
                for j in jets:
                    folder.append(d+i+'/'+str(j))
        
        #self.book('jobInfo', "jobInfo", "jobInfo", "inputfilename/C:events/l", type=pytree.PyTree)

       
        
        for f in folder: 
            #print f, f.startswith('os'), f.startswith('ss')
            if f.startswith('os') or f.startswith('ss') :
                self.book(f, "evtInfo", "evtInfo", "run/l:lumi/l:evt/l", type=pytree.PyTree)
            
            self.book(f,"e1Pt", "e1 p_{T}", 200, 0, 200)
            self.book(f,"e1Phi", "e1 phi",  100, -3.2, 3.2)
            self.book(f,"e1Eta", "e1 eta", 50, -2.5, 2.5)
            self.book(f,"e2Pt", "e2 p_{T}", 200, 0, 200)
            self.book(f,"e2Phi", "e2 phi",  100, -3.2, 3.2)
            self.book(f,"e2Eta", "e2 eta", 50, -2.5, 2.5)
            
            self.book(f, "e1e2_DeltaPhi", "e1-e2 DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e1e2_DeltaR", "e1-e2 DeltaR" , 50, 0, 3.2)
                
            self.book(f, "e1e2Mass",  "e1e2 Inv Mass",  320, 0, 320)
                
            self.book(f, "pfMetEt", "pfMetEt",  50, 0, 100)
            self.book(f, "type1_pfMetEt", "type1_pfMetEt",  50, 0, 100)
            self.book(f, "mvaMetEt", "mvaMetEt", 50, 0, 100)
            self.book(f, "pfMetPhi", "pfMetPhi", 100, -3.2, 3.2)
            self.book(f, "type1_pfMetPhi", "type1_pfMetPhi", 100, -3.2, 3.2)
            self.book(f, "mvaMetPhi", "mvaMetPhi", 100, -3.2, 3.2)

            self.book(f, "type1_pfMetEt_par", "type1_pfMetEt_par", 100, -100, 100)
            self.book(f, "type1_pfMetEt_perp", "type1_pfMetEt_perp", 50, 0, 100)
            self.book(f, "pfMetEt_par", "pfMetEt_par", 100, -100, 100)
            self.book(f, "pfMetEt_perp", "pfMetEt_perp", 50, 0, 100)
            self.book(f, "mvaMetEt_par", "mvaMetEt_par", 100, -100, 100)
            self.book(f, "mvaMetEt_perp", "mvaMetEt_perp", 50, 0, 100)
                
               
            self.book(f, "e1PFMET_DeltaPhi_Ty1", "e1-type1PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e1PFMET_Mt_Ty1", "e1-type1PFMET M_{T}" , 200, 0, 200)
            self.book(f, "e1PFMET_DeltaPhi", "e1-PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e1PFMET_Mt", "e1-PFMET M_{T}" , 200, 0, 200)
            self.book(f, "e1MVAMET_DeltaPhi", "e1-MVAMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e1MVAMET_Mt", "e1-MVAMET M_{T}" , 200, 0, 200)
            
            self.book(f, "e2PFMET_DeltaPhi_Ty1", "e2-type1PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e2PFMET_Mt_Ty1", "e2-type1PFMET M_{T}" , 200, 0, 200)
            self.book(f, "e2PFMET_DeltaPhi", "e2-PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e2PFMET_Mt", "e2-PFMET M_{T}" , 200, 0, 200)
            self.book(f, "e2MVAMET_DeltaPhi", "e2-MVAMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "e2MVAMET_Mt", "e2-MVAMET M_{T}" , 200, 0, 200)


            self.book(f,"pfMetEt_jes","pfMetEt_jes", 50, 0, 100)
            self.book(f,"pfMetPhi_jes","pfMetPhi_jes", 100, -3.2, 3.2)
            self.book(f,"pfMetEt_mes","pfMetEt_mes",50, 0, 100 )
            self.book(f,"pfMetPhi_mes","pfMetPhi_mes", 100, -3.2, 3.2)
            self.book(f,"pfMetEt_tes","pfMetEt_tes", 50, 0, 100)
            self.book(f,"pfMetPhi_tes","pfMetPhi_tes", 100, -3.2, 3.2)
            self.book(f,"pfMetEt_ees","pfMetEt_ees", 50, 0, 100)
            self.book(f,"pfMetPhi_ees","pfMetPhi_ees", 100, -3.2, 3.2)
            self.book(f,"pfMetEt_ues","pfMetEt_ues", 50, 0, 100)
            self.book(f,"pfMetPhi_ues", "pfMetPhi_ues", 100, -3.2, 3.2)

            self.book(f,"pfMetEt_perp_jes","pfMetEt_perp_jes", 50, 0, 100)
            self.book(f,"pfMetEt_perp_mes","pfMetEt_perp_mes", 50, 0, 100 )
            self.book(f,"pfMetEt_perp_tes","pfMetEt_perp_tes", 50, 0, 100)
            self.book(f,"pfMetEt_perp_ees","pfMetEt_perp_ees", 50, 0, 100)
            self.book(f,"pfMetEt_perp_ues","pfMetEt_perp_ues", 50, 0, 100)

            self.book(f,"pfMetEt_par_jes","pfMetEt_par_jes", 100, -100, 100)
            self.book(f,"pfMetEt_par_mes","pfMetEt_par_mes",100, -100, 100 )
            self.book(f,"pfMetEt_par_tes","pfMetEt_par_tes", 100, -100, 100)
            self.book(f,"pfMetEt_par_ees","pfMetEt_par_ees", 100, -100, 100)
            self.book(f,"pfMetEt_par_ues","pfMetEt_par_ues", 100, -100, 100)
            
            self.book(f,"e2PFMET_Mt_jes","e2PFMET_Mt_jes",  200, 0, 200 )
            self.book(f,"e2PFMET_Mt_mes","e2PFMET_Mt_mes",  200, 0, 200 )
            self.book(f,"e2PFMET_Mt_ees","e2PFMET_Mt_ees",  200, 0, 200 )
            self.book(f,"e2PFMET_Mt_tes","e2PFMET_Mt_tes",  200, 0, 200 )
            self.book(f,"e2PFMET_Mt_ues","e2PFMET_Mt_ues",  200, 0, 200 )

 
            self.book(f,"e1PFMET_Mt_jes","e1PFMET_Mt_jes",  200, 0, 200 )
            self.book(f,"e1PFMET_Mt_mes","e1PFMET_Mt_mes",  200, 0, 200 )
            self.book(f,"e1PFMET_Mt_ees","e1PFMET_Mt_ees",  200, 0, 200 )
            self.book(f,"e1PFMET_Mt_tes","e1PFMET_Mt_tes",  200, 0, 200 )
            self.book(f,"e1PFMET_Mt_ues","e1PFMET_Mt_ues",  200, 0, 200 )

            self.book(f,"e2PFMET_DeltaPhi_jes","e2PFMET_DeltaPhi_jes",  50, 0, 3.2 )
            self.book(f,"e2PFMET_DeltaPhi_mes","e2PFMET_DeltaPhi_mes",  50, 0, 3.2 )
            self.book(f,"e2PFMET_DeltaPhi_ees","e2PFMET_DeltaPhi_ees",  50, 0, 3.2 )
            self.book(f,"e2PFMET_DeltaPhi_tes","e2PFMET_DeltaPhi_tes",  50, 0, 3.2 )
            self.book(f,"e2PFMET_DeltaPhi_ues","e2PFMET_DeltaPhi_ues",  50, 0, 3.2 )

 

            self.book(f,"e1PFMET_DeltaPhi_jes","e1PFMET_DeltaPhi_jes",  50, 0, 3.2 )
            self.book(f,"e1PFMET_DeltaPhi_mes","e1PFMET_DeltaPhi_mes",  50, 0, 3.2 )
            self.book(f,"e1PFMET_DeltaPhi_ees","e1PFMET_DeltaPhi_ees",  50, 0, 3.2 )
            self.book(f,"e1PFMET_DeltaPhi_tes","e1PFMET_DeltaPhi_tes",  50, 0, 3.2 )
            self.book(f,"e1PFMET_DeltaPhi_ues","e1PFMET_DeltaPhi_ues",  50, 0, 3.2 )

            
            self.book(f, "nPV", "N of vertices", 100, 0, 100)
            self.book(f, "nPV_unweighted", "unweighted N of vertices", 100, 0, 100)
                    
    def fill_histos(self, row, f='os', fakeRate = False):
        
        weight = self.event_weight(row) # pu central, +1 sigma, -1sigma

        histos = self.histograms
        pudir =['']
        if row.run < 2: pudir.extend( ['p1s/', 'm1s/', 'trp1s/', 'trm1s/', 'e1idp1s/','e1idm1s/',  'e1isop1s/','e1isom1s/','e2idp1s/','e2idm1s/',  'e2isop1s/','e2isom1s/' ])
        
        for n,d  in enumerate(pudir) :
        
            folder = d+f

            histos[folder+'/e1Pt'].Fill(row.e1Pt, weight[n])
            histos[folder+'/e1Eta'].Fill(row.e1Eta, weight[n])
            histos[folder+'/e1Phi'].Fill(row.e1Phi, weight[n]) 

            histos[folder+'/e2Pt'].Fill(row.e2Pt, weight[n])
            histos[folder+'/e2Eta'].Fill(row.e2Eta, weight[n])
            histos[folder+'/e2Phi'].Fill(row.e2Phi, weight[n])

            histos[folder+'/e1e2_DeltaPhi'].Fill(deltaPhi(row.e1Phi, row.e2Phi), weight[n])
            histos[folder+'/e1e2_DeltaR'].Fill(row.e1_e2_DR, weight[n])

            histos[folder+'/e1e2Mass'].Fill(row.e1_e2_Mass, weight[n])
            
            histos[folder+'/e1PFMET_DeltaPhi'].Fill(deltaPhi(row.e1Phi, row.pfMetPhi), weight[n])
            histos[folder+'/e1MVAMET_DeltaPhi'].Fill(deltaPhi(row.e1Phi, row.mva_metPhi), weight[n])
            histos[folder+'/e1PFMET_DeltaPhi_Ty1'].Fill(deltaPhi(row.e1Phi, row.type1_pfMetPhi), weight[n])
            histos[folder+'/e1PFMET_Mt'].Fill(row.e1MtToPFMET, weight[n])
            histos[folder+'/e1MVAMET_Mt'].Fill(row.e1MtToMVAMET, weight[n])
            histos[folder+'/e1PFMET_Mt_Ty1'].Fill(row.e1MtToPfMet_Ty1, weight[n])

            histos[folder+'/type1_pfMetEt'].Fill(row.type1_pfMetEt, weight[n])
            histos[folder+'/pfMetEt'].Fill(row.pfMetEt, weight[n])
            histos[folder+'/mvaMetEt'].Fill(row.mva_metEt, weight[n])

            histos[folder+'/type1_pfMetPhi'].Fill(row.type1_pfMetPhi, weight[n])
            histos[folder+'/pfMetPhi'].Fill(row.pfMetPhi, weight[n])
            histos[folder+'/mvaMetPhi'].Fill(row.mva_metPhi, weight[n])

            zphi = Z(row).Phi()
            histos[folder+'/type1_pfMetEt_par'].Fill(row.pfMetEt*cos(deltaPhi(zphi, row.type1_pfMetPhi)), weight[n])
            histos[folder+'/type1_pfMetEt_perp'].Fill(row.pfMetEt*sin(deltaPhi(zphi, row.type1_pfMetPhi)), weight[n])
            histos[folder+'/pfMetEt_par'].Fill(row.pfMetEt*cos(deltaPhi(zphi, row.pfMetPhi)), weight[n])
            histos[folder+'/pfMetEt_perp'].Fill(row.pfMetEt*sin(deltaPhi(zphi, row.pfMetPhi)), weight[n])
            histos[folder+'/mvaMetEt_par'].Fill(row.mva_metEt*cos(deltaPhi(zphi, row.mva_metPhi)), weight[n])
            histos[folder+'/mvaMetEt_perp'].Fill(row.mva_metEt*sin(deltaPhi(zphi, row.mva_metPhi)), weight[n])
            
            histos[folder+'/e2PFMET_DeltaPhi'].Fill(deltaPhi(row.e2Phi, row.pfMetPhi), weight[n])
            histos[folder+'/e2MVAMET_DeltaPhi'].Fill(deltaPhi(row.e2Phi, row.mva_metPhi), weight[n])
            histos[folder+'/e2PFMET_DeltaPhi_Ty1'].Fill(deltaPhi(row.e2Phi, row.type1_pfMetPhi), weight[n])
            histos[folder+'/e2PFMET_Mt'].Fill(row.e2MtToPFMET, weight[n])
            histos[folder+'/e2PFMET_Mt_Ty1'].Fill(row.e2MtToPfMet_Ty1, weight[n])
            histos[folder+'/e2MVAMET_Mt'].Fill(row.e2MtToMVAMET, weight[n])

            histos[folder+'/pfMetEt_jes'].Fill(row.pfMet_jes_Et, weight[n])
            histos[folder+'/pfMetPhi_jes'].Fill(row.pfMet_jes_Phi, weight[n])
            histos[folder+'/pfMetEt_mes'].Fill(row.pfMet_mes_Et, weight[n])
            histos[folder+'/pfMetPhi_mes'].Fill(row.pfMet_mes_Phi, weight[n])
            histos[folder+'/pfMetEt_tes'].Fill(row.pfMet_tes_Et, weight[n])
            histos[folder+'/pfMetPhi_tes'].Fill(row.pfMet_tes_Phi, weight[n])
            histos[folder+'/pfMetEt_ees'].Fill(row.pfMet_ees_Et, weight[n])
            histos[folder+'/pfMetPhi_ees'].Fill(row.pfMet_ees_Phi, weight[n])
            histos[folder+'/pfMetEt_ues'].Fill(row.pfMet_ues_Et, weight[n])
            histos[folder+'/pfMetPhi_ues'].Fill(row.pfMet_ues_Phi, weight[n])
            
            histos[folder+'/pfMetEt_par_jes'].Fill(row.pfMet_jes_Et*cos(deltaPhi(zphi, row.pfMet_jes_Phi)), weight[n])
            histos[folder+'/pfMetEt_par_mes'].Fill(row.pfMet_mes_Et*cos(deltaPhi(zphi, row.pfMet_mes_Phi)), weight[n])
            histos[folder+'/pfMetEt_par_tes'].Fill(row.pfMet_tes_Et*cos(deltaPhi(zphi, row.pfMet_tes_Phi)), weight[n])
            histos[folder+'/pfMetEt_par_ees'].Fill(row.pfMet_ees_Et*cos(deltaPhi(zphi, row.pfMet_ees_Phi)), weight[n])
            histos[folder+'/pfMetEt_par_ues'].Fill(row.pfMet_ues_Et*cos(deltaPhi(zphi, row.pfMet_ues_Phi)), weight[n])

            histos[folder+'/pfMetEt_perp_jes'].Fill(row.pfMet_jes_Et*sin(deltaPhi(zphi, row.pfMet_jes_Phi)), weight[n])
            histos[folder+'/pfMetEt_perp_mes'].Fill(row.pfMet_mes_Et*sin(deltaPhi(zphi, row.pfMet_mes_Phi)), weight[n])
            histos[folder+'/pfMetEt_perp_tes'].Fill(row.pfMet_tes_Et*sin(deltaPhi(zphi, row.pfMet_tes_Phi)), weight[n])
            histos[folder+'/pfMetEt_perp_ees'].Fill(row.pfMet_ees_Et*sin(deltaPhi(zphi, row.pfMet_ees_Phi)), weight[n])
            histos[folder+'/pfMetEt_perp_ues'].Fill(row.pfMet_ues_Et*sin(deltaPhi(zphi, row.pfMet_ues_Phi)), weight[n])

  
            histos[folder+'/e2PFMET_Mt_jes'].Fill(row.e2MtToPfMet_jes, weight[n])
            histos[folder+'/e2PFMET_Mt_mes'].Fill(row.e2MtToPfMet_mes, weight[n])
            histos[folder+'/e2PFMET_Mt_ees'].Fill(row.e2MtToPfMet_ees, weight[n])
            histos[folder+'/e2PFMET_Mt_tes'].Fill(row.e2MtToPfMet_tes, weight[n])
            histos[folder+'/e2PFMET_Mt_ues'].Fill(row.e2MtToPfMet_ues, weight[n])

            histos[folder+'/e1PFMET_Mt_jes'].Fill(row.e1MtToPfMet_jes, weight[n])
            histos[folder+'/e1PFMET_Mt_mes'].Fill(row.e1MtToPfMet_mes, weight[n])
            histos[folder+'/e1PFMET_Mt_ees'].Fill(row.e1MtToPfMet_ees, weight[n])
            histos[folder+'/e1PFMET_Mt_tes'].Fill(row.e1MtToPfMet_tes, weight[n])
            histos[folder+'/e1PFMET_Mt_ues'].Fill(row.e1MtToPfMet_ues, weight[n])
            
            histos[folder+'/e2PFMET_DeltaPhi_jes'].Fill(deltaPhi(row.e2Phi, row.pfMet_jes_Phi), weight[n])
            histos[folder+'/e2PFMET_DeltaPhi_mes'].Fill(deltaPhi(row.e2Phi, row.pfMet_mes_Phi), weight[n])
            histos[folder+'/e2PFMET_DeltaPhi_ees'].Fill(deltaPhi(row.e2Phi, row.pfMet_ees_Phi), weight[n])
            histos[folder+'/e2PFMET_DeltaPhi_tes'].Fill(deltaPhi(row.e2Phi, row.pfMet_tes_Phi), weight[n])
            histos[folder+'/e2PFMET_DeltaPhi_ues'].Fill(deltaPhi(row.e2Phi, row.pfMet_ues_Phi), weight[n])

            histos[folder+'/e1PFMET_DeltaPhi_jes'].Fill(deltaPhi(row.e1Phi, row.pfMet_jes_Phi), weight[n])
            histos[folder+'/e1PFMET_DeltaPhi_mes'].Fill(deltaPhi(row.e1Phi, row.pfMet_mes_Phi), weight[n])
            histos[folder+'/e1PFMET_DeltaPhi_ees'].Fill(deltaPhi(row.e1Phi, row.pfMet_ees_Phi), weight[n])
            histos[folder+'/e1PFMET_DeltaPhi_tes'].Fill(deltaPhi(row.e1Phi, row.pfMet_tes_Phi), weight[n])
            histos[folder+'/e1PFMET_DeltaPhi_ues'].Fill(deltaPhi(row.e1Phi, row.pfMet_ues_Phi), weight[n])

## still some type1 met histos missing

            histos[folder+'/nPV'].Fill(row.nvtx, weight[n])
            histos[folder+'/nPV_unweighted'].Fill(row.nvtx)
        
         

        


    def process(self):
        event = ()
        filename = 'unnamed'
        evts_processed = 0
        #set_trace()
        for row in self.tree:
        #for i, row in enumerate(self.tree):
         #   if  i >= 100:
          #      return

            #current_file = self.tree.inputfilename
            #if filename != 'unnamed' and current_file <> filename:
            #     self.histograms['jobInfo/jobInfo'].Fill([[i for i in filename], evts_processed])
            #if not filename == 'unnamed' or current_file <> filename:
            #    filename = current_file
            evts_processed += 1
            #self.histograms[folder+'/jobInfo'].Fill(row)
 
            
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            sign = 'ss' if row.e1_e2_SS else 'os'
            #if row.run > 2 : #apply the trigger to data only (MC triggers enter in the scale factors)
            if not bool(row.singleE27WP80Pass) : continue
            if  not  bool(row.e1MatchesSingleE27WP80) and   not bool(row.e1MatchesSingleE27WP80) : continue
           
            if jn != 0 and row.bjetCSVVeto30!=0 : continue 
            
            if row.e1Pt < 30 : continue
            if row.e2Pt < 30 : continue
            if not abs(row.e1_e2_Mass-91.2) < 20: continue
            
            if not selections.eSelection(row, 'e1'): continue
            if not selections.lepton_id_iso(row, 'e1', 'eid13Tight_etauiso01'): continue
            if abs(row.e1Eta) > 1.4442 and abs(row.e1Eta) < 1.566 : continue
            
            if not selections.eSelection(row, 'e2'): continue
            if not selections.lepton_id_iso(row, 'e2', 'eid13Tight_etauiso01'): continue
            if abs(row.e2Eta) > 1.4442 and abs(row.e2Eta) < 1.566 : continue
           

            #if not selections.vetos(row) : continue
            if row.muVetoPt5IsoIdVtx : continue
            if row.eVetoCicLooseIso : continue # change it with Loose
            if row.tauVetoPt20EleTight3MuLoose : continue           
#            if row.tauHpsVetoPt20 : continue
        

            folder = sign
            self.histograms[folder+'/evtInfo'].Fill(row)
            new_event = (row.run, row.lumi, row.evt) 
            if event != new_event :
                event = new_event 
                self.fill_histos(row, folder)
                folder = sign+'/'+str(int(jn)) 
                self.fill_histos(row, folder)
             
        if filename != 'unnamed':
            self.histograms['jobInfo/jobInfo'].Fill([[i for i in filename], evts_processed])
            
    def finish(self):
        self.write_histos()


