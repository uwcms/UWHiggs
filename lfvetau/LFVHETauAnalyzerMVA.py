from ETauTree import ETauTree
import os
import ROOT
import math
import glob
import array
import mcCorrections
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, cos
from fakerate_functions import fakerate_central_histogram, fakerate_p1s_histogram, fakerate_m1s_histogram
from FinalStateAnalysis.Utilities.struct import struct



def collmass(row, met, metPhi):
    ptnu =abs(met*cos(deltaPhi(metPhi, row.tPhi)))
    visfrac = row.tPt/(row.tPt+ptnu)
    #print met, cos(deltaPhi(metPhi, row.tPhi)), ptnu, visfrac
    return (row.e_t_Mass / sqrt(visfrac))

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

class LFVHETauAnalyzerMVA(MegaBase):
    tree = 'et/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel='ET'
        super(LFVHETauAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        self.tree = ETauTree(tree)
        self.out=outfile
        self.histograms = {}
        self.pucorrector = mcCorrections.make_puCorrector('singlee')
        self.pucorrectorUp = mcCorrections.make_puCorrectorUp('singlee')
        self.pucorrectorDown = mcCorrections.make_puCorrectorDown('singlee')
        target = os.path.basename(os.environ['megatarget'])
        self.is_data = target.startswith('data_')
        self.is_embedded = ('Embedded' in target)
        self.is_mc = not (self.is_data or self.is_embedded)

    @staticmethod 
    def tau_veto(row):
        if not row.tAntiMuonLoose2 or not row.tAntiElectronMVA5Tight or not row.tDecayFinding :
            return False

    @staticmethod
    def obj1_matches_gen(row):
        return row.eGenPdgId == -1*row.eCharge*11
    @staticmethod 
    def obj3_matches_gen(row):
        return t.genDecayMode != -2 

    
    def event_weight(self, row):
        if self.is_data: #FIXME! add tight ID correction
            return [1.]

        allmcCorrections=    mcCorrections.get_electronId_corrections13_MVA(row, 'e') * \
                                 mcCorrections.get_electronIso_corrections13_MVA(row, 'e')
        trUp_mcCorrections = 1.
        trDown_mcCorrections = 1.
        eidUp_mcCorrections=  mcCorrections.get_electronId_corrections13_p1s_MVA(row, 'e') *\
                              mcCorrections.get_electronIso_corrections13_MVA(row, 'e') 
        eidDown_mcCorrections= mcCorrections.get_electronId_corrections13_m1s_MVA(row, 'e') * \
                               mcCorrections.get_electronIso_corrections13_MVA(row, 'e')
        eisoUp_mcCorrections=    mcCorrections.get_electronId_corrections13_MVA(row, 'e') * \
                                 mcCorrections.get_electronIso_corrections13_p1s_MVA(row, 'e') 
        eisoDown_mcCorrections= mcCorrections.get_electronId_corrections13_m1s_MVA(row, 'e') * \
                                mcCorrections.get_electronIso_corrections13_p1s_MVA(row, 'e') 
  

        if self.is_embedded:
            allmcCorrections= allmcCorrections*mcCorrections.get_trigger_efficiency_MVA(row,'e') 
            trUp_mcCorrections =    allmcCorrections*mcCorrections.get_trigger_efficiency_p1s_MVA(row,'e') 
            trDown_mcCorrections  = allmcCorrections*mcCorrections.get_trigger_efficiency_m1s_MVA(row,'e') 
            eidUp_mcCorrections=  eidUp_mcCorrections*  mcCorrections.get_trigger_efficiency_MVA(row,'e') 
            eidDown_mcCorrections= eidDown_mcCorrections*  mcCorrections.get_trigger_efficiency_MVA(row,'e') 
            eisoUp_mcCorrections=   eisoUp_mcCorrections * mcCorrections.get_trigger_efficiency_MVA(row,'e') 
            eisoDown_mcCorrections= eisoDown_mcCorrections* mcCorrections.get_trigger_efficiency_MVA(row,'e') 
            

            return [allmcCorrections, allmcCorrections, allmcCorrections, trUp_mcCorrections, trDown_mcCorrections ,eidUp_mcCorrections, eidDown_mcCorrections, eisoUp_mcCorrections,eisoDown_mcCorrections]


            

                       
        else:
            allmcCorrections=    allmcCorrections * mcCorrections.get_trigger_corrections_MVA(row,'e') 
            
            
            trUp_mcCorrections =    allmcCorrections*  mcCorrections.get_trigger_corrections_p1s_MVA(row,'e') 
            trDown_mcCorrections =  allmcCorrections*  mcCorrections.get_trigger_corrections_m1s_MVA(row,'e') 
            
            eidUp_mcCorrections=  eidUp_mcCorrections*  mcCorrections.get_trigger_corrections_MVA(row,'e') 
            eidDown_mcCorrections= eidDown_mcCorrections*  mcCorrections.get_trigger_corrections_MVA(row,'e') 
            eisoUp_mcCorrections=   eisoUp_mcCorrections * mcCorrections.get_trigger_corrections_MVA(row,'e') 
            eisoDown_mcCorrections= eisoDown_mcCorrections* mcCorrections.get_trigger_corrections_MVA(row,'e') 
            
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
            
            
            weight_eid_up =  self.pucorrector(row.nTruePU) *\
                             eidUp_mcCorrections
            weight_eid_down =  self.pucorrector(row.nTruePU) *\
                               eidDown_mcCorrections
            weight_eiso_up =  self.pucorrector(row.nTruePU) *\
                            eisoUp_mcCorrections
            weight_eiso_down =  self.pucorrector(row.nTruePU) *\
                                eisoDown_mcCorrections
        
            #if row.evt == 76 :
            #    print row.evt, weight, self.pucorrector(row.nTruePU) , allmcCorrections, mcCorrections.get_trigger_corrections_MVA(row,'e') , mcCorrections.get_electronId_corrections13_MVA(row, 'e'),  mcCorrections.get_electronIso_corrections13_MVA(row, 'e')
                    
            return [weight, weight_up, weight_down, weight_tr_up,  weight_tr_down, weight_eid_up, weight_eid_down, weight_eiso_up,  weight_eiso_down]


## 
    def begin(self):

        processtype=['gg']
        threshold=['ept30']
        sign=['os', 'ss']
        jetN = ['0','0_jes_plus','0_jes_minus', '1','1_jes_plus','1_jes_minus', '2','2_jes_plus','2_jes_minus', '3','3_jes_plus','3_jes_minus']
        folder=[]
        pudir = ['','p1s/', 'm1s/','trp1s/', 'trm1s/', 'eidp1s/','eidm1s/',  'eisop1s/','eisom1s/', 'mVetoUp/', 'mVetoDown/', 'eVetoUp/', 'eVetoDown/', 'tVetoUp/', 'tVetoDown/',  'tLoose/','tLooseUp/','tLooseDown/', 'tLooseUnweight/']

        for d  in pudir :
            for i in sign:
                for j in processtype:
                    for k in threshold:
                        #folder.append(d+i+'/'+j+'/'+k)
                        for jn in jetN: 

                            folder.append(d+i+'/'+j+'/'+k +'/'+jn)
                            folder.append(d+i+'/'+j+'/'+k +'/'+jn+'/selected')

            self.book(d+'os/gg/ept30/', "h_collmass_pfmet" , "h_collmass_pfmet",  32, 0, 320)
            self.book(d+'os/gg/ept30/', "h_vismass",  "h_vismass",  32, 0, 320)
            

                        
        for f in folder: 
            self.book(
                f,
                'evtInfo', 'evtInfo',
                'run/l:lumi/l:evt/l:weight/D',
                type=pytree.PyTree
            )
            
            self.book(f,"tPt", "tau p_{T}", 200, 0, 200)
            self.book(f,"tPhi", "tau phi", 100, -3.2, 3.2)
            self.book(f,"tEta", "tau eta",  50, -2.5, 2.5)
            
            self.book(f,"ePt", "e p_{T}", 200, 0, 200)
            self.book(f,"ePhi", "e phi",  100, -3.2, 3.2)
            self.book(f,"eEta", "e eta", 50, -2.5, 2.5)
            
            self.book(f, "et_DeltaPhi", "e-tau DeltaPhi" , 50, 0, 3.2)
            self.book(f, "et_DeltaR", "e-tau DeltaR" , 50, 0, 3.2)
            
            self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)
            self.book(f, "h_collmass_mvamet",  "h_collmass_mvamet",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_Ty1",  "h_collmass_pfmet_Ty1",  32, 0, 320)

            self.book(f, "h_collmass_pfmet_jes_plus", "h_collmass_pfmet_jes_plus",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_mes_plus", "h_collmass_pfmet_mes_plus",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_tes_plus", "h_collmass_pfmet_tes_plus",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_ees_plus", "h_collmass_pfmet_ees_plus",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_ues_plus", "h_collmass_pfmet_ues_plus",  32, 0, 320)

            self.book(f, "h_collmass_pfmet_jes_minus", "h_collmass_pfmet_jes_minus",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_mes_minus", "h_collmass_pfmet_mes_minus",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_tes_minus", "h_collmass_pfmet_tes_minus",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_ees_minus", "h_collmass_pfmet_ees_minus",  32, 0, 320)
            self.book(f, "h_collmass_pfmet_ues_minus", "h_collmass_pfmet_ues_minus",  32, 0, 320)

            self.book(f, "h_collmassSpread_pfmet",  "h_collmassSpread_pfmet",  40, -100, 100)
            self.book(f, "h_collmassSpread_mvamet",  "h_collmassSpread_mvamet",  40, -100, 100)
            self.book(f, "h_collmassSpread_lowPhi_pfmet",  "h_collmassSpread_lowPhi_pfmet",  40, -100, 100)
            self.book(f, "h_collmassSpread_lowPhi_mvamet",  "h_collmassSpread_lowPhi_mvamet", 40, -100, 100)
            self.book(f, "h_collmassSpread_highPhi_pfmet",  "h_collmassSpread_highPhi_pfmet", 40, -100, 100)
            self.book(f, "h_collmassSpread_highPhi_mvamet",  "h_collmassSpread_highPhi_mvamet", 40, -100, 100)
            self.book(f, "h_collmass_lowPhi_pfmet",  "h_collmass_lowPhi_pfmet",  32, 0, 320)
            self.book(f, "h_collmass_lowPhi_mvamet",  "h_collmass_lowPhi_mvamet",  32, 0, 320)
            self.book(f, "h_collmass_highPhi_pfmet",  "h_collmass_highPhi_pfmet",  32, 0, 320)
            self.book(f, "h_collmass_highPhi_mvamet", "h_collmass_highPhi_mvamet",  32, 0, 320)
            self.book(f, "h_collmass_vs_dPhi_pfmet",  "h_collmass_vs_dPhi_pfmet", 50, 0, 3.2, 32, 0, 320, type=ROOT.TH2F)
            self.book(f, "h_collmass_vs_dPhi_mvamet",  "h_collmass_vs_dPhi_mvamet", 50, 0, 3.2, 32, 0, 320, type=ROOT.TH2F)
            self.book(f, "h_collmassSpread_vs_dPhi_pfmet",  "h_collmassSpread_vs_dPhi_pfmet", 50, 0, 3.2, 20, -100, 100, type=ROOT.TH2F)
            self.book(f, "h_collmassSpread_vs_dPhi_mvamet",  "h_collmassSpread_vs_dPhi_mvamet", 50, 0, 3.2, 20, -100, 100, type=ROOT.TH2F)

                
            
            self.book(f, "h_vismass",  "h_vismass",  32, 0, 320)
            
            self.book(f, "type1_pfMet_Et", "PFMet", 200, 0, 200)
            self.book(f, "pfMet_Et", "PFMet", 200, 0, 200)
            self.book(f, "type1_pfMet_Phi", "PFMet #phi", 100, -3.2, 3.2)
            self.book(f, "pfMet_Phi", "PFMet #phi", 100, -3.2, 3.2)
            
            self.book(f, "pfMet_Et_ees_minus",  "pfMet_Et_ees_minus",  200, 0, 200)
            self.book(f, "pfMet_Et_jes_minus",  "pfMet_Et_jes_minus",  200, 0, 200)
            self.book(f, "pfMet_Et_mes_minus",  "pfMet_Et_mes_minus",  200, 0, 200)
            self.book(f, "pfMet_Et_tes_minus",  "pfMet_Et_tes_minus",  200, 0, 200)
            self.book(f, "pfMet_Et_ues_minus",  "pfMet_Et_ues_minus",  200, 0, 200)


            self.book(f, "pfMet_Et_jes_plus",   "pfMet_Et_jes_plus",   200, 0, 200)
            self.book(f, "pfMet_Et_ees_plus",   "pfMet_Et_ees_plus",   200, 0, 200)
            self.book(f, "pfMet_Et_mes_plus",   "pfMet_Et_mes_plus",   200, 0, 200)
            self.book(f, "pfMet_Et_tes_plus",   "pfMet_Et_tes_plus",   200, 0, 200)
            self.book(f, "pfMet_Et_ues_plus",   "pfMet_Et_ues_plus",   200, 0, 200)

            self.book(f, "pfMet_Phi_ees_plus",  "pfMet_Phi_ees_plus",  100, -3.2, 3.2)
            self.book(f, "pfMet_Phi_jes_plus",  "pfMet_Phi_jes_plus",  100, -3.2, 3.2)
            self.book(f, "pfMet_Phi_mes_plus",  "pfMet_Phi_mes_plus",  100, -3.2, 3.2)
            self.book(f, "pfMet_Phi_tes_plus",  "pfMet_Phi_tes_plus",  100, -3.2, 3.2)
            self.book(f, "pfMet_Phi_ues_plus",  "pfMet_Phi_ues_plus",  100, -3.2, 3.2)

            self.book(f, "pfMet_Phi_ees_minus", "pfMet_Phi_ees_minus", 100, -3.2, 3.2)
            self.book(f, "pfMet_Phi_jes_minus", "pfMet_Phi_jes_minus", 100, -3.2, 3.2)
            self.book(f, "pfMet_Phi_mes_minus", "pfMet_Phi_mes_minus", 100, -3.2, 3.2)
            self.book(f, "pfMet_Phi_tes_minus", "pfMet_Phi_tes_minus", 100, -3.2, 3.2)
            self.book(f, "pfMet_Phi_ues_minus", "pfMet_Phi_ues_minus", 100, -3.2, 3.2)

            #self.book(f, "pfMet_jes_Et",        "pfMet_jes_Et",        200, 0, 200)
            #self.book(f, "pfMet_jes_Phi",       "pfMet_jes_Phi",       100, -3.2, 3.2)
            #self.book(f, "pfMet_ues_AtanToPhi", "pfMet_ues_AtanToPhi", 100, -3.2, 3.2)

 
            self.book(f, "type1_pfMet_Et_vs_dPhi", "PFMet vs #Delta#phi(#tau,PFMet)", 50, 0, 3.2, 64, 0, 320, type=ROOT.TH2F)
            self.book(f, "mvaMet_Et_vs_dPhi", "MVAMet vs #Delta#phi(#tau,MVAMet)", 50, 0, 3.2, 64, 0, 320, type=ROOT.TH2F)

            self.book(f, "tPFMET_DeltaPhi", "tau-PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "tPFMET_Mt", "tau-PFMET M_{T}" , 200, 0, 200)
            self.book(f, "tPFMET_DeltaPhi_Ty1", "tau-type1PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "tPFMET_Mt_Ty1", "tau-type1PFMET M_{T}" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_jes_plus', "tau-MVAMET M_{T} JES_plus" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_mes_plus', "tau-MVAMET M_{T} JES_plus" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_ees_plus', "tau-MVAMET M_{T} JES_plus" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_tes_plus', "tau-MVAMET M_{T} JES_plus" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_ues_plus', "tau-MVAMET M_{T} JES_plus" , 200, 0, 200)

            self.book(f, 'tPFMET_Mt_jes_minus', "tau-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_mes_minus', "tau-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_ees_minus', "tau-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_tes_minus', "tau-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'tPFMET_Mt_ues_minus', "tau-MVAMET M_{T} JES_minus" , 200, 0, 200)
            
            self.book(f, "tMVAMET_DeltaPhi", "tau-MVAMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "tMVAMET_Mt", "tau-MVAMET M_{T}" , 200, 0, 200)
               
            self.book(f, "ePFMET_DeltaPhi_Ty1", "e-type1PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "ePFMET_Mt_Ty1", "e-type1PFMET M_{T}" , 200, 0, 200)
            self.book(f, "ePFMET_DeltaPhi", "e-PFMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "ePFMET_Mt", "e-PFMET M_{T}" , 200, 0, 200)
            #self.book(f, 'ePFMET_Mt_jes', "e-MVAMET M_{T} JES" , 200, 0, 200)
            #self.book(f, 'ePFMET_Mt_mes', "e-MVAMET M_{T} JES" , 200, 0, 200)
            #self.book(f, 'ePFMET_Mt_ees', "e-MVAMET M_{T} JES" , 200, 0, 200)
            #self.book(f, 'ePFMET_Mt_tes', "e-MVAMET M_{T} JES" , 200, 0, 200)
            #self.book(f, 'ePFMET_Mt_ues', "e-MVAMET M_{T} JES" , 200, 0, 200)
            #self.book(f, "ePFMET_Mt_Ty1_ues_minus", "e-type1PFMET M_{T} ues_minus" , 200, 0, 200)
            #self.book(f, "ePFMET_Mt_Ty1_ues_plus", "e-type1PFMET M_{T} ues_plus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_jes_minus', "e-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_mes_minus', "e-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_ees_minus', "e-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_tes_minus', "e-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_ues_minus', "e-MVAMET M_{T} JES_minus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_jes_plus', "e-MVAMET M_{T} JES_plus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_mes_plus', "e-MVAMET M_{T} JES_plus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_ees_plus', "e-MVAMET M_{T} JES_plus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_tes_plus', "e-MVAMET M_{T} JES_plus" , 200, 0, 200)
            self.book(f, 'ePFMET_Mt_ues_plus', "e-MVAMET M_{T} JES_plus" , 200, 0, 200)

            self.book(f, "eMVAMET_DeltaPhi", "e-MVAMET DeltaPhi" , 50, 0, 3.2)
            self.book(f, "eMVAMET_Mt", "e-MVAMET M_{T}" , 200, 0, 200)
            
            self.book(f, "jetN_20", "Number of jets, p_{T}>20", 10, -0.5, 9.5) 
            self.book(f, "jetN_30", "Number of jets, p_{T}>30", 10, -0.5, 9.5) 

    def fakerate_weights(self, tEta, central_weights, p1s_weights, m1s_weights):
        frweight=[1.,1.,1.]

        #central_weights = fakerate_central_histogram(25,0, 2.5)
        #p1s_weights = fakerate_central_histogram(25,0, 2.5)
        #m1s_weights = fakerate_central_histogram(25,0, 2.5)

        for n,w in enumerate( central_weights ):
            if abs(tEta) < w[1]:
                break
            frweight[0] = w[0]
            frweight[1] = p1s_weights[n][0]
            frweight[2] = m1s_weights[n][0]
 
        
        return  frweight;

    
                    
    def fill_histos(self, row, f='os/gg/ept0/0',  isTauTight=False, frw=[1.,1.,1.]):
        weight = self.event_weight(row)
        histos = self.histograms
        pudir =['']
        if self.is_data == False : pudir.extend( ['p1s/', 'm1s/', 'trp1s/', 'trm1s/', 'eidp1s/','eidm1s/',  'eisop1s/','eisom1s/'])
        looseList = ['tLoose/', 'tLooseUp/', 'tLooseDown/', 'tLooseUnweight/']
        
        
        if bool(isTauTight) == False:               
            if f.startswith('os') or  f.startswith('ss')  :
                frweight_bv = frw[0]/(1.-frw[0])
                #err = 0.3*abs(2-frw[0]/pow(1-frw[0], 2)) #tau pog told to mu-tau group to use 30% uncertainty on tau fake rate.
                err=0.3
                frweight_p1s = frweight_bv*(1+err)
                frweight_m1s = frweight_bv*(1-err)
                fr_weights = [frweight_bv, frweight_p1s, frweight_m1s]
            
                for n, l in enumerate(looseList) :
                    frweight = weight[0]*fr_weights[n] if n < len(looseList)-1  else weight[0]
                    folder = l+f
                    frweight = row.EmbPtWeight*frweight
                    if f=='os/gg/ept30' :
                        histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), frweight)
                        histos[folder+'/h_vismass'].Fill(row.e_t_Mass, frweight)
                        continue

                    
                    histos[folder+'/tPt'].Fill(row.tPt, frweight)
                    histos[folder+'/tEta'].Fill(row.tEta, frweight)
                    histos[folder+'/tPhi'].Fill(row.tPhi, frweight) 
                    histos[folder+'/ePt'].Fill(row.ePt, frweight)
                    histos[folder+'/eEta'].Fill(row.eEta, frweight)
                    histos[folder+'/ePhi'].Fill(row.ePhi, frweight)
                    histos[folder+'/et_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.tPhi), frweight)
                    histos[folder+'/et_DeltaR'].Fill(row.e_t_DR, frweight)
                    histos[folder+'/h_collmass_vs_dPhi_pfmet'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), frweight)
                    histos[folder+'/h_collmass_vs_dPhi_mvamet'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), collmass(row, row.mva_metEt, row.mva_metPhi), frweight)
                    histos[folder+'/h_collmassSpread_vs_dPhi_pfmet'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, frweight)
                    histos[folder+'/h_collmassSpread_vs_dPhi_mvamet'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, frweight)
                    if deltaPhi(row.tPhi, row.pfMetPhi) > 1.57 :  
                        histos[folder+'/h_collmass_highPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), frweight)
                        histos[folder+'/h_collmassSpread_highPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, frweight)
                    if deltaPhi(row.tPhi, row.pfMetPhi) < 1.57 :  
                        histos[folder+'/h_collmass_lowPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), frweight)
                        histos[folder+'/h_collmassSpread_lowPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, frweight)
                    if deltaPhi(row.tPhi, row.mva_metPhi) > 1.57 :  
                        histos[folder+'/h_collmass_highPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), frweight)
                        histos[folder+'/h_collmassSpread_highPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, frweight)
                    if deltaPhi(row.tPhi, row.mva_metPhi) < 1.57 :  
                        histos[folder+'/h_collmass_lowPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), frweight)
                        histos[folder+'/h_collmassSpread_lowPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, frweight)

                    histos[folder+'/h_collmassSpread_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, frweight)
                    histos[folder+'/h_collmassSpread_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, frweight)
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.pfMetEt, row.pfMetPhi), frweight)
                    histos[folder+'/h_collmass_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), frweight)
                    histos[folder+'/h_collmass_pfmet_Ty1'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), frweight)

                    #histos[folder+'/h_collmass_pfmet_jes'].Fill(collmass(row, row.pfMet_jes_Et, row.pfMet_jes_Phi), frweight)
                    #histos[folder+'/h_collmass_pfmet_mes'].Fill(collmass(row, row.pfMet_mes_Et, row.pfMet_mes_Phi), frweight)
                    #histos[folder+'/h_collmass_pfmet_tes'].Fill(collmass(row, row.pfMet_tes_Et, row.pfMet_tes_Phi), frweight)
                    #histos[folder+'/h_collmass_pfmet_ees'].Fill(collmass(row, row.pfMet_ees_Et, row.pfMet_ees_Phi), frweight)
                    #histos[folder+'/h_collmass_pfmet_ues'].Fill(collmass(row, row.pfMet_ues_Et, row.pfMet_ues_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_ees_minus'].Fill(collmass(row, row.pfMet_ees_minus_Et, row.pfMet_ees_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_ees_plus' ].Fill(collmass(row, row.pfMet_ees_plus_Et , row.pfMet_ees_plus_Phi ), frweight)
                    histos[folder+'/h_collmass_pfmet_jes_minus'].Fill(collmass(row, row.pfMet_jes_minus_Et, row.pfMet_jes_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_jes_plus' ].Fill(collmass(row, row.pfMet_jes_plus_Et , row.pfMet_jes_plus_Phi ), frweight)
                    histos[folder+'/h_collmass_pfmet_mes_minus'].Fill(collmass(row, row.pfMet_mes_minus_Et, row.pfMet_mes_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_mes_plus' ].Fill(collmass(row, row.pfMet_mes_plus_Et , row.pfMet_mes_plus_Phi ), frweight)
                    histos[folder+'/h_collmass_pfmet_tes_minus'].Fill(collmass(row, row.pfMet_tes_minus_Et, row.pfMet_tes_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_tes_plus' ].Fill(collmass(row, row.pfMet_tes_plus_Et , row.pfMet_tes_plus_Phi ), frweight)
                    histos[folder+'/h_collmass_pfmet_ues_minus'].Fill(collmass(row, row.pfMet_ues_minus_Et, row.pfMet_ues_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_ues_plus' ].Fill(collmass(row, row.pfMet_ues_plus_Et , row.pfMet_ues_plus_Phi ), frweight)

                    
                    histos[folder+'/h_vismass'].Fill(row.e_t_Mass, frweight)
                    
                    histos[folder+'/type1_pfMet_Et'].Fill(row.type1_pfMetEt, frweight)
                    histos[folder+'/pfMet_Et'].Fill(row.pfMetEt, frweight)
                    histos[folder+'/type1_pfMet_Phi'].Fill(row.type1_pfMetPhi, frweight)
                    histos[folder+'/pfMet_Phi'].Fill(row.pfMetPhi, frweight)

                    histos[folder+'/pfMet_Et_ees_minus'].Fill( row.pfMet_ees_minus_Et ,  frweight)
                    histos[folder+'/pfMet_Et_jes_minus'].Fill( row.pfMet_jes_minus_Et ,  frweight)
                    histos[folder+'/pfMet_Et_mes_minus'].Fill( row.pfMet_mes_minus_Et ,  frweight)
                    histos[folder+'/pfMet_Et_tes_minus'].Fill( row.pfMet_tes_minus_Et ,  frweight)
                    histos[folder+'/pfMet_Et_ues_minus'].Fill( row.pfMet_ues_minus_Et ,  frweight)

                    histos[folder+'/pfMet_Et_ees_plus'].Fill(  row.pfMet_ees_plus_Et  ,  frweight)
                    histos[folder+'/pfMet_Et_jes_plus'].Fill(  row.pfMet_jes_plus_Et  ,  frweight)
                    histos[folder+'/pfMet_Et_mes_plus'].Fill(  row.pfMet_mes_plus_Et  ,  frweight)
                    histos[folder+'/pfMet_Et_tes_plus'].Fill(  row.pfMet_tes_plus_Et  ,  frweight)
                    histos[folder+'/pfMet_Et_ues_plus'].Fill(  row.pfMet_ues_plus_Et  ,  frweight)
                    
                    
                    histos[folder+'/pfMet_Phi_ees_minus'].Fill(row.pfMet_ees_minus_Phi,  frweight)
                    histos[folder+'/pfMet_Phi_jes_minus'].Fill(row.pfMet_jes_minus_Phi,  frweight)
                    histos[folder+'/pfMet_Phi_mes_minus'].Fill(row.pfMet_mes_minus_Phi,  frweight)
                    histos[folder+'/pfMet_Phi_tes_minus'].Fill(row.pfMet_tes_minus_Phi,  frweight)
                    histos[folder+'/pfMet_Phi_ues_minus'].Fill(row.pfMet_ues_minus_Phi,  frweight)

                    histos[folder+'/pfMet_Phi_ees_plus'].Fill( row.pfMet_ees_plus_Phi ,  frweight)
                    histos[folder+'/pfMet_Phi_jes_plus'].Fill( row.pfMet_jes_plus_Phi ,  frweight)
                    histos[folder+'/pfMet_Phi_mes_plus'].Fill( row.pfMet_mes_plus_Phi ,  frweight)
                    histos[folder+'/pfMet_Phi_tes_plus'].Fill( row.pfMet_tes_plus_Phi ,  frweight)
                    histos[folder+'/pfMet_Phi_ues_plus'].Fill( row.pfMet_ues_plus_Phi ,  frweight)

                    #histos[folder+'/pfMet_jes_Et'].Fill(       row.pfMet_jes_Et       ,  frweight)
                    #histos[folder+'/pfMet_jes_Phi'].Fill(      row.pfMet_jes_Phi      ,  frweight)

                    #histos[folder+'/pfMet_ues_AtanToPhi'].Fill(rowpfMet_ues_AtanToPhi, frweight)
 
                    histos[folder+'/type1_pfMet_Et_vs_dPhi'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), row.type1_pfMetEt, frweight)
                    histos[folder+'/mvaMet_Et_vs_dPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), row.mva_metEt, frweight)
                        
                    histos[folder+'/ePFMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.pfMetPhi), frweight)
                    histos[folder+'/ePFMET_DeltaPhi_Ty1'].Fill(deltaPhi(row.ePhi, row.type1_pfMetPhi), frweight)
                    histos[folder+'/eMVAMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mva_metPhi), frweight)
                    histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPFMET, frweight)
                    histos[folder+'/ePFMET_Mt_Ty1'].Fill(row.eMtToPfMet_Ty1, frweight)
                    histos[folder+'/eMVAMET_Mt'].Fill(row.eMtToMVAMET, frweight)
                    
                        
                    histos[folder+'/tPFMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.pfMetPhi), frweight)
                    histos[folder+'/tPFMET_DeltaPhi_Ty1'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), frweight)
                    histos[folder+'/tMVAMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), frweight)
                    histos[folder+'/tPFMET_Mt'].Fill(row.tMtToPFMET, frweight)
                    histos[folder+'/tMVAMET_Mt'].Fill(row.tMtToMVAMET, frweight)
                    
                    histos[folder+'/h_collmass_pfmet_ees_minus'].Fill(collmass(row, row.pfMet_ees_minus_Et, row.pfMet_ees_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_ees_plus' ].Fill(collmass(row, row.pfMet_ees_plus_Et , row.pfMet_ees_plus_Phi ), frweight)
                    histos[folder+'/h_collmass_pfmet_jes_minus'].Fill(collmass(row, row.pfMet_jes_minus_Et, row.pfMet_jes_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_jes_plus' ].Fill(collmass(row, row.pfMet_jes_plus_Et , row.pfMet_jes_plus_Phi ), frweight)
                    histos[folder+'/h_collmass_pfmet_mes_minus'].Fill(collmass(row, row.pfMet_mes_minus_Et, row.pfMet_mes_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_mes_plus' ].Fill(collmass(row, row.pfMet_mes_plus_Et , row.pfMet_mes_plus_Phi ), frweight)
                    histos[folder+'/h_collmass_pfmet_tes_minus'].Fill(collmass(row, row.pfMet_tes_minus_Et, row.pfMet_tes_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_tes_plus' ].Fill(collmass(row, row.pfMet_tes_plus_Et , row.pfMet_tes_plus_Phi ), frweight)
                    histos[folder+'/h_collmass_pfmet_ues_minus'].Fill(collmass(row, row.pfMet_ues_minus_Et, row.pfMet_ues_minus_Phi), frweight)
                    histos[folder+'/h_collmass_pfmet_ues_plus' ].Fill(collmass(row, row.pfMet_ues_plus_Et , row.pfMet_ues_plus_Phi ), frweight)

                                        
                    histos[folder+'/ePFMET_Mt_jes_minus'].Fill(row.eMtToPfMet_jes_minus,frweight)
                    histos[folder+'/ePFMET_Mt_jes_plus'].Fill(row.eMtToPfMet_jes_plus,  frweight)
                    
                    histos[folder+'/ePFMET_Mt_mes_minus'].Fill(row.eMtToPfMet_mes_minus,frweight) 
                    histos[folder+'/ePFMET_Mt_mes_plus'].Fill(row.eMtToPfMet_mes_plus,  frweight) 
                    
                    histos[folder+'/ePFMET_Mt_ees_minus'].Fill(row.eMtToPfMet_ees_minus,frweight) 
                    histos[folder+'/ePFMET_Mt_ees_plus'].Fill(row.eMtToPfMet_ees_plus,  frweight) 
                    
                    histos[folder+'/ePFMET_Mt_tes_minus'].Fill(row.eMtToPfMet_tes_minus,frweight) 
                    histos[folder+'/ePFMET_Mt_tes_plus'].Fill(row.eMtToPfMet_tes_plus,  frweight) 
                
                    histos[folder+'/ePFMET_Mt_ues_minus'].Fill(row.eMtToPfMet_ues_minus,frweight) 
                    histos[folder+'/ePFMET_Mt_ues_plus'].Fill(row.eMtToPfMet_ues_plus,  frweight) 
                    
                    histos[folder+'/tPFMET_Mt_jes_plus'].Fill(row.tMtToPfMet_jes_plus,  frweight)
                    histos[folder+'/tPFMET_Mt_mes_plus'].Fill(row.tMtToPfMet_mes_plus,  frweight)
                    histos[folder+'/tPFMET_Mt_ees_plus'].Fill(row.tMtToPfMet_ees_plus,  frweight)
                    histos[folder+'/tPFMET_Mt_tes_plus'].Fill(row.tMtToPfMet_tes_plus,  frweight)
                    histos[folder+'/tPFMET_Mt_ues_plus'].Fill(row.tMtToPfMet_ues_plus,  frweight)
                    histos[folder+'/tPFMET_Mt_jes_minus'].Fill(row.tMtToPfMet_jes_minus,frweight)
                    histos[folder+'/tPFMET_Mt_mes_minus'].Fill(row.tMtToPfMet_mes_minus,frweight)
                    histos[folder+'/tPFMET_Mt_ees_minus'].Fill(row.tMtToPfMet_ees_minus,frweight)
                    histos[folder+'/tPFMET_Mt_tes_minus'].Fill(row.tMtToPfMet_tes_minus,frweight)
                    histos[folder+'/tPFMET_Mt_ues_minus'].Fill(row.tMtToPfMet_ues_minus,frweight)
                


                    histos[folder+'/jetN_20'].Fill(row.jetVeto20, frweight) 
                    histos[folder+'/jetN_30'].Fill(row.jetVeto30, frweight) 
                    
        else: # if it is TauTight
            if not f.startswith('os') and not  f.startswith('ss') : # if the dir name start with mVeto, eVeto or tVeto I don't want the different weight of MC corrections
                pudir = ['']
                    
            for n,d  in enumerate(pudir) :
                if 'minus' in f  or 'plus' in f :
                    if n >0 : break
                    
                folder = d+f
                weight[n] = row.EmbPtWeight*weight[n]
                if f=='os/gg/ept30' :
                    histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                    histos[folder+'/h_vismass'].Fill(row.e_t_Mass, weight[n])
                    continue

                histos[folder+'/evtInfo'].Fill( struct(run=row.run,lumi=row.lumi,evt=row.evt,weight=weight[n]))

                histos[folder+'/tPt'].Fill(row.tPt, weight[n])
                histos[folder+'/tEta'].Fill(row.tEta, weight[n])
                histos[folder+'/tPhi'].Fill(row.tPhi, weight[n]) 
                histos[folder+'/ePt'].Fill(row.ePt, weight[n])
                histos[folder+'/eEta'].Fill(row.eEta, weight[n])
                histos[folder+'/ePhi'].Fill(row.ePhi, weight[n])
                histos[folder+'/et_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.tPhi), weight[n])
                histos[folder+'/et_DeltaR'].Fill(row.e_t_DR, weight[n])
                    
                histos[folder+'/h_collmass_vs_dPhi_pfmet'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                histos[folder+'/h_collmass_vs_dPhi_mvamet'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), collmass(row, row.mva_metEt, row.mva_metPhi), weight[n])
                histos[folder+'/h_collmassSpread_vs_dPhi_pfmet'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, weight[n])
                histos[folder+'/h_collmassSpread_vs_dPhi_mvamet'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, weight[n])
                if deltaPhi(row.tPhi, row.pfMetPhi) > 1.57 :  
                    histos[folder+'/h_collmass_highPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                    histos[folder+'/h_collmassSpread_highPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, weight[n])
                if deltaPhi(row.tPhi, row.pfMetPhi) < 1.57 :  
                    histos[folder+'/h_collmass_lowPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                    histos[folder+'/h_collmassSpread_lowPhi_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, weight[n])
                if deltaPhi(row.tPhi, row.mva_metPhi) > 1.57 :  
                    histos[folder+'/h_collmass_highPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), weight[n])
                    histos[folder+'/h_collmassSpread_highPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, weight[n])
                if deltaPhi(row.tPhi, row.mva_metPhi) < 1.57 :  
                    histos[folder+'/h_collmass_lowPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), weight[n])
                    histos[folder+'/h_collmassSpread_lowPhi_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, weight[n])
                histos[folder+'/h_collmassSpread_pfmet'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi)-125.0, weight[n])
                histos[folder+'/h_collmassSpread_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi)-125.0, weight[n])


                histos[folder+'/h_collmass_pfmet'].Fill(collmass(row, row.pfMetEt, row.pfMetPhi), weight[n])
                histos[folder+'/h_collmass_mvamet'].Fill(collmass(row, row.mva_metEt, row.mva_metPhi), weight[n])
                histos[folder+'/h_collmass_pfmet_Ty1'].Fill(collmass(row, row.type1_pfMetEt, row.type1_pfMetPhi), weight[n])
                
                 
                histos[folder+'/h_vismass'].Fill(row.e_t_Mass, weight[n])

                histos[folder+'/type1_pfMet_Et'].Fill(row.type1_pfMetEt, weight[n])
                histos[folder+'/pfMet_Et'].Fill(row.pfMetEt, weight[n])
                histos[folder+'/type1_pfMet_Phi'].Fill(row.type1_pfMetPhi, weight[n])
                histos[folder+'/pfMet_Phi'].Fill(row.pfMetPhi,weight[n])
                
                histos[folder+'/pfMet_Et_ees_minus'].Fill( row.pfMet_ees_minus_Et ,  weight[n])
                histos[folder+'/pfMet_Et_jes_minus'].Fill( row.pfMet_jes_minus_Et ,  weight[n])
                histos[folder+'/pfMet_Et_mes_minus'].Fill( row.pfMet_mes_minus_Et ,  weight[n])
                histos[folder+'/pfMet_Et_tes_minus'].Fill( row.pfMet_tes_minus_Et ,  weight[n])
                histos[folder+'/pfMet_Et_ues_minus'].Fill( row.pfMet_ues_minus_Et ,  weight[n])
                
                histos[folder+'/pfMet_Et_ees_plus'].Fill(  row.pfMet_ees_plus_Et  ,  weight[n])
                histos[folder+'/pfMet_Et_jes_plus'].Fill(  row.pfMet_jes_plus_Et  ,  weight[n])
                histos[folder+'/pfMet_Et_mes_plus'].Fill(  row.pfMet_mes_plus_Et  ,  weight[n])
                histos[folder+'/pfMet_Et_tes_plus'].Fill(  row.pfMet_tes_plus_Et  ,  weight[n])
                histos[folder+'/pfMet_Et_ues_plus'].Fill(  row.pfMet_ues_plus_Et  ,  weight[n])
                    
                
                histos[folder+'/pfMet_Phi_ees_minus'].Fill(row.pfMet_ees_minus_Phi,  weight[n])
                histos[folder+'/pfMet_Phi_jes_minus'].Fill(row.pfMet_jes_minus_Phi,  weight[n])
                histos[folder+'/pfMet_Phi_mes_minus'].Fill(row.pfMet_mes_minus_Phi,  weight[n])
                histos[folder+'/pfMet_Phi_tes_minus'].Fill(row.pfMet_tes_minus_Phi,  weight[n])
                histos[folder+'/pfMet_Phi_ues_minus'].Fill(row.pfMet_ues_minus_Phi,  weight[n])

                histos[folder+'/pfMet_Phi_ees_plus'].Fill( row.pfMet_ees_plus_Phi ,  weight[n])
                histos[folder+'/pfMet_Phi_jes_plus'].Fill( row.pfMet_jes_plus_Phi ,  weight[n])
                histos[folder+'/pfMet_Phi_mes_plus'].Fill( row.pfMet_mes_plus_Phi ,  weight[n])
                histos[folder+'/pfMet_Phi_tes_plus'].Fill( row.pfMet_tes_plus_Phi ,  weight[n])
                histos[folder+'/pfMet_Phi_ues_plus'].Fill( row.pfMet_ues_plus_Phi ,  weight[n])

                #histos[folder+'/pfMet_jes_Et'].Fill(       row.pfMet_jes_Et       , weight[n])
                #histos[folder+'/pfMet_jes_Phi'].Fill(      row.pfMet_jes_Phi      , weight[n])

                #histos[folder+'/pfMet_ues_AtanToPhi'].Fill(row.pfMet_ues_AtanToPhi, weight[n])


                histos[folder+'/type1_pfMet_Et_vs_dPhi'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), row.type1_pfMetEt, weight[n])
                histos[folder+'/mvaMet_Et_vs_dPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), row.mva_metEt, weight[n])
                
                histos[folder+'/ePFMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.type1_pfMetPhi), weight[n])
                histos[folder+'/ePFMET_DeltaPhi_Ty1'].Fill(deltaPhi(row.ePhi, row.type1_pfMetPhi), weight[n])
                histos[folder+'/eMVAMET_DeltaPhi'].Fill(deltaPhi(row.ePhi, row.mva_metPhi), weight[n])
                
                histos[folder+'/ePFMET_Mt'].Fill(row.eMtToPFMET, weight[n])
                histos[folder+'/ePFMET_Mt_Ty1'].Fill(row.eMtToPfMet_Ty1, weight[n])
                histos[folder+'/eMVAMET_Mt'].Fill(row.eMtToMVAMET, weight[n])
                
                histos[folder+'/tPFMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.pfMetPhi), weight[n])
                histos[folder+'/tPFMET_DeltaPhi_Ty1'].Fill(deltaPhi(row.tPhi, row.type1_pfMetPhi), weight[n])
                histos[folder+'/tMVAMET_DeltaPhi'].Fill(deltaPhi(row.tPhi, row.mva_metPhi), weight[n])
                
                histos[folder+'/tPFMET_Mt'].Fill(row.tMtToPFMET, weight[n])
                histos[folder+'/tPFMET_Mt_Ty1'].Fill(row.tMtToPfMet_Ty1, weight[n])
                histos[folder+'/tMVAMET_Mt'].Fill(row.tMtToMVAMET, weight[n])
                    
                histos[folder+'/jetN_20'].Fill(row.jetVeto20, weight[n]) 
                histos[folder+'/jetN_30'].Fill(row.jetVeto30, weight[n]) 
                

                #                if n == 0: # if I'm in the dir starting with os or ss I want also the energy scale uncertainties
                histos[folder+'/h_collmass_pfmet_ees_minus'].Fill(collmass(row, row.pfMet_ees_minus_Et, row.pfMet_ees_minus_Phi), weight[n])
                histos[folder+'/h_collmass_pfmet_ees_plus' ].Fill(collmass(row, row.pfMet_ees_plus_Et , row.pfMet_ees_plus_Phi ), weight[n])
                histos[folder+'/h_collmass_pfmet_jes_minus'].Fill(collmass(row, row.pfMet_jes_minus_Et, row.pfMet_jes_minus_Phi), weight[n])
                histos[folder+'/h_collmass_pfmet_jes_plus' ].Fill(collmass(row, row.pfMet_jes_plus_Et , row.pfMet_jes_plus_Phi ), weight[n])
                histos[folder+'/h_collmass_pfmet_mes_minus'].Fill(collmass(row, row.pfMet_mes_minus_Et, row.pfMet_mes_minus_Phi), weight[n])
                histos[folder+'/h_collmass_pfmet_mes_plus' ].Fill(collmass(row, row.pfMet_mes_plus_Et , row.pfMet_mes_plus_Phi ), weight[n])
                histos[folder+'/h_collmass_pfmet_tes_minus'].Fill(collmass(row, row.pfMet_tes_minus_Et, row.pfMet_tes_minus_Phi), weight[n])
                histos[folder+'/h_collmass_pfmet_tes_plus' ].Fill(collmass(row, row.pfMet_tes_plus_Et , row.pfMet_tes_plus_Phi ), weight[n])
                histos[folder+'/h_collmass_pfmet_ues_minus'].Fill(collmass(row, row.pfMet_ues_minus_Et, row.pfMet_ues_minus_Phi), weight[n])
                histos[folder+'/h_collmass_pfmet_ues_plus' ].Fill(collmass(row, row.pfMet_ues_plus_Et , row.pfMet_ues_plus_Phi ), weight[n])

                                        
                histos[folder+'/ePFMET_Mt_jes_minus'].Fill(row.eMtToPfMet_jes_minus, weight[n])
                histos[folder+'/ePFMET_Mt_jes_plus'].Fill(row.eMtToPfMet_jes_plus,   weight[n])
                
                histos[folder+'/ePFMET_Mt_mes_minus'].Fill(row.eMtToPfMet_mes_minus,weight[n]) 
                histos[folder+'/ePFMET_Mt_mes_plus'].Fill(row.eMtToPfMet_mes_plus,  weight[n]) 
                
                histos[folder+'/ePFMET_Mt_ees_minus'].Fill(row.eMtToPfMet_ees_minus,weight[n]) 
                histos[folder+'/ePFMET_Mt_ees_plus'].Fill(row.eMtToPfMet_ees_plus,  weight[n]) 
                
                histos[folder+'/ePFMET_Mt_tes_minus'].Fill(row.eMtToPfMet_tes_minus,weight[n]) 
                histos[folder+'/ePFMET_Mt_tes_plus'].Fill(row.eMtToPfMet_tes_plus,  weight[n]) 
                
                histos[folder+'/ePFMET_Mt_ues_minus'].Fill(row.eMtToPfMet_ues_minus,weight[n]) 
                histos[folder+'/ePFMET_Mt_ues_plus'].Fill(row.eMtToPfMet_ues_plus,  weight[n]) 
                
                histos[folder+'/tPFMET_Mt_jes_plus'].Fill(row.tMtToPfMet_jes_plus,weight[n])
                histos[folder+'/tPFMET_Mt_mes_plus'].Fill(row.tMtToPfMet_mes_plus,weight[n])
                histos[folder+'/tPFMET_Mt_ees_plus'].Fill(row.tMtToPfMet_ees_plus,weight[n])
                histos[folder+'/tPFMET_Mt_tes_plus'].Fill(row.tMtToPfMet_tes_plus,weight[n])
                histos[folder+'/tPFMET_Mt_ues_plus'].Fill(row.tMtToPfMet_ues_plus,weight[n])
                histos[folder+'/tPFMET_Mt_jes_minus'].Fill(row.tMtToPfMet_jes_minus,weight[n])
                histos[folder+'/tPFMET_Mt_mes_minus'].Fill(row.tMtToPfMet_mes_minus,weight[n])
                histos[folder+'/tPFMET_Mt_ees_minus'].Fill(row.tMtToPfMet_ees_minus,weight[n])
                histos[folder+'/tPFMET_Mt_tes_minus'].Fill(row.tMtToPfMet_tes_minus,weight[n])
                histos[folder+'/tPFMET_Mt_ues_minus'].Fill(row.tMtToPfMet_ues_minus,weight[n])
                


 
            


    def process(self):
        
        

        central_weights = fakerate_central_histogram(25,0, 2.5)
        
        p1s_weights = fakerate_p1s_histogram(25,0, 2.5)#fakerate_p1s_histogram(25,0, 2.5)
        
        m1s_weights = fakerate_m1s_histogram(25,0, 2.5)#fakerate_m1s_histogram(25,0, 2.5)
        
                
        
        frw = []
        myevent =()
        for row in self.tree:
        #for n, row in enumerate(self.tree):
            
            sign = 'ss' if row.e_t_SS else 'os'
            processtype = '' ## use a line as for sign when the vbf when selections are defined            
            ptthreshold = [30]
            processtype ='gg'##changed from 20
            jn = row.jetVeto30
            if jn > 3 : jn = 3
            jn_jes_plus = row.jetVeto30jes_plus
            jn_jes_minus = row.jetVeto30jes_minus

            if jn_jes_plus >3 :  jn_jes_plus=3
            if jn_jes_minus >3 : jn_jes_minus=3

            #if row.run > 2 : #apply the trigger to data only (MC triggers enter in the scale factors)
            
            if self.is_embedded :

                if not bool(row.doubleMuPass) : continue
            else: 
                if not bool(row.singleE27WP80Pass) : continue
                if  not  bool(row.eMatchesSingleE27WP80): continue
                        
            if not selections.eSelection(row, 'e'): continue
               
            if not selections.lepton_id_iso(row, 'e', 'eid13Tight_etauiso01'): continue
                        
            if abs(row.eEta) > 1.4442 and abs(row.eEta) < 1.566 : continue
            if not selections.tauSelection(row, 't'): continue
                        
            if row.ePt < 30 : continue
            #if row.eMtToPFMET <40 : continue
            
            if not row.tAntiElectronMVA5Tight : continue
            if not row.tAntiMuon2Loose : continue
            if not row.tLooseIso3Hits : continue
            
            #isTauTight = False
            frw=self.fakerate_weights(row.tEta, central_weights, p1s_weights, m1s_weights )

            #if row.tauVetoPt20EleTight3MuLoose : continue 
            #if row.muVetoPt5IsoIdVtx : continue
            #if row.eVetoCicLooseIso : continue # change it with Loose

            if row.tauVetoPt20EleTight3MuLoose  and row.tauVetoPt20EleTight3MuLoose_tes_minus and row.tauVetoPt20EleTight3MuLoose_tes_plus: continue 
            if row.muVetoPt5IsoIdVtx and row.muVetoPt5IsoIdVtx_mes_minus and row.muVetoPt5IsoIdVtx_mes_plus : continue
            if row.eVetoCicLooseIso and row.eVetoCicLooseIso_ees_minus and row.eVetoCicLooseIso_ees_plus : continue
            
            standardSelection=True
            tesminus =True
            tesplus  =True
            mesminus =True
            mesplus  =True
            eesminus =True
            eesplus  =True        
            
            if  row.tauVetoPt20EleTight3MuLoose or row.muVetoPt5IsoIdVtx or  row.eVetoCicLooseIso : standardSelection =  False             

            if  row.tauVetoPt20EleTight3MuLoose_tes_minus or row.muVetoPt5IsoIdVtx or  row.eVetoCicLooseIso : tesminus =  False 
            if  row.tauVetoPt20EleTight3MuLoose_tes_plus or row.muVetoPt5IsoIdVtx or  row.eVetoCicLooseIso  : tesplus  =  False 
            if  row.tauVetoPt20EleTight3MuLoose or row.muVetoPt5IsoIdVtx_mes_minus or  row.eVetoCicLooseIso : mesminus =  False 
            if  row.tauVetoPt20EleTight3MuLoose or row.muVetoPt5IsoIdVtx_mes_plus or  row.eVetoCicLooseIso  : mesplus  =  False 
            if  row.tauVetoPt20EleTight3MuLoose or row.muVetoPt5IsoIdVtx or  row.eVetoCicLooseIso_ees_minus : eesminus =  False 
            if  row.tauVetoPt20EleTight3MuLoose or row.muVetoPt5IsoIdVtx or  row.eVetoCicLooseIso_ees_plus  : eesplus  =  False 
            
            dirpaths = [(standardSelection, sign+'/'+processtype), (mesplus, 'mVetoUp/'+sign+'/'+processtype),  (mesminus, 'mVetoDown/'+sign+'/'+processtype), \
                        (eesplus, 'eVetoUp/'+sign+'/'+processtype),  (eesminus, 'eVetoDown/'+sign+'/'+processtype), \
                        (tesplus, 'tVetoUp/'+sign+'/'+processtype),  (tesminus, 'tVetoDown/'+sign+'/'+processtype)]

            
            
            if (row.run, row.lumi, row.evt)==myevent: continue
            if myevent!=() and (row.run, row.lumi, row.evt)==(myevent[0], myevent[1], myevent[2]): print row.ePt, row.tPt
            
            myevent=(row.run, row.lumi, row.evt)

            isTauTight = bool(row.tTightIso3Hits)
            folder = dirpaths[0][1]+'/ept30'
            if dirpaths[0][0] == True and sign=='os':
                self.fill_histos(row, folder,isTauTight, frw)
            
            if row.bjetCSVVeto30!=0 : continue # it is better vetoing on b-jets  after the histo for the DY embedded

            for n,dirpath in enumerate(dirpaths):
                jetlist = [(int(jn), str(int(jn)))]
                if  dirpath[0]==False : continue 
                if n==0:
                    jetlist.extend([(int(jn_jes_plus), str(int(jn_jes_plus))+'_jes_plus'), (int(jn_jes_minus), str(int(jn_jes_plus))+'_jes_minus')])
                for jet in jetlist:
                    #for j in ptthreshold:
                    folder = dirpath[1]+'/ept30/'+jet[1]
                    #print folder
                                        
                        #print row.tLooseIso3Hits, row.tTightIso3Hits, isTauTight
                                                
                    self.fill_histos(row, folder,isTauTight, frw)
                
                    if jet[0] == 0 :
                        if row.tPt < 35: continue 
                        if row.ePt < 40 : continue
                        if deltaPhi(row.ePhi, row.tPhi) < 2.7 : continue
                        if row.tMtToPFMET > 50 : continue
                    if jet[0] == 1 :
                        if row.tPt < 40: continue 
                        if row.ePt < 35 : continue
                        if row.tMtToPFMET > 35 : continue
                    if jet[0] == 2 :
                        if row.tPt < 40: continue 
                        if row.ePt < 30 : continue # no cut as only electrons with pt>30 are in the ntuples
                        if row.tMtToPFMET > 35 : continue
                        if row.vbfMass < 550 : continue
                        if row.vbfDeta < 3.5 : continue
                    folder = dirpath[1]+'/ept30/'+jet[1]+'/selected'
                    self.fill_histos(row, folder, isTauTight,frw)
                
                 
                    
             
            
    def finish(self):
        self.write_histos()


