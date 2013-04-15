'''

Analyze EET events for the WH analysis

'''

from EETauTree import EETauTree
import os
from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
import ROOT
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits
import ROOT
import math

#mtr = frfits.mt_likelihood_ratio
################################################################################
#### Analysis logic ############################################################
################################################################################

def logic_cut_1(row, weight):
    if row.e1_e2_Mass > 81 and row.e1_e2_Mass < 101:
        if row.type1_pfMetEt < 70:
            return (0.,weight)
        else:
            (1.,weight)
    else:
        return (1.,weight)


class WHAnalyzeEET(WHAnalyzerBase):
    tree = 'eet/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeEET, self).__init__(tree, outfile, EETauTree, **kwargs)
        self.hfunc['subMTMass'] = lambda row, weight: (row.e2_t_Mass, weight) if row.e1MtToMET > row.e2MtToMET else (row.e1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.hfunc['pt_ratio' ] = lambda row, weight: (row.e2Pt/row.e1Pt, weight)
        self.hfunc['mass'     ] = lambda row, weight: (inv_mass(\
                                                                (row.e1Pt, row.e1Eta, row.e1Phi, 0.),\
                                                                (row.e2Pt, row.e2Eta, row.e2Phi, 0.),\
                                                                (row.tPt , row.tEta , row.tPhi , 0.),\
                                                                ), weight)
##         self.hfunc['logic_cut_1'] = lambda row, weight: \
##             (0.,weight) if row.e1_e2_Mass > 81 and row.e1_e2_Mass < 101 and row.metEt < 70 else \
##             (1.,weight)
        
##         self.hfunc['lepRecoil_wMET'] = lambda row, weight: ( \
##                                                         quad( (row.e1Pt*math.cos(row.e1Phi) + row.e2Pt*math.cos(row.e2Phi) + row.metEt*math.cos(row.metPhi) ), \
##                                                               (row.e1Pt*math.sin(row.e1Phi) + row.e2Pt*math.sin(row.e2Phi) + row.metEt*math.sin(row.metPhi) ), ),\
##                                                         weight)
        self.hfunc["e*1_e2_Mass"] = lambda row, weight: ( frfits.mass_scaler( row.e1_e2_Mass), weight)
        self.hfunc["e*1_t_Mass" ] = lambda row, weight: ( frfits.mass_scaler( row.e1_t_Mass ), weight)
        self.hfunc["e1_e*2_Mass"] = lambda row, weight: ( frfits.mass_scaler( row.e1_e2_Mass), weight)
        self.hfunc["e*2_t_Mass" ] = lambda row, weight: ( frfits.mass_scaler( row.e2_t_Mass ), weight)
        self.hfunc["_recoilDaught" ] = lambda row, weight: (math.sqrt(row.recoilDaught) , weight)
        self.hfunc["_recoilWithMet"] = lambda row, weight: (math.sqrt(row.recoilWithMet), weight)
        self.hfunc["e1eta_on_z_peak"] = lambda row, weight: ( row.e1AbsEta, weight) if row.e1_e2_Mass > 80 and row.e1_e2_Mass < 100 else (-10,0.)
        self.hfunc["e1pt_on_z_peak" ] = lambda row, weight: ( row.e1Pt    , weight) if row.e1_e2_Mass > 80 and row.e1_e2_Mass < 100 else (-10,0.)
        self.hfunc["e2eta_on_z_peak"] = lambda row, weight: ( row.e2AbsEta, weight) if row.e1_e2_Mass > 80 and row.e1_e2_Mass < 100 else (-10,0.)
        self.hfunc["e2pt_on_z_peak" ] = lambda row, weight: ( row.e2Pt    , weight) if row.e1_e2_Mass > 80 and row.e1_e2_Mass < 100 else (-10,0.)
        

        
        #MC ONLY
        ## self.hfunc['higgsLMtToMet'] = lambda row, weight: ((row.e1MtToMET,row.e2MtToMET), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2MtToMET,row.e1MtToMET), weight)
        ## self.hfunc['higgsLIso']     = lambda row, weight: ((row.e1RelPFIsoDB,row.e2RelPFIsoDB), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2RelPFIsoDB,row.e1RelPFIsoDB), weight)
        ## self.hfunc['higgsLPt']      = lambda row, weight: ((row.e1Pt,row.e2Pt), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2Pt,row.e1Pt), weight)

        ## self.hfunc['higgsLMtToMet_1d'] = lambda row, weight: ((row.e1MtToMET-row.e2MtToMET), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2MtToMET-row.e1MtToMET), weight)
        ## self.hfunc['higgsLIso_1d']     = lambda row, weight: ((row.e1RelPFIsoDB-row.e2RelPFIsoDB), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2RelPFIsoDB-row.e1RelPFIsoDB), weight)
        ## self.hfunc['higgsLPt_1d']      = lambda row, weight: ((row.e1Pt-row.e2Pt), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2Pt-row.e1Pt), weight)
        ## self.hfunc['higgsTDR_1d']   = lambda row, weight: ((row.e1_t_DR-row.e2_t_DR), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2_t_DR-row.e1_t_DR), weight)
        ## self.hfunc['higgsTPt_1d']   = lambda row, weight: ((row.e1_t_Pt-row.e2_t_Pt), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2_t_Pt-row.e1_t_Pt), weight)
        ## self.hfunc['higgsDPhiMet']  = lambda row, weight: ((row.e1ToMETDPhi-row.e2ToMETDPhi), weight) if bool(row.e1ComesFromHiggs)  else ((row.e2ToMETDPhi-row.e1ToMETDPhi), weight)

        ## self.hfunc['H_LMtToMet'] = lambda row, weight: (row.e1MtToMET, weight) if bool(row.e1ComesFromHiggs)  else (row.e2MtToMET, weight)
        ## self.hfunc['H_LIso']     = lambda row, weight: (row.e1RelPFIsoDB, weight) if bool(row.e1ComesFromHiggs)  else (row.e2RelPFIsoDB, weight)
        ## self.hfunc['H_LPt']      = lambda row, weight: (row.e1Pt, weight) if bool(row.e1ComesFromHiggs)  else (row.e2Pt, weight)
        ## self.hfunc['W_LMtToMet'] = lambda row, weight: (row.e2MtToMET, weight) if bool(row.e1ComesFromHiggs)  else (row.e1MtToMET, weight)
        ## self.hfunc['W_LIso']     = lambda row, weight: (row.e2RelPFIsoDB, weight) if bool(row.e1ComesFromHiggs)  else (row.e1RelPFIsoDB, weight)
        ## self.hfunc['W_LPt']      = lambda row, weight: (row.e2Pt, weight) if bool(row.e1ComesFromHiggs)  else (row.e1Pt, weight)
        ## self.hfunc['true_mass']  = lambda row, weight: (row.e1_t_Mass, weight) if bool(row.e1ComesFromHiggs)  else (row.e2_t_Mass, weight)


        self.pucorrector = mcCorrectors.make_puCorrector('doublee')

    def book_histos(self, folder):
        self.book(folder, "weight", "Event weight", 100, 0, 5)
        #self.book(folder, "weight_nopu", "Event weight without PU", 100, 0, 5)
        self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)

        self.book(folder, "e1Pt", "E 1 Pt", 100, 0, 100)
        self.book(folder, "e2Pt", "E 2 Pt", 100, 0, 100)

        self.book(folder, "e1AbsEta", "Muon 1 AbsEta", 100, 0, 2.5)
        self.book(folder, "e2AbsEta", "Muon 2 AbsEta", 100, 0, 2.5)

        self.book(folder, "e1_e2_Mass", "E 1-2 Mass", 120, 0, 120)
        self.book(folder, "e1_t_Mass", "leadingMass", 200, 0, 200)
        self.book(folder, "e2_t_Mass", "subleadingMass", 200, 0, 200)

        self.book(folder, "e1_e2_CosThetaStar", "subleadingMass", 200, -1., 1.)
        self.book(folder, "e1_t_CosThetaStar" , "subleadingMass", 200, -1., 1.)
        self.book(folder, "e2_t_CosThetaStar" , "subleadingMass", 200, -1., 1.)

        self.book(folder, "e2RelPFIsoDB", "e2RelPFIsoDB", 100, 0, 0.3)
        self.book(folder, "tPt", "tPt", 100, 0,100)
        self.book(folder, "tAbsEta", "tAbsEta", 100, 0, 2.3)
        #self.book(folder, "metSignificance", "MET significance", 100, 0, 15)
        self.book(folder, "LT", "L_T", 100, 0, 300)

        #Charge mis-id special histograms
        if 'c1' in folder:
            self.book(folder, "e*1_e2_Mass", "E 1-2 Mass with misid sclaing correction", 120, 0, 120)
            self.book(folder, "e*1_t_Mass", "leadingMass with misid sclaing correction", 200, 0, 200)
        elif 'c2' in folder:
            self.book(folder, "e1_e*2_Mass", "E 1-2 Mass with misid sclaing correction", 120, 0, 120)
            self.book(folder, "e*2_t_Mass", "subleadingMass with misid sclaing correction", 200, 0, 200)

        if 'f3' in folder:
            self.book(folder, "e1eta_on_z_peak", "Muon 1 AbsEta", 100, 0, 2.4)
            self.book(folder, "e1pt_on_z_peak" , "E 1 Pt", 100, 0, 100)
            self.book(folder, "e2eta_on_z_peak", "Muon 1 AbsEta", 100, 0, 2.4)
            self.book(folder, "e2pt_on_z_peak" , "E 1 Pt", 100, 0, 100)


        #let's look for osme other possible selections
        self.book(folder, "mass"          , "mass"          , 800, 0, 800 )
        self.book(folder, "pt_ratio"      , "pt_ratio"      , 100, 0, 1)
        self.book(folder, "tToMETDPhi"    , "tToMETDPhi"    , 100, 0, 4)
        self.book(folder, "_recoilDaught"  , "recoilDaught"  , 600, 0, 8000)
        self.book(folder, "_recoilWithMet" , "recoilWithMet" , 600, 0, 8000)
        self.book(folder, "e1_e2_Pt"       , "lepRecoil"     , 600, 0, 8000)
        self.book(folder, "e1_e2_DR"       , "e1_e2_DR"      , 500, 0, 10)
        self.book(folder, "m1_m2_CosThetaStar", "m1_m2_CosThetaStar", 200, -1, 1)
        self.book(folder, "m1_t_CosThetaStar" , "m1_t_CosThetaStar" , 200, -1, 1)
        self.book(folder, "m2_t_CosThetaStar" , "m2_t_CosThetaStar" , 200, -1, 1)
        #self.book(folder, "lepRecoil_wMET", "lepRecoil_wMET", 600, 0, 8000)
        self.book(folder, "type1_pfMetEt"         , "metEt"         , 300, 0, 2000)
        self.book(folder, "type1_pfMetEt#e1_e2_Mass", "metEt#e1_e2_Mass", 100, 0, 300, 120, 0, 120, type=ROOT.TH2F)
        #self.book(folder, "logic_cut_1" ,"logic_cut_1", 2, 0.,2.)

        #Book additial histograms for signal MC
        ## if 'VH' in os.environ['megatarget'] and folder == 'ss/p1p2p3' and 'VHTests' in os.environ and os.environ['VHTests'] == 'YES':
        ##     self.book(folder, "true_mass", "True Mass", 200, 0, 200)
        ##     self.book(folder, "higgsLPt"        , "p_{T} lepton from higgs vs p_{T} lepton from W", 100, 0, 100, 100, 0, 100, type=ROOT.TH2F)
        ##     self.book(folder, "higgsLIso"       , "Isolation lepton from higgs vs Isolation lepton from W", 100, 0, 0.3, 100, 0, 0.3, type=ROOT.TH2F)
        ##     self.book(folder, "higgsLMtToMet"   , "M_{T} lepton from higgs vs M_{T} lepton from W", 100, 0, 200, 100, 0, 200, type=ROOT.TH2F)
        ##     self.book(folder, 'higgsLMtToMet_1d', "difference between lepton coming from higgs and the one from W", 100, -200, 200)
        ##     self.book(folder, 'higgsLIso_1d'    , "difference between lepton coming from higgs and the one from W", 100, -0.3, 0.3)
        ##     self.book(folder, 'higgsLPt_1d'     , "difference between lepton coming from higgs and the one from W", 100, -100, 100)
        ##     self.book(folder, 'higgsTDR_1d', "", 100, -10, 10)
        ##     self.book(folder, 'higgsTPt_1d', "", 100, -100, 100)
        ##     self.book(folder, 'higgsDPhiMet', "", 100, -7,7)
        ##     self.book(folder, 'H_LMtToMet', "", 100, 0, 200)
        ##     self.book(folder, 'H_LIso'    , "", 100, 0, 0.3) 
        ##     self.book(folder, 'H_LPt'     , "", 100, 0, 100) 
        ##     self.book(folder, 'W_LMtToMet', "", 100, 0, 200)
        ##     self.book(folder, 'W_LIso'    , "", 100, 0, 0.3) 
        ##     self.book(folder, 'W_LPt'     , "", 100, 0, 100) 
            
    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def preselection(row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        #object basic selections
        if not row.doubleEPass:                   return False
        if row.e1Pt < 20:                         return False
        if not selections.eSelection(row, 'e1'):  return False
        if not selections.eSelection(row, 'e2'):  return False
        if not selections.tauSelection(row, 't'): return False
        if row.e1_e2_SS and row.e1_t_SS         : return False #remove three SS leptons

        if row.e1_e2_Mass < 20:       return False
            #if row.metSignificance < 2.5: return False
        ## if row.e1_e2_Mass > 81 \
        ##     and row.e1_e2_Mass < 101: return False
        if row.LT < 80:              return False
        if not selections.vetos(row): return False #applies mu bjet e additional tau vetoes

        if not row.tAntiMuonLoose:   return False
        #'t_ElectronOverlapWP95 < 0.5',
        return True

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def sign_cut(row):
        ''' Returns true if muons are SS '''
        return bool(row.e1_e2_SS)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj1_id(row):
        return selections.h2tau_eid(row, 'e1') #bool(row.e1MVAIDH2TauWP) and bool( row.e1RelPFIsoDB < 0.1 or (row.e1RelPFIsoDB < 0.15 and row.e1AbsEta < 1.479))

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj2_id(row):
        return selections.h2tau_eid(row, 'e2') #bool(row.e2MVAIDH2TauWP) and bool( row.e2RelPFIsoDB < 0.1 or (row.e2RelPFIsoDB < 0.15 and row.e2AbsEta < 1.479))

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj3_id(row):
        return bool(row.tLooseIso3Hits)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def anti_wz(row):
        if row.e_t_Zcompat < 20:
            if not row.tAntiElectronMVA3Tight:
                return False
        elif not row.tAntiElectronMVA3Medium:
            return False
        return True

    def enhance_wz(self, row):
        # Require the "tau" to be a electron, and require the third electron
        # to have M_Z +- 20 
        if self.anti_wz(row):
            return False
        # Make sure any Z is from e1
        e2_good_Z = bool(71 < row.e2_t_Mass < 111)
        return not e2_good_Z

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_electron_corrections(row,'e1','e2')

   
    def obj1_weight(self, row):
        return frfits.highpt_ee_fr(max(row.e1JetPt, row.e1Pt))

    def obj2_weight(self, row):
        return frfits.lowpt_ee_fr(max(row.e2JetPt, row.e2Pt))

    def obj3_weight(self, row):
        return frfits.tau_fr(row.tPt)

    def obj1_qcd_weight(self, row):
        return frfits.highpt_ee_qcd_fr(max(row.e1JetPt, row.e1Pt))

    def obj2_qcd_weight(self, row):
        return frfits.lowpt_ee_qcd_fr(max(row.e2JetPt, row.e2Pt))

    def obj3_qcd_weight(self, row):
        return frfits.tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return not row.e1_t_SS

    def obj1_charge_flip(self, row):
        return frfits.e_charge_flip(row.e1AbsEta,row.e1Pt) #highpt_e_charge_flip

    def obj2_charge_flip(self, row):
        return frfits.e_charge_flip(row.e2AbsEta,row.e2Pt) #lowpt_e_charge_flip

    def obj1_charge_flip_sysup(self, row):
        return frfits.e_charge_flip_up(row.e1AbsEta,row.e1Pt) #highpt_e_charge_flip

    def obj2_charge_flip_sysup(self, row):
        return frfits.e_charge_flip_up(row.e2AbsEta,row.e2Pt) #lowpt_e_charge_flip







