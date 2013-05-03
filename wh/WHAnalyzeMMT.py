'''

Analyze MMT events for the WH analysis

'''

from MuMuTauTree import MuMuTauTree
import ROOT
import os
from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits
import math

################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeMMT(WHAnalyzerBase):
    tree = 'mmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeMMT, self).__init__(tree, outfile, MuMuTauTree, **kwargs)
        self.hfunc['subMTMass'] = lambda row, weight: (row.m2_t_Mass, weight) if row.m1MtToMET > row.m2MtToMET else (row.m1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.hfunc['pt_ratio' ] = lambda row, weight: (row.m2Pt/row.m1Pt, weight)
##         self.hfunc['mass'     ] = lambda row, weight: (inv_mass(\
##                                                                 (row.m1Pt, row.m1Eta, row.m1Phi, 0.),\
##                                                                 (row.m2Pt, row.m2Eta, row.m2Phi, 0.),\
##                                                                 (row.tPt , row.tEta , row.tPhi , 0.),\
##                                                                 ), weight)
##         self.hfunc['lepRecoil'] = lambda row, weight: ( \
##                                                         quad( (row.m1Pt*math.cos(row.m1Phi) + row.m2Pt*math.cos(row.m2Phi) ), \
##                                                               (row.m1Pt*math.sin(row.m1Phi) + row.m2Pt*math.sin(row.m2Phi) ), ),\
##                                                         weight)
##         self.hfunc['lepRecoil_wMET'] = lambda row, weight: ( \
##                                                         quad( (row.m1Pt*math.cos(row.m1Phi) + row.m2Pt*math.cos(row.m2Phi) + row.metEt*math.cos(row.metPhi) ), \
##                                                               (row.m1Pt*math.sin(row.m1Phi) + row.m2Pt*math.sin(row.m2Phi) + row.metEt*math.sin(row.metPhi) ), ),\
##                                                         weight)
        self.hfunc["_recoilDaught" ] = lambda row, weight: (math.sqrt(row.recoilDaught) , weight)
        self.hfunc["_recoilWithMet"] = lambda row, weight: (math.sqrt(row.recoilWithMet), weight)

        #MC ONLY
        self.hfunc['higgsLMtToMet'] = lambda row, weight: ((row.m1MtToMET,row.m2MtToMET), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2MtToMET,row.m1MtToMET), weight)
        self.hfunc['higgsLIso']     = lambda row, weight: ((row.m1RelPFIsoDB,row.m2RelPFIsoDB), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2RelPFIsoDB,row.m1RelPFIsoDB), weight)
        self.hfunc['higgsLPt']      = lambda row, weight: ((row.m1Pt,row.m2Pt), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2Pt,row.m1Pt), weight)
        #self.hfunc['higgsLPt']      = lambda row, weight: (((row.m1Pt-row.m2Pt),(row.m1RelPFIsoDB-row.m2RelPFIsoDB)), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2Pt,row.m1Pt), weight)

        self.hfunc['higgsLMtToMet_1d'] = lambda row, weight: ((row.m1MtToMET-row.m2MtToMET), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2MtToMET-row.m1MtToMET), weight)
        self.hfunc['higgsLIso_1d']     = lambda row, weight: ((row.m1RelPFIsoDB-row.m2RelPFIsoDB), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2RelPFIsoDB-row.m1RelPFIsoDB), weight)
        self.hfunc['higgsLPt_1d']      = lambda row, weight: ((row.m1Pt-row.m2Pt), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2Pt-row.m1Pt), weight)

        self.hfunc['higgsTDR_1d']   = lambda row, weight: ((row.m1_t_DR-row.m2_t_DR), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2_t_DR-row.m1_t_DR), weight)
        self.hfunc['higgsTPt_1d']   = lambda row, weight: ((row.m1_t_Pt-row.m2_t_Pt), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2_t_Pt-row.m1_t_Pt), weight)
        self.hfunc['higgsDPhiMet']  = lambda row, weight: ((row.m1ToMETDPhi-row.m2ToMETDPhi), weight) if bool(row.m1ComesFromHiggs)  else ((row.m2ToMETDPhi-row.m1ToMETDPhi), weight)

        self.hfunc['H_LMtToMet'] = lambda row, weight: (row.m1MtToMET, weight) if bool(row.m1ComesFromHiggs)  else (row.m2MtToMET, weight)
        self.hfunc['H_LIso']     = lambda row, weight: (row.m1RelPFIsoDB, weight) if bool(row.m1ComesFromHiggs)  else (row.m2RelPFIsoDB, weight)
        self.hfunc['H_LPt']      = lambda row, weight: (row.m1Pt, weight) if bool(row.m1ComesFromHiggs)  else (row.m2Pt, weight)
        self.hfunc['W_LMtToMet'] = lambda row, weight: (row.m2MtToMET, weight) if bool(row.m1ComesFromHiggs)  else (row.m1MtToMET, weight)
        self.hfunc['W_LIso']     = lambda row, weight: (row.m2RelPFIsoDB, weight) if bool(row.m1ComesFromHiggs)  else (row.m1RelPFIsoDB, weight)
        self.hfunc['W_LPt']      = lambda row, weight: (row.m2Pt, weight) if bool(row.m1ComesFromHiggs)  else (row.m1Pt, weight)
        self.hfunc['true_mass']  = lambda row, weight: (row.m1_t_Mass, weight) if bool(row.m1ComesFromHiggs)  else (row.m2_t_Mass, weight)
        
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')

    def book_histos(self, folder):
        self.book(folder, "weight", "Event weight", 100, 0, 5)
        #self.book(folder, "weight_nopu", "Event weight without PU", 100, 0, 5) #Booked but not filled??
        self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
        self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)
        self.book(folder, "m1Pt", "Muon 1 Pt", 100, 0, 100)
        self.book(folder, "m1JetPt", "Muon 1 Jet Pt", 100, 0, 200)
        self.book(folder, "m2Pt", "Muon 2 Pt", 100, 0, 100)
        self.book(folder, "m2JetBtag", "Muon 2 Pt", 100, -10, 3.3)
        self.book(folder, "m2JetPt", "Muon 2 Jet Pt", 100, 0, 200)
        self.book(folder, "m1AbsEta", "Muon 1 AbsEta", 100, 0, 2.4)
        self.book(folder, "m2AbsEta", "Muon 2 AbsEta", 100, 0, 2.4)
        self.book(folder, "m1_m2_Mass", "Muon 1-2 Mass", 120, 0, 120)
        self.book(folder, "m2_t_Mass", "subleadingMass", 200, 0, 200)
        self.book(folder, "m1_t_Mass", "leadingMass", 200, 0, 200)
        # Rank muons by less MT to MET, for WZ control region
        self.book(folder, "subMTMass", "subMTMass", 200, 0, 200)
        self.book(folder, "m2RelPFIsoDB", "m2Iso", 100, 0, 0.3)
        self.book(folder, "tPt", "Tau Pt", 100, 0, 100)
        self.book(folder, "tAbsEta", "Tau AbsEta", 100, 0, 2.3)
        #self.book(folder, "tDecayMode", "Tau AbsEta", 15, -0.5, 14.5)
        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
        self.book(folder, "m1DZ",  "m1DZ", 100, 0., 1)
        self.book(folder, "m2DZ",  "m2DZ", 100, 0., 1)
        self.book(folder, "tDZ" ,  "tDZ" , 100, 0., 1)
        self.book(folder, "LT" ,  "LT" , 100, 0., 500)

        #let's look for osme other possible selections
        self.book(folder, "pt_ratio"      , "pt_ratio"      , 100, 0, 1)
        self.book(folder, "tToMETDPhi"    , "tToMETDPhi"    , 100, 0, 4)
        self.book(folder, "_recoilDaught"  , "recoilDaught"  , 600, 0, 8000)
        self.book(folder, "_recoilWithMet" , "recoilWithMet" , 600, 0, 8000)
        #self.book(folder, "lepRecoil"     , "lepRecoil"     , 600, 0, 8000)
        #self.book(folder, "lepRecoil_wMET", "lepRecoil_wMET", 600, 0, 8000)
        #self.book(folder, "mass"          , "mass"          , 800, 0, 800 )
        self.book(folder, "type1_pfMetEt"         , "metEt"         , 300, 0, 2000)
       
        ## self.book(folder, "tToMETDPhi#metEt", 100, 0, 4, 100, 0, 800, type=ROOT.TH2F)
        ## self.book(folder, "recoilWithMet#metEt", "recoilWithMet#metEt", 100, 0, 800, 100, 0, 800, type=ROOT.TH2F)
        
        #Book additial histograms for signal MC
##         if 'VH' in os.environ['megatarget'] and folder == 'ss/p1p2p3' and 'VHTests' in os.environ and os.environ['VHTests'] == 'YES':
##             self.book(folder, "true_mass", "True Mass", 200, 0, 200)
##             self.book(folder, "higgsLPt", "p_{T} lepton from higgs vs p_{T} lepton from W", 100, 0, 100, 100, 0, 100, type=ROOT.TH2F)
##             self.book(folder, "higgsLIso", "Isolation lepton from higgs vs Isolation lepton from W", 100, 0, 0.3, 100, 0, 0.3, type=ROOT.TH2F)
##             self.book(folder, "higgsLMtToMet", "M_{T} lepton from higgs vs M_{T} lepton from W", 100, 0, 200, 100, 0, 200, type=ROOT.TH2F)
##             self.book(folder, 'higgsLMtToMet_1d', "difference between lepton coming from higgs and the one from W", 100, -200, 200)
##             self.book(folder, 'higgsLIso_1d'    , "difference between lepton coming from higgs and the one from W", 100, -0.3, 0.3)
##             self.book(folder, 'higgsLPt_1d'     , "difference between lepton coming from higgs and the one from W", 100, -100, 100)
##             self.book(folder, 'higgsTDR_1d', "", 100, -10, 10)
##             self.book(folder, 'higgsTPt_1d', "", 100, -100, 100)
##             self.book(folder, 'higgsDPhiMet', "", 100, -7,7)
##             self.book(folder, 'H_LMtToMet', "", 100, 0, 200)
##             self.book(folder, 'H_LIso'    , "", 100, 0, 0.3) 
##             self.book(folder, 'H_LPt'     , "", 100, 0, 100) 
##             self.book(folder, 'W_LMtToMet', "", 100, 0, 200)
##             self.book(folder, 'W_LIso'    , "", 100, 0, 0.3) 
##             self.book(folder, 'W_LPt'     , "", 100, 0, 100)

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.doubleMuPass:                  return False
        if row.m1Pt < row.m2Pt:                   return False
        if row.m1Pt < 20:                         return False
        if not selections.muSelection(row, 'm1'): return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        if not selections.muSelection(row, 'm2'): return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        if not selections.tauSelection(row, 't'): return False #applies basic selection (eta, pt > 20, DZ)
        if row.m1_m2_SS and row.m1_t_SS         : return False #remove three SS leptons

        if row.m1_m2_Mass < 20:                    return False
        if row.LT < selections.lt_lower_threshold: return False

        if not selections.vetos(row):              return False #applies mu bjet e additional tau vetoes

        if not row.tAntiElectronMVA3Loose: return False

        if not self.trigger_match_m1(row): return False
        if not self.trigger_match_m2(row): return False

        return True

    @staticmethod
    def trigger_match_m1(row):
        return True
        if row.m1DiMuonL3p5PreFiltered8  > 0 or \
           row.m1DiMuonL3PreFiltered7  > 0 or \
           row.m1SingleMu13L3Filtered13  > 0 or \
           row.m1SingleMu13L3Filtered17  > 0 or \
           row.m1DiMuonMu17Mu8DzFiltered0p2  > 0 or \
           row.m1L3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17:
            return True

    @staticmethod
    def trigger_match_m2(row):
        return True
        if row.m2DiMuonL3p5PreFiltered8  > 0 or \
           row.m2DiMuonL3PreFiltered7  > 0 or \
           row.m2SingleMu13L3Filtered13  > 0 or \
           row.m2SingleMu13L3Filtered17  > 0 or \
           row.m2DiMuonMu17Mu8DzFiltered0p2  > 0 or \
           row.m2L3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17:
            return True

    @staticmethod
    def sign_cut(row):
        ''' Returns true if muons are SS '''
        return bool(row.m1_m2_SS)

    @staticmethod
    def obj1_id(row):

    @staticmethod
    def obj2_id(row):
        return selections.leading_lepton_id_iso(row, 'm1')

    @staticmethod
    def obj3_id(row):
        return selections.subleading_lepton_id_iso(row, 'e2')

    @staticmethod
    def anti_wz(row):
        return row.tAntiMuonTight # and not row.tMuOverlap

    def enhance_wz(self, row):
        # Require the "tau" to be a muon, and require the third muon
        # to have M_Z +- 20
        if self.anti_wz(row):
            return False
        # Cut on m2 PT > 20
        #if row.m2Pt < 20:
            #return False
        # Make sure any Z is from m1
        m2_good_Z = bool(71 < row.m2_t_Mass < 111)
        return not m2_good_Z

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m1','m2') * \
            mcCorrectors.double_muon_trigger(row,'m1','m2')

    def obj1_weight(self, row):
        return frfits.highpt_mu_fr(max(row.m1JetPt, row.m1Pt))

    def obj2_weight(self, row):
        return frfits.lowpt_mu_fr(max(row.m2JetPt, row.m2Pt))

    def obj3_weight(self, row):
        return frfits.tau_fr(row.tPt)

    def obj1_qcd_weight(self, row):
        return frfits.highpt_mu_qcd_fr(max(row.m1JetPt, row.m1Pt))

    def obj2_qcd_weight(self, row):
        return frfits.lowpt_mu_qcd_fr(max(row.m2JetPt, row.m2Pt))

    def obj3_qcd_weight(self, row):
        return frfits.tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return row.m1_t_SS

    ## def obj1_charge_flip(self, row):
    ##     return 0

    ## def obj2_charge_flip(self, row):
    ##     return 0
