'''

Analyze EMT events for the WH analysis

'''

from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
from EMuTauTree import EMuTauTree
import glob
import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits
import ROOT
import math

mtr = frfits.mt_likelihood_ratio
################################################################################
#### Analysis logic ############################################################
################################################################################
is7TeV = bool('7TeV' in os.environ['jobid'])
use_iso_trigger = not is7TeV

class WHAnalyzeEMT(WHAnalyzerBase):
    tree = 'emt/final/Ntuple'

    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeEMT, self).__init__(tree, outfile, EMuTauTree, **kwargs)
        #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.hfunc['subMass']   = lambda row, weight: (row.e_t_Mass, weight)    if row.ePt < row.mPt else (row.m_t_Mass, weight) 
        self.hfunc['tLeadDR']   = lambda row, weight: (row.m_t_DR,   weight)    if row.ePt < row.mPt else (row.e_t_DR,   weight) 
        self.hfunc['tSubDR']    = lambda row, weight: (row.e_t_DR,   weight)    if row.ePt < row.mPt else (row.m_t_DR,   weight) 
        self.hfunc['pt_ratio' ] = lambda row, weight: (row.ePt/row.mPt, weight) if row.ePt < row.mPt else (row.mPt/row.ePt, weight)
        self.hfunc['mass'     ] = lambda row, weight: (inv_mass(\
                                                                (row.mPt, row.mEta, row.mPhi, 0.),\
                                                                (row.ePt, row.eEta, row.ePhi, 0.),\
                                                                (row.tPt, row.tEta, row.tPhi, 0.),\
                                                                ), weight)
        self.hfunc['lepRecoil'] = lambda row, weight: ( \
                                                        quad( (row.mPt*math.cos(row.mPhi) + row.ePt*math.cos(row.ePhi) ), \
                                                              (row.mPt*math.sin(row.mPhi) + row.ePt*math.sin(row.ePhi) ), ),\
                                                        weight)
        self.hfunc['lepRecoil_wMET'] = lambda row, weight: ( \
                                                        quad( (row.mPt*math.cos(row.mPhi) + row.ePt*math.cos(row.ePhi) + row.metEt*math.cos(row.metPhi) ), \
                                                              (row.mPt*math.sin(row.mPhi) + row.ePt*math.sin(row.ePhi) + row.metEt*math.sin(row.metPhi) ), ),\
                                                        weight)
        self.hfunc["e*_t_Mass"] = lambda row, weight: ( frfits.mass_scaler( row.e_t_Mass), weight)
        self.hfunc["e*_m_Mass"] = lambda row, weight: ( frfits.mass_scaler( row.e_m_Mass), weight)
        self.hfunc["subMass*" ] = lambda row, weight: ( frfits.mass_scaler( row.e_t_Mass), weight)    if row.ePt < row.mPt else (row.m_t_Mass, weight)
        self.hfunc["_recoilDaught" ] = lambda row, weight: (math.sqrt(row.recoilDaught) , weight)
        self.hfunc["_recoilWithMet"] = lambda row, weight: (math.sqrt(row.recoilWithMet), weight)


        
        #MC ONLY
        self.hfunc['higgsLMtToMet'] = lambda row, weight: ((row.eMtToMET,row.mMtToMET), weight) if bool(row.eComesFromHiggs)  else ((row.mMtToMET,row.eMtToMET), weight)
        self.hfunc['higgsLIso']     = lambda row, weight: ((row.eRelPFIsoDB,row.mRelPFIsoDB), weight) if bool(row.eComesFromHiggs)  else ((row.mRelPFIsoDB,row.eRelPFIsoDB), weight)
        self.hfunc['higgsLPt']      = lambda row, weight: ((row.ePt,row.mPt), weight) if bool(row.eComesFromHiggs)  else ((row.mPt,row.ePt), weight)
        self.hfunc['higgsLMtToMet_1d'] = lambda row, weight: ((row.eMtToMET-row.mMtToMET), weight) if bool(row.eComesFromHiggs)  else ((row.mMtToMET-row.eMtToMET), weight)
        self.hfunc['higgsLIso_1d']     = lambda row, weight: ((row.eRelPFIsoDB-row.mRelPFIsoDB), weight) if bool(row.eComesFromHiggs)  else ((row.mRelPFIsoDB-row.eRelPFIsoDB), weight)
        self.hfunc['higgsLPt_1d']      = lambda row, weight: ((row.ePt-row.mPt), weight) if bool(row.eComesFromHiggs)  else ((row.mPt-row.ePt), weight)
        self.hfunc['higgsMtRatio_1d']  = lambda row, weight: ((mtr(row.eMtToMET)-mtr(row.mMtToMET)), weight) if bool(row.eComesFromHiggs)  else ((mtr(row.mMtToMET)-mtr(row.eMtToMET)), weight)
        self.hfunc['H_LMtToMet'] = lambda row, weight: (row.eMtToMET, weight) if bool(row.eComesFromHiggs)  else (row.mMtToMET, weight)
        self.hfunc['H_LIso']     = lambda row, weight: (row.eRelPFIsoDB, weight) if bool(row.eComesFromHiggs)  else (row.mRelPFIsoDB, weight)
        self.hfunc['H_LPt']      = lambda row, weight: (row.ePt, weight) if bool(row.eComesFromHiggs)  else (row.mPt, weight)
        self.hfunc['W_LMtToMet'] = lambda row, weight: (row.mMtToMET, weight) if bool(row.eComesFromHiggs)  else (row.eMtToMET, weight)
        self.hfunc['W_LIso']     = lambda row, weight: (row.mRelPFIsoDB, weight) if bool(row.eComesFromHiggs)  else (row.eRelPFIsoDB, weight)
        self.hfunc['W_LPt']      = lambda row, weight: (row.mPt, weight) if bool(row.eComesFromHiggs)  else (row.ePt, weight)
        self.hfunc['higgsTDR_1d']   = lambda row, weight: ((row.e_t_DR-row.m_t_DR), weight) if bool(row.eComesFromHiggs)  else ((row.m_t_DR-row.e_t_DR), weight)
        self.hfunc['higgsTPt_1d']   = lambda row, weight: ((row.e_t_Pt-row.m_t_Pt), weight) if bool(row.eComesFromHiggs)  else ((row.m_t_Pt-row.e_t_Pt), weight)
        self.hfunc['higgsDPhiMet']  = lambda row, weight: ((row.eToMETDPhi-row.mToMETDPhi), weight) if bool(row.eComesFromHiggs)  else ((row.mToMETDPhi-row.eToMETDPhi), weight)
        self.hfunc['true_mass']     = lambda row, weight: (row.e_t_Mass, weight) if bool(row.eComesFromHiggs)  else (row.m_t_Mass, weight)


        self.pucorrector = mcCorrectors.make_puCorrector('mueg')

    def book_histos(self, folder):
        self.book(folder, "mPt", "Muon Pt", 100, 0, 100)
        self.book(folder, "ePt", "Electron Pt", 100, 0, 100)
        self.book(folder, "tPt", "Tau Pt", 100, 0, 100)
        self.book(folder, "mAbsEta", "Muon AbsEta", 100, 0, 2.4)
        self.book(folder, "eAbsEta", "Electron AbsEta", 100, 0, 2.5)
        self.book(folder, "tAbsEta", "Tau AbsEta", 100, 0, 2.3)

        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
        self.book(folder, "e_m_Mass", "Electron-Muon Mass", 200, 0, 200)
        self.book(folder, "eChargeIdTight", "Elec charge ID tight", 2, -0.5, 1.5)
        #self.book(folder, "eChargeIdMedium", "Elec charge ID medium", 2, -0.5, 1.5)
        self.book(folder, "e_t_Mass", "Electron-Tau Mass", 200, 0, 200)
        self.book(folder, "subMass", "Subleading Mass", 200, 0, 200)
        #self.book(folder, "bCSVVeto", "BjetCSV", 10, -0.5, 9.5)
        #self.book(folder, "metSig", "MET significance", 100, 0, 15)
        self.book(folder, "tLeadDR", "DR between leading lepton and tau", 100, 0, 5)
        self.book(folder, "tSubDR", "DR between subleading lepton and tau", 100, 0, 5)
        self.book(folder, "LT", "L_T", 100, 0, 300)
        #Charge mis-id special histograms
        if 'c1' in folder:
            self.book(folder, "e*_t_Mass", "Electron-Tau Mass", 200, 0, 200)
            self.book(folder, "subMass*", "Subleading Mass", 200, 0, 200)
            self.book(folder, "e*_m_Mass", "Electron-Muon Mass", 200, 0, 200)
        
        #let's look for osme other possible selections
        self.book(folder, "mass"          , "mass"          , 800, 0, 800 )
        self.book(folder, "pt_ratio"      , "pt_ratio"      , 100, 0, 1)
        self.book(folder, "tToMETDPhi"    , "tToMETDPhi"    , 100, 0, 4)
        self.book(folder, "_recoilDaught"  , "recoilDaught"  , 600, 0, 8000)
        self.book(folder, "_recoilWithMet" , "recoilWithMet" , 600, 0, 8000)
        self.book(folder, "lepRecoil"     , "lepRecoil"     , 600, 0, 8000)
        self.book(folder, "lepRecoil_wMET", "lepRecoil_wMET", 600, 0, 8000)
        self.book(folder, "metEt"         , "metEt"         , 300, 0, 2000)

        #Book additial histograms for signal MC
        if 'VH' in os.environ['megatarget'] and folder == 'ss/p1p2p3' and 'VHTests' in os.environ and os.environ['VHTests'] == 'YES':
            self.book(folder, "true_mass", "True Mass", 200, 0, 200)
            self.book(folder, "higgsLPt", "p_{T} lepton from higgs vs p_{T} lepton from W", 100, 0, 100, 100, 0, 100, type=ROOT.TH2F)
            self.book(folder, "higgsLIso", "Isolation lepton from higgs vs Isolation lepton from W", 100, 0, 0.3, 100, 0, 0.3, type=ROOT.TH2F)
            self.book(folder, "higgsLMtToMet", "M_{T} lepton from higgs vs M_{T} lepton from W", 100, 0, 200, 100, 0, 200, type=ROOT.TH2F)
            self.book(folder, 'higgsLMtToMet_1d', "difference between lepton coming from higgs and the one from W", 100, -200, 200)
            self.book(folder, 'higgsLIso_1d'    , "difference between lepton coming from higgs and the one from W", 100, -0.3, 0.3)
            self.book(folder, 'higgsLPt_1d'     , "difference between lepton coming from higgs and the one from W", 100, -100, 100)
            self.book(folder, 'higgsMtRatio_1d' , "", 100, -10, 10)
            self.book(folder, 'higgsTDR_1d', "", 100, -10, 10)
            self.book(folder, 'higgsTPt_1d', "", 100, -100, 100)
            self.book(folder, 'higgsDPhiMet', "", 100, -7,7)
            self.book(folder, 'H_LMtToMet', "", 100, 0, 200)
            self.book(folder, 'H_LIso'    , "", 100, 0, 0.3) 
            self.book(folder, 'H_LPt'     , "", 100, 0, 100) 
            self.book(folder, 'W_LMtToMet', "", 100, 0, 200)
            self.book(folder, 'W_LIso'    , "", 100, 0, 0.3) 
            self.book(folder, 'W_LPt'     , "", 100, 0, 100) 


    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def preselection( row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not use_iso_trigger:
            if not row.mu17ele8Pass:
                return False
        elif not row.mu17ele8isoPass:
            return False
        if row.mPt < 20:            return False
        if not selections.muSelection(row, 'm'):  return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        if not selections.eSelection(row, 'e'):   return False #applies basic selection (eta, pt > 10, DZ, missingHits, jetBTag, HasConversion and chargedIdTight)
        if not selections.tauSelection(row, 't'): return False #applies basic selection (eta, pt > 20, DZ)
        if not selections.vetos(row):             return False #applies mu bjet e additional tau vetoes
        if row.eJetBtag > 3.3:                    return False
        if row.e_m_SS and row.e_t_SS            : return False #remove three SS leptons
            
        if row.LT < 80:            return False
        if row.tMuOverlap:         return False
        if not row.tAntiMuonTight: return False
        #'t_ElectronOverlapWP95 < 0.5',
        return True

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def sign_cut( row):
        ''' Returns true if muons are SS '''
        return bool(row.e_m_SS)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj1_id( row):
        #return bool(row.mPFIDTight) and bool(row.mRelPFIsoDB < 0.2)
        return bool(row.mPFIDTight) and (
            row.mRelPFIsoDB < 0.1 or
            (row.mRelPFIsoDB < 0.15 and row.mAbsEta < 1.479))

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj2_id( row):
        return bool(row.eMVAIDH2TauWP) and bool(
            row.eRelPFIsoDB < 0.1 or
            (row.eRelPFIsoDB < 0.15 and row.eAbsEta < 1.479))

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj3_id( row):
        return bool(row.tLooseMVAIso)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def anti_wz( row):
        if row.tCiCTightElecOverlap:
            return False
        if row.e_t_Zcompat < 20:
            if not row.tAntiElectronMVA:
                return False
        elif not row.tAntiElectronLoose:
            return False
        return True

    def enhance_wz(self, row):
        if row.e_t_Zcompat < 15 and not row.tAntiElectronMVA:
            return True
        return False

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m') * \
            mcCorrectors.get_electron_corrections(row,'e') * \
            mcCorrectors.correct_mueg_mu(row.mPt, row.mAbsEta) * \
            mcCorrectors.correct_mueg_e(row.ePt, row.eAbsEta)

    def obj1_weight(self, row):
        return frfits.highpt_mu_fr(max(row.mJetPt, row.mPt))

    def obj2_weight(self, row):
        return frfits.lowpt_e_fr(max(row.eJetPt, row.ePt))

    def obj3_weight(self, row):
        return frfits.tau_fr(row.tPt)

    def obj1_qcd_weight(self, row):
        return frfits.highpt_mu_qcd_fr(max(row.mJetPt, row.mPt))

    def obj2_qcd_weight(self, row):
        return frfits.lowpt_e_qcd_fr(max(row.eJetPt, row.ePt))

    def obj3_qcd_weight(self, row):
        return frfits.tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return row.m_t_SS

    def obj1_charge_flip(self, row):
        return 0

    def obj2_charge_flip(self, row):
        return frfits.e_charge_flip(row.eAbsEta,row.ePt)

