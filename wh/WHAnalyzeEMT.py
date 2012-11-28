'''

Analyze EMT events for the WH analysis

'''

import WHAnalyzerBase
from EMuTauTree import EMuTauTree
import glob
import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits

################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeEMT(WHAnalyzerBase.WHAnalyzerBase):
    tree = 'emt/final/Ntuple'

    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeEMT, self).__init__(tree, outfile, EMuTauTree, **kwargs)
        #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.hfunc['subMass'] = lambda row, weight: (row.e_t_Mass, weight) if row.ePt < row.mPt else (row.m_t_Mass, weight) 
        self.hfunc['tLeadDR'] = lambda row, weight: (row.m_t_DR,   weight) if row.ePt < row.mPt else (row.e_t_DR,   weight) 
        self.hfunc['tSubDR']  = lambda row, weight: (row.e_t_DR,   weight) if row.ePt < row.mPt else (row.m_t_DR,   weight) 
        self.pucorrector = mcCorrectors.make_puCorrector('mueg')

    def book_histos(self, folder):
        self.book(folder, "mPt", "Muon Pt", 100, 0, 100)
        self.book(folder, "ePt", "Electron Pt", 100, 0, 100)
        self.book(folder, "tPt", "Tau Pt", 100, 0, 100)
        self.book(folder, "mAbsEta", "Muon AbsEta", 100, 0, 2.4)
        self.book(folder, "eAbsEta", "Electron AbsEta", 100, 0, 2.5)
        self.book(folder, "tAbsEta", "Tau AbsEta", 100, 0, 2.3)

        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
        self.book(folder, "emMass", "Electron-Muon Mass", 200, 0, 200)
        self.book(folder, "eChargeIdTight", "Elec charge ID tight",
                  2, -0.5, 1.5)
        self.book(folder, "eChargeIdMedium", "Elec charge ID medium",
                  2, -0.5, 1.5)
        self.book(folder, "etMass", "Electron-Tau Mass", 200, 0, 200)
        self.book(folder, "subMass", "Subleading Mass", 200, 0, 200)
        self.book(folder, "bCSVVeto", "BjetCSV", 10, -0.5, 9.5)
        self.book(folder, "metSig", "MET significance", 100, 0, 15)
        self.book(folder, "tLeadDR", "DR between leading lepton and tau",
                  100, 0, 5)
        self.book(folder, "tSubDR", "DR between subleading lepton and tau",
                  100, 0, 5)


    def preselection(self, row):
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
            
        if row.LT < 80:            return False
        if row.tMuOverlap:         return False
        if not row.tAntiMuonTight: return False
        #'t_ElectronOverlapWP95 < 0.5',
        return True

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return bool(row.e_m_SS)

    def obj1_id(self, row):
        #return bool(row.mPFIDTight) and bool(row.mRelPFIsoDB < 0.2)
        return bool(row.mPFIDTight) and (
            row.mRelPFIsoDB < 0.1 or
            (row.mRelPFIsoDB < 0.15 and row.mAbsEta < 1.479))

    def obj2_id(self, row):
        return bool(row.eMVAIDH2TauWP) and bool(
            row.eRelPFIsoDB < 0.1 or
            (row.eRelPFIsoDB < 0.15 and row.eAbsEta < 1.479))

    def obj3_id(self, row):
        return bool(row.tLooseMVAIso)

    def anti_wz(self, row):
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
        return not row.e_t_SS

    def obj1_charge_flip(self, row):
        if row.eAbsEta < 1.5:
            return 0.003
        return 0.02
