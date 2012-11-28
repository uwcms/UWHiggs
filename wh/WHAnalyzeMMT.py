'''

Analyze MMT events for the WH analysis

'''

from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import fnmatch
import glob
from MuMuTauTree import MuMuTauTree
import os
import FinalStateAnalysis.MetaData.data_views as data_views
import logging
data_views.log.setLevel(logging.INFO)
#import sys
#logging.basicConfig(stream=sys.stderr, level=logging.INFO)
import WHAnalyzerBase
import ROOT
from TwoDimFakeRate import TwoDimFakeRate
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits

################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeMMT(WHAnalyzerBase.WHAnalyzerBase):
    tree = 'mmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeMMT, self).__init__(tree, outfile, MuMuTauTree, **kwargs)
        self.hfunc['subMTMass'] = lambda row, weight: (row.m2_t_Mass, weight) if row.m1MtToMET > row.m2MtToMET else (row.m1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
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
        self.book(folder, "tDecayMode", "Tau AbsEta", 15, -0.5, 14.5)
        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
        self.book(folder, "m1DZ",  "m1DZ", 100, 0., 1)
        self.book(folder, "m2DZ",  "m2DZ", 100, 0., 1)
        self.book(folder, "tDZ" ,  "tDZ" , 100, 0., 1)
        self.book(folder, "LT" ,  "LT" , 100, 0., 500)

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

        if row.m1_m2_Mass < 20:    return False
        if row.LT < 80:            return False

        if not selections.vetos(row): return False #applies mu bjet e additional tau vetoes

        if not row.tAntiElectronLoose: return False
        if row.tCiCTightElecOverlap:   return False

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

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return bool(row.m1_m2_SS)

    def obj1_id(self, row):
        return bool(row.m1PFIDTight) and (
            row.m1RelPFIsoDB < 0.1 or
            (row.m1RelPFIsoDB < 0.15 and row.m1AbsEta < 1.479))

    def obj2_id(self, row):
        return bool(row.m2PFIDTight) and (
            row.m2RelPFIsoDB < 0.1 or
            (row.m2RelPFIsoDB < 0.15 and row.m2AbsEta < 1.479))

    def obj3_id(self, row):
        return bool(row.tLooseMVAIso)

    def anti_wz(self, row):
        return row.tAntiMuonTight and not row.tMuOverlap

    def enhance_wz(self, row):
        # Require the "tau" to be a muon, and require the third muon
        # to have M_Z +- 20
        if row.tAntiMuonTight or not row.tMuOverlap:
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
        return not row.m1_t_SS

    def obj1_charge_flip(self, row):
        return 0
