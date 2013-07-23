'''

Analyze EEMT events for the ZH analysis

'''

#from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import glob
from EEMuTauTree import EEMuTauTree
import os
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import baseSelections as selections
import mcCorrectors
import ZHAnalyzerBase
import ROOT
import fake_rate_functions as fr_fcn

################################################################################
#### Analysis logic ############################################################
################################################################################

class ZHAnalyzeEEMT(ZHAnalyzerBase.ZHAnalyzerBase):
    tree = 'eemt/final/Ntuple'
    name = 5
    def __init__(self, tree, outfile, **kwargs):
        super(ZHAnalyzeEEMT, self).__init__(tree, outfile, EEMuTauTree, 'MT', **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        target = os.environ['megatarget']
        self.pucorrector = mcCorrectors.make_puCorrector('doublee')
        ## if 'HWW3l' in target:
        ##     print "HACK using S6 PU weights for HWW3l"
        ##     mcCorrectors.force_pu_distribution('S6')

    def Z_decay_products(self):
        return ('e1','e2')

    def H_decay_products(self):
        return ('m','t')

    def book_histos(self, folder):
        self.book_general_histos(folder)
        self.book_kin_histos(folder, 'e1')
        self.book_kin_histos(folder, 'e2')
        self.book_kin_histos(folder, 'm')
        self.book_kin_histos(folder, 't')
        self.book_Z_histos(folder)
        self.book_H_histos(folder)

    def leg3_id(self, row):
        return bool(row.mPFIDTight) and selections.muIsoLoose(row, 'm') 

    def leg4_id(self, row):
        return bool(row.tLooseIso3Hits) ##Why not tMediumMVAIso

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        #Z Selection
        if not selections.ZEESelection(row): return False
        if selections.overlap(row, 'e1','e2','m','t') : return False
        if not selections.signalTauSelection(row,'t'): return False
        if not bool(row.tAntiMuonTight2): return False
        if not bool(row.tAntiElectronLoose): return False
        #if row.LT < 45: return False
        if row.mPt + row.tPt < 35: return False
        return selections.signalMuonSelection(row,'m')

    def sign_cut(self, row):
        ''' Returns true if muon and tau are OS '''
        return not bool(row.m_t_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m') * \
            mcCorrectors.get_electron_corrections(row, 'e1', 'e2')

    def leg3_weight(self, row):
        return fr_fcn.mu_loose_jetpt_fr( row.mJetPt ) / (1 - fr_fcn.mu_loose_jetpt_fr( row.mJetPt ))

    def leg4_weight(self, row):
        return fr_fcn.tau_jetpt_fr( row.tJetPt ) / (1 -  fr_fcn.tau_jetpt_fr( row.tJetPt ))
