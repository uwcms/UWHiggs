'''

Analyze EEEM events for the ZH analysis

'''

#from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import glob
from EEEMuTree import EEEMuTree
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

class ZHAnalyzeEEEM(ZHAnalyzerBase.ZHAnalyzerBase):
    tree = 'eeem/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ZHAnalyzeEEEM, self).__init__(tree, outfile, EEEMuTree, 'EM', **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        target = os.environ['megatarget']
        self.pucorrector = mcCorrectors.make_puCorrector('doublee')
        ## if 'HWW3l' in target:
        ##     print "HACK using S6 PU weights for HWW3l"
        ##     mcCorrectors.force_pu_distribution('S6')

    def Z_decay_products(self):
        return ('e1','e2')

    def H_decay_products(self):
        return ('e3','m')

    def book_histos(self, folder):
        self.book_general_histos(folder)
        self.book_kin_histos(folder, 'e1')
        self.book_kin_histos(folder, 'e2')
        self.book_kin_histos(folder, 'm')
        self.book_kin_histos(folder, 'e3')
        self.book_Z_histos(folder)
        self.book_H_histos(folder)

    def probe1_id(self, row):
        return selections.eleID(row, 'e3') and selections.elIsoLoose(row, 'e3')

    def probe2_id(self, row):
        return bool(row.mPFIDTight) and selections.muIsoLoose(row, 'm')

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not selections.ZEESelection(row): return False
        if selections.overlap(row, 'e1','e2','e3','m') : return False
        if not selections.signalMuonSelection(row,'m'): return False
        if not selections.signalElectronSelection(row,'e3'): return False
        return True

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return not bool(row.e3_m_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m') * \
            mcCorrectors.get_electron_corrections(row, 'e1','e2')

    def obj1_weight(self, row):
        return fr_fcn.e_loose_fr( row.e3Pt ) / (1 - fr_fcn.e_loose_fr( row.e3Pt ))

    def obj2_weight(self, row):
        return fr_fcn.mu_loose_fr( row.mPt) / (1 - fr_fcn.mu_loose_fr(row.mPt))
