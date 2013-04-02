'''

Analyze MMMT events for the ZH analysis

'''

import glob
from MuMuTauTauTree import MuMuTauTauTree
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

class ZHAnalyzeMMTT(ZHAnalyzerBase.ZHAnalyzerBase):
    tree = 'mmtt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ZHAnalyzeMMTT, self).__init__(tree, outfile, MuMuTauTauTree, 'TT', **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        target = os.environ['megatarget']
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')
        ## if 'HWW3l' in target:
        ##     print "HACK using S6 PU weights for HWW3l"
        ##     mcCorrectors.force_pu_distribution('S10')

    def Z_decay_products(self):
        return ('m1','m2')

    def H_decay_products(self):
        return ('t1','t2')

    def book_histos(self, folder):
        self.book_general_histos(folder)
        self.book_kin_histos(folder, 'm1')
        self.book_kin_histos(folder, 'm2')
        self.book_kin_histos(folder, 't1')
        self.book_kin_histos(folder, 't2')
        self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)
        self.book(folder, "leadingTauPt", "Leading #tau_{h} p_{T};p_{T} (GeV);Counts", 100, 0, 100)
        self.hfunc['leadingTauPt'] = lambda row, weight: max( row.t1Pt, row.t2Pt ) #Will NOT be weighted!!
        self.book_Z_histos(folder)
        self.book_H_histos(folder)

    def probe1_id(self, row):
        return bool(row.t1TightIso) ##THIS SEEMS too low

    def probe2_id(self, row):
        return bool(row.t2TightIso) ##SHOULD BE TIGHT!!!

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not selections.ZMuMuSelection(row): return False
        if selections.overlap(row, 'm1','m2','t1','t2') : return False
        if not selections.signalTauSelection(row,'t1'): return False
        if not selections.signalTauSelection(row,'t2'): return False
        if not bool(row.t1AntiMuonTight): return False
        if not bool(row.t1AntiElectronMedium): return False
        if not bool(row.t2AntiMuonTight): return False
        if not bool(row.t2AntiElectronMedium): return False
        if row.t1Pt < row.t2Pt: return False #Avoid double counting
        return True

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return not bool(row.t1_t2_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m1','m2') * \
            mcCorrectors.double_muon_trigger(row,'m1','m2')

    def obj1_weight(self, row):
        return fr_fcn.tau_tight_fr( row.t1Pt ) / (1- fr_fcn.tau_tight_fr( row.t1Pt ))

    def obj2_weight(self, row):
        return fr_fcn.tau_tight_fr( row.t2Pt ) / (1 - fr_fcn.tau_tight_fr( row.t2Pt ))

    ## def dump(self, row):
    ##     'debugging / sync helper function'
