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
        if 'HWW3l' in target:
            print "HACK using S6 PU weights for HWW3l"
            mcCorrectors.force_pu_distribution('S6')

    def book_histos(self, folder):
        super(ZHAnalyzeMMTT, self).book_general_histos(folder)
        super(ZHAnalyzeMMTT, self).book_kin_histos(folder, 'm1')
        super(ZHAnalyzeMMTT, self).book_kin_histos(folder, 'm2')
        super(ZHAnalyzeMMTT, self).book_kin_histos(folder, 't1')
        super(ZHAnalyzeMMTT, self).book_kin_histos(folder, 't2')
        super(ZHAnalyzeMMTT, self).book_mass_histos(folder, 'm1','m2','t1','t2')
        self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)

    def probe1_id(self, row):
        return bool(row.t1MediumIso) ##THIS SEEMS too low

    def probe2_id(self, row):
        return bool(row.t2MediumIso) ##SHOULD BE TIGHT!!!

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not selections.ZMuMuSelection(row): return False
        if not selections.overlap(row, 'm1','m2','t1','t2') : return False
        if not selections.signalTauSelection(row,'t1'): return False
        if not selections.signalTauSelection(row,'t2'): return False
        if bool(row.t1AntiMuonTight): return False
        if bool(row.t1AntiElectronMedium): return False
        if bool(row.t2AntiMuonTight): return False
        if bool(row.t2AntiElectronMedium): return False
        return True

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return not bool(row.t1_t2_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return meCorrectors.pu_corrector(row.nTruePU) * \
            meCorrectors.get_muon_corrections(row,'m1','m2') * \
            double_muon_trigger(row,'m1','m2')

    def obj1_weight(self, row):
        return fr_fcn.tau_fr(max(row.t1JetPt, row.t1Pt))
        #return highpt_mu_fr(row.m1Pt)

    def obj2_weight(self, row):
        return fr_fcn.tau_fr(max(row.t2JetPt, row.t2Pt))
        #return lowpt_mu_fr(row.m2Pt)
