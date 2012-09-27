'''

Analyze EEET events for the ZH analysis

'''

#from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import glob
from EEETauTree import EEETauTree
import os
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import mcCorrectors
import ZHAnalyzerBase
import baseSelections as selections
import ROOT
import fake_rate_functions as fr_fcn

################################################################################
#### Analysis logic ############################################################
################################################################################

class ZHAnalyzeEEET(ZHAnalyzerBase.ZHAnalyzerBase):
    tree = 'eeet/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ZHAnalyzeEEET, self).__init__(tree, outfile, EEETauTree, 'ET', **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        target = os.environ['megatarget']
        if 'HWW3l' in target:
            print "HACK using S6 PU weights for HWW3l"
            mcCorrectors.force_pu_distribution('S6')

    def book_histos(self, folder):
        super(ZHAnalyzeEEET, self).book_general_histos(folder)
        super(ZHAnalyzeEEET, self).book_kin_histos(folder, 'e1')
        super(ZHAnalyzeEEET, self).book_kin_histos(folder, 'e2')
        super(ZHAnalyzeEEET, self).book_kin_histos(folder, 'm')
        super(ZHAnalyzeEEET, self).book_kin_histos(folder, 't')
        super(ZHAnalyzeEEET, self).book_mass_histos(folder, 'e1','e2','m','t')

    def probe1_id(self, row):
        return bool(row.eRelPFIsoDB < 0.10) ##THIS SEEMS too low

    def probe2_id(self, row):
        return bool(row.tMediumIso) ##Why not tMediumMVAIso

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        #Z Selection
        if not selections.ZEESelection(row): return False
        if not selections.overlap(row, 'e1','e2','e','t') : return False
        if not selections.signalTauSelection(row,'t'): return False
        if bool(row.tAntiMuonLoose): return False
        if bool(row.tAntiElectronMVA): return False
        if bool(row.e_t_SS): return False
        return selections.signalElectronSelection(row,'e')

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return not bool(row.e3_t_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return meCorrectors.pu_corrector(row.nTruePU) * \
            get_electron_corrections(row, 'e1','e2','e3')

    def obj1_weight(self, row):
        return fr_fcn.e_fr(max(row.e3JetPt, row.e3Pt))
        #return highpt_mu_fr(row.m1Pt)

    def obj2_weight(self, row):
        return fr_fcn.tau_fr(max(row.tJetPt, row.tPt))
        #return lowpt_mu_fr(row.m2Pt)
