'''

Analyze MMMT events for the ZH analysis

'''

#from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import glob
from MuMuMuTauTree import MuMuMuTauTree
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

class ZHAnalyzeMMMT(ZHAnalyzerBase.ZHAnalyzerBase):
    tree = 'mmmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ZHAnalyzeMMMT, self).__init__(tree, outfile, MuMuMuTauTree, "MT", **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        target = os.environ['megatarget']
        if 'HWW3l' in target:
            print "HACK using S6 PU weights for HWW3l"
            mcCorrectors.force_pu_distribution('S6')

##     def get_channel(self):
##         return 'MT'
    @staticmethod
    def Z_decay_products():
        return ('m1','m2')

    @staticmethod
    def H_decay_products():
        return ('m3','t')

    def book_histos(self, folder):
        self.book_general_histos(folder)
        self.book_kin_histos(folder, 'm1')
        self.book_kin_histos(folder, 'm2')
        self.book_kin_histos(folder, 'm3')
        self.book_kin_histos(folder, 't')
        self.book_Z_histos(folder)
        self.book_H_histos(folder)
        self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)

    def probe1_id(self, row):
        return bool(row.m3RelPFIsoDB < 0.15) ##THIS SEEMS too low

    def probe2_id(self, row):
        return bool(row.tMediumIso) ##Why not tMediumMVAIso

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not selections.ZMuMuSelection(row): return False
        if selections.overlap(row, 'm1','m2','m3','t') : return False
        if not selections.signalTauSelection(row,'t'): return False
        if not bool(row.tAntiMuonTight): return False
        if not bool(row.tAntiElectronLoose): return False
        return selections.signalMuonSelection(row,'m3')

    def sign_cut(self, row):
        ''' Returns true if the probes are OS '''
        return not bool(row.m3_t_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return mcCorrectors.pu_corrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m1','m2','m3') * \
            mcCorrectors.double_muon_trigger(row,'m1','m2')

    def obj1_weight(self, row):
        return fr_fcn.mu_fr(max(row.m3JetPt, row.m3Pt))
        #return highpt_mu_fr(row.m1Pt)

    def obj2_weight(self, row):
        return fr_fcn.tau_fr(max(row.tJetPt, row.tPt))
        #return lowpt_mu_fr(row.m2Pt)


if __name__ == "__main__":
    import pprint
    pprint.pprint(ZHAnalyzeMMMT.build_zh_folder_structure())
