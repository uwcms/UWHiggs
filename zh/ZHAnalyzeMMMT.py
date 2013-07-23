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
    name = 1
    def __init__(self, tree, outfile, **kwargs):
        super(ZHAnalyzeMMMT, self).__init__(tree, outfile, MuMuMuTauTree, "MT", **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        target = os.environ['megatarget']
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')
        ## if 'HWW3l' in target:
        ##     print "HACK using S6 PU weights for HWW3l"
        ##     mcCorrectors.force_pu_distribution('S10')

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

        self.book(folder, 'Pt',"total pt",100,0,100) # stephane
        #self.book(folder, '"Pt/(m1_m2_Pt+m3_t_Pt)', "ZH system kinematic ratio",20,0,1) 

        self.book_Z_histos(folder)
        self.book_H_histos(folder)
        self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)

    def leg3_id(self, row):
        return bool(row.m3PFIDTight) and selections.muIsoLoose(row, 'm3') ##THIS SEEMS too low

    def leg4_id(self, row):
        return bool(row.tLooseIso3Hits) ##Why not tMediumMVAIso

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not selections.ZMuMuSelection(row): return False
        if selections.overlap(row, 'm1','m2','m3','t') : return False
        if not selections.signalTauSelection(row,'t'): return False
        if not bool(row.tAntiMuonTight2): return False
        if not bool(row.tAntiElectronLoose): return False
       # if row.LT < 45: return False
        if row.m3Pt + row.tPt < 35: return False
        return selections.signalMuonSelection(row,'m3')

    def sign_cut(self, row):
        ''' Returns true if the probes are OS '''
        return not bool(row.m3_t_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m1','m2','m3') * \
            mcCorrectors.double_muon_trigger(row,'m1','m2')

    def leg3_weight(self, row):
        return fr_fcn.mu_loose_jetpt_fr( row.m3JetPt)/(1 - fr_fcn.mu_loose_jetpt_fr( row.m3JetPt))

    def leg4_weight(self, row):
        return fr_fcn.tau_jetpt_fr(row.tJetPt)/(1 - fr_fcn.tau_jetpt_fr(row.tJetPt))


if __name__ == "__main__":
    import pprint
    pprint.pprint(ZHAnalyzeMMMT.build_zh_folder_structure())
