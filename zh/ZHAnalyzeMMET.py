'''

Analyze MMET events for the ZH analysis

'''

#from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import glob
from MuMuETauTree import MuMuETauTree
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

class ZHAnalyzeMMET(ZHAnalyzerBase.ZHAnalyzerBase):
    tree = 'emmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ZHAnalyzeMMET, self).__init__(tree, outfile, MuMuETauTree, 'ET', **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')
        target = os.environ['megatarget']
        ## if 'HWW3l' in target:
        ##     print "HACK using S6 PU weights for HWW3l"
        ##     mcCorrectors.force_pu_distribution('S6')
            
    def Z_decay_products(self):
        return ('m1','m2')

    def H_decay_products(self):
        return ('e','t')

    def book_histos(self, folder):
        self.book_general_histos(folder)
        self.book_kin_histos(folder, 'm1')
        self.book_kin_histos(folder, 'm2')
        self.book_kin_histos(folder, 'e')
        self.book_kin_histos(folder, 't')
        self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)
        self.book_Z_histos(folder)
        self.book_H_histos(folder)

    def probe1_id(self, row):
        return selections.eleID(row, 'e') and selections.elIsoTight(row, 'e') ##THIS SEEMS too low

    def probe2_id(self, row):
        return bool(row.tMediumIso) ##Why not tMediumMVAIso

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        #Z Selection
        if not selections.ZMuMuSelection(row): return False
        if selections.overlap(row, 'm1','m2','e','t') : return False
        if not selections.signalTauSelection(row,'t'): return False
        if not bool(row.tAntiMuonLoose): return False
        if not bool(row.tAntiElectronMVA): return False
        return selections.signalElectronSelection(row,'e')

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return not bool(row.e_t_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m1','m2') * \
            mcCorrectors.get_electron_corrections(row, 'e') * \
            mcCorrectors.double_muon_trigger(row,'m1','m2')

    def obj1_weight(self, row):
        return fr_fcn.e_tight_jetpt_fr( row.eJetPt ) / (1 - fr_fcn.e_tight_jetpt_fr( row.eJetPt ))
    #return fr_fcn.e_tight_fr( row.ePt )

    def obj2_weight(self, row):
        return fr_fcn.tau_medium_fr( row.tPt ) / (1 - fr_fcn.tau_medium_fr( row.tPt ))
