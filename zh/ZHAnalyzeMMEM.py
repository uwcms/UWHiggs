'''

Analyze MMEM events for the ZH analysis

'''

#from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import glob
from EMuMuMuTree import EMuMuMuTree
import os
import baseSelections as selections
import mcCorrectors
import ZHAnalyzerBase
import ROOT
import fake_rate_functions as fr_fcn

################################################################################
#### Analysis logic ############################################################
################################################################################

class ZHAnalyzeMMEM(ZHAnalyzerBase.ZHAnalyzerBase):
    tree = 'emmm/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ZHAnalyzeMMEM, self).__init__(tree, outfile, EMuMuMuTree, 'EM', **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')
        target = os.environ['megatarget']
        ## if 'HWW3l' in target:
        ##     print "HACK using S6 PU weights for HWW3l"
        ##     mcCorrectors.force_pu_distribution('S6')

    def Z_decay_products(self):
        return ('m1','m2')

    def H_decay_products(self):
        return ('e','m3')

    def book_histos(self, folder):
        self.book_general_histos(folder)
        self.book_kin_histos(folder, 'm1')
        self.book_kin_histos(folder, 'm2')
        self.book_kin_histos(folder, 'e')
        self.book_kin_histos(folder, 'm3')
        self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)
        self.book_Z_histos(folder)
        self.book_H_histos(folder)

    def probe1_id(self, row):
        return selections.eleID(row, 'e') and selections.elIsoLoose(row, 'e')

    def probe2_id(self, row):
        return bool(row.m3PFIDTight) and selections.muIsoLoose(row, 'm3')

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not selections.ZMuMuSelection(row): return False
        if selections.overlap(row, 'm1','m2','e','m3') : return False
        if not selections.signalMuonSelection(row,'m3'): return False
        if not selections.signalElectronSelection(row,'e'): return False
        return True

    def sign_cut(self, row):
        ''' Returns true if e and mu are OS '''
        return not bool(row.e_m3_SS)

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m1','m2', 'm3') * \
            mcCorrectors.double_muon_trigger(row,'m1','m2')

    def obj1_weight(self, row):
        return fr_fcn.e_loose_fr( row.ePt ) / (1 - fr_fcn.e_loose_fr( row.ePt ))

    def obj2_weight(self, row):
        return fr_fcn.mu_loose_fr( row.m3Pt) / (1 -  fr_fcn.mu_loose_fr( row.m3Pt));
