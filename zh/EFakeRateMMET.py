'''
Base analyzer for muon id fake-rate estimation: Z->e e
'''

import EMUFakeRatesBase
import baseSelections as selections
from MuMuETree import MuMuETree

class EFakeRateMMET(EMUFakeRatesBase.EMUFakeRatesBase):
    tree = 'emm/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(EFakeRateMMET, self).__init__(tree, outfile, MuMuETree, **kwargs)
        self.lepton   = 'electron'
        self.branchId = 'e'
        
    def zSelection(self, row):
        if not selections.ZMuMuSelectionNoVetos(row): return False
        if selections.overlap(row, 'm1','m2','e') : return False
        if bool(row.eVetoMVAIso):       return False
        if bool(row.bjetCSVVeto):       return False
        ## if row.ePt     < 10:            return False
        ## if row.eAbsEta > 2.5:           return False
        ## if row.eDZ     > 0.1:           return False
        ## return True
        return selections.signalElectronSelection(row,'e')
    
    def lepton_passes_tight_iso(self, row):
        return bool(row.eRelPFIsoDB < 0.10) and bool( row.eMVAIDH2TauWP ) ##THIS SEEMS too low        

    def lepton_passes_loose_iso(self, row):
        return bool(row.eRelPFIsoDB < 0.25)  and bool( row.eMVAIDH2TauWP ) ##THIS SEEMS too low        
