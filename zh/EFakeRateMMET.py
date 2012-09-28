'''
Base analyzer for muon id fake-rate estimation: Z->e e
'''

import EMUFakeRatesBase
import baseSelections as selections
from MuMuETauTree import MuMuETauTree

class EFakeRateMMET(EMUFakeRatesBase.EMUFakeRatesBase):
    tree = 'emmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(EFakeRateMMET, self).__init__(tree, outfile, MuMuETauTree, **kwargs)
        self.lepton   = 'electron'
        self.branchId = 'e'
        
    def zSelection(self, row):
        if not selections.ZMuMuSelection(row): return False
        if selections.overlap(row, 'm1','m2','e','t') : return False
        if not selections.signalTauSelection(row,'t'): return False
        if not bool(row.tLooseIso): return False
        return selections.signalElectronSelection(row,'e')
    
    def lepton_passes_iso(self, row):
        return bool(row.eRelPFIsoDB < 0.15) ##THIS SEEMS too low        
