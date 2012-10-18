'''
Base analyzer for muon id fake-rate estimation: Z->e e
'''

import EMUFakeRatesBase
from EEETauTree import EEETauTree
import baseSelections as selections

class EFakeRateEEET(EMUFakeRatesBase.EMUFakeRatesBase):
    tree = 'eeet/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        #print 'Called: EFakeRateEEET'
        super(EFakeRateEEET, self).__init__(tree, outfile, EEETauTree, **kwargs)
        self.lepton   = 'electron'
        self.branchId = 'e3'
        
    def zSelection(self, row):
        if not selections.ZEESelectionNoVetos(row): return False
        if selections.overlap(row, 'e1','e2','e3','t') : return False
        if bool(row.eVetoMVAIso):       return False
        return selections.signalElectronSelection(row,'e3')

    def lepton_passes_tight_iso(self, row):
        return bool(row.e3RelPFIsoDB < 0.10) ##THIS SEEMS too low        

    def lepton_passes_loose_iso(self, row):
        return bool(row.e3RelPFIsoDB < 0.25) ##THIS SEEMS too low        

    
