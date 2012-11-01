'''
Base analyzer for muon id fake-rate estimation: Z->e e
'''

import EMUFakeRatesBase
from EEETree import EEETree
import baseSelections as selections

class EFakeRateEEET(EMUFakeRatesBase.EMUFakeRatesBase):
    tree = 'eee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        #print 'Called: EFakeRateEEET'
        super(EFakeRateEEET, self).__init__(tree, outfile, EEETree, **kwargs)
        self.lepton   = 'electron'
        self.branchId = 'e3'
        
    def zSelection(self, row):
        if not selections.ZEESelectionNoVetos(row): return False
        if selections.overlap(row, 'e1','e2','e3') : return False
        if bool(row.eVetoMVAIso):       return False
        if bool(row.bjetCSVVeto):       return False            
        ## if row.e3Pt     < 10:            return False
        ## if row.e3AbsEta > 2.5:           return False
        ## if row.e3DZ     > 0.1:           return False
        ## return True
        return selections.signalElectronSelection(row,'e3')

    def lepton_passes_tight_iso(self, row):
        return bool(row.e3RelPFIsoDB < 0.10) and bool( row.e3MVAIDH2TauWP ) ##THIS SEEMS too low        

    def lepton_passes_loose_iso(self, row):
        return bool(row.e3RelPFIsoDB < 0.25) and bool( row.e3MVAIDH2TauWP ) ##THIS SEEMS too low        

    
