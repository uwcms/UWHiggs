'''
Base analyzer for muon id fake-rate estimation: Z->mu mu
'''

import EMUFakeRatesBase
import baseSelections as selections
from MuMuMuTauTree import MuMuMuTauTree

class MUFakeRateMMMT(EMUFakeRatesBase.EMUFakeRatesBase):
    tree = 'mmmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(MUFakeRateMMMT, self).__init__(tree, outfile, MuMuMuTauTree, **kwargs)
        self.lepton   = 'muon'
        self.branchId = 'm3'
        
    def zSelection(self, row):
        if not selections.ZMuMuSelectionNoVetos(row): return False
        if bool(row.muGlbIsoVetoPt10):  return False
        if selections.overlap(row, 'm1','m2','m3','t') : return False
        return selections.signalMuonSelection(row,'m3') #Include PFIdTight

    def lepton_passes_iso(self, row):
        return bool(row.m3RelPFIsoDB < 0.15)# and bool(getattr(row, 'm3PFIDTight') ) ##THIS SEEMS too low        

