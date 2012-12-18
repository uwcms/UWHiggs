'''

Measure fake rates in dielectron events.

NB that the triggers we are using always have the same iso on the leading
and subleading legs.

We measure W+jet (iso electron) control regions.

The layout of output is:

    region/denom_tag/var1
    region/denom_tag/var2
    region/denom_tag/num_tag/var1
    region/denom_tag/num_tag/var2

Author: Evan K. Friis, UW

'''

import EETree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import baseSelections as selections
import math
import os
import ROOT
import array

barrelThr = 1.479
accThr    = 2.5

class ChargeFlipProbabilityEE(MegaBase):
    tree = 'ee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ChargeFlipProbabilityEE, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EETree.EETree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):
        # Charge mis-ID measurements
        binsEndcap = 5
        binsBarrel = 10
        etabins    = [(i*barrelThr/binsBarrel) for i in range(binsBarrel+1)]+[ (barrelThr+i*(accThr-barrelThr)/binsEndcap) for i in range(1,binsEndcap+1)] 
        ptbins     = [10.,40.,130.]
        self.book('charge', 'flipped_electrons', 'Flipped electrons distribution; |#eta|; p_{T} [GeV]',     len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
        self.book('charge', 'matched_electrons', 'Gen-matched electrons distribution; |#eta|; p_{T} [GeV]', len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
            
        #print self.histograms.keys()

    def process(self):

        def preselection(row):
            if not row.doubleEPass: return False
            if not row.e1Pt > 20: return False
            if not selections.eSelection(row, 'e1'): return False
            if not row.e1MVAIDH2TauWP: return False
            if not selections.eSelection(row, 'e2'): return False
            if not selections.vetos(row): return False
            if  row.e1_e2_SS and \
                row.e1RelPFIsoDB < 0.15 and \
                row.e1MtToMET > 50 and \
                row.metSignificance > 3: return False
            if  row.e1_e2_SS and \
                row.e1RelPFIsoDB > 0.3 and \
                row.metSignificance < 3: return False
            return (row.e1RelPFIsoDB < 0.1 and row.e2RelPFIsoDB < 0.1 and row.e2MVAIDH2TauWP and row.e1MVAIDH2TauWP and not any([ row.muVetoPt5, row.tauVetoPt20, row.eVetoCicTightIso,]))

        #def fill(the_histos, row):

        histos = self.histograms
        for row in self.tree:
            if not preselection(row):
                continue

            if 'Zjets' in os.environ['megatarget']:
                for el in ['e1','e2']:
                    pdgid = getattr(row,el+'GenPdgId')
                    histos['charge/matched_electrons'].Fill( getattr(row,el+'AbsEta'), getattr(row,el+'Pt') )
                    if pdgid == -999 and row.e1_e2_SS: #FIXME make the electrons matched without caring about the charge   #e- --> -11; e+ --> 11 MISMEASURED ELECTRON # getattr(row,el+'Charge')*(11) == 
                        histos['charge/flipped_electrons'].Fill( getattr(row,el+'AbsEta'), getattr(row,el+'Pt') )

    def finish(self):
        self.write_histos()
