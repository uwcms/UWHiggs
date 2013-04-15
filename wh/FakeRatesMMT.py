'''

Measure tau fake rates in Z->mumu + tau jet events

Author: Evan K. Friis, UW Madison

The layout of output is:

    zmm/denom_tag/var1
    zmm/denom_tag/var2
    zmm/denom_tag/num_tag/var1
    zmm/denom_tag/num_tag/var2

'''

import MuMuTauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import os

class FakeRatesMMT(MegaBase):
    tree = 'mmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(FakeRatesMMT, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuMuTauTree.MuMuTauTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):
        # Book histograms
        region = 'ztt'
        for denom in ['pt20']:
            denom_key = (region, denom)
            denom_histos = {}
            self.histograms[denom_key] = denom_histos

            for numerator in ['mvaloose',
                              'hpsloose',
                              'mvamedium',
                              'hpsmedium',
                              ]:
                num_key = (region, denom, numerator)
                num_histos = {}
                self.histograms[num_key] = num_histos

                def book_histo(name, *args):
                    # Helper to book a histogram
                    if name not in denom_histos:
                        denom_histos[name] = self.book(os.path.join(
                            region, denom), name, *args)
                    num_histos[name] = self.book(os.path.join(
                        region, denom, numerator), name, *args)

                book_histo('tauPt', 'Tau Pt', 100, 0, 100)
                book_histo('tauJetPt', 'Tau Jet Pt', 100, 0, 100)
                book_histo('tauAbsEta', 'Tau Abs Eta', 100, -2.5, 2.5)

    def process(self):
        # Generic filler function to fill plots after selection
        def fill(the_histos, row):
            weight = 1.0
            the_histos['tauPt'].Fill(row.tPt, weight)
            the_histos['tauJetPt'].Fill(row.tJetPt, weight)
            the_histos['tauAbsEta'].Fill(row.tAbsEta, weight)

        def preselection(row):
            if not row.m1RelPFIsoDB < 0.25: return False
            if not row.m2RelPFIsoDB < 0.25: return False
            if not row.m1PFIDTight: return False
            if not row.m2PFIDTight: return False
            if not abs(row.m1_m2_Mass-91.2) < 10: return False
            if not row.tPt > 20: return False
            if not row.tAbsEta < 2.3: return False
            if row.tMuOverlap: return False
            if row.tCiCTightElecOverlap: return False
            if not row.tAntiMuonTight: return False
            if not row.tAntiElectronMVA3Tight: return False
            if not row.tMtToMET < 30: return False
            return True

        histos = self.histograms

        # Denominator regions (dicts of histograms)
        pt20 = histos[('ztt', 'pt20')]

        # Analyze data.  Select events with a good Z.
        for row in self.tree:
            if not preselection(row):
                continue

            # Fill denominator
            fill(pt20, row)

            # Fill numerators
            if row.tLooseIso3Hits:
                fill(histos[('ztt', 'pt20', 'hpsloose')], row)
            if row.tLooseMVAIso:
                fill(histos[('ztt', 'pt20', 'mvaloose')], row)
            if row.tMediumIso3Hits:
                fill(histos[('ztt', 'pt20', 'hpsmedium')], row)
            if row.tMediumMVAIso:
                fill(histos[('ztt', 'pt20', 'mvamedium')], row)


    def finish(self):
        self.write_histos()
