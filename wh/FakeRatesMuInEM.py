'''

Measure *muon* fake rates in the E+Mu channel

We measure in QCD (anti-iso e) and W+jet (iso e + MT) control regions.

The layout of output is:
    region/denom_tag/var1
    region/denom_tag/var2
    region/denom_tag/num_tag/var1
    region/denom_tag/num_tag/var2

Author: Evan K. Friis, UW

'''

from array import array
import EMuTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import os
import ROOT

def control_region(row):
    # Figure out what control region we are in.
    if row.eRelPFIsoDB < 0.15 and row.eMtToMET > 40 and row.mMtToMET < 30:
        return 'wejets'
    elif row.eRelPFIsoDB > 0.3 and row.metEt < 25:
        return 'qcd'
    else:
        return None

class FakeRatesMuInEM(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(FakeRatesMuInEM, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EMuTree.EMuTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):
        for region in ['wejets', 'qcd']:
            for denom in ['pt10', 'pt20', 'pt10b', 'pt10t', 'pt10f']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos

                for numerator in ['pfid', 'iso03', 'pfidiso03',
                                  'pfidiso02', 'pfidiso01']:
                    num_key = (region, denom, numerator)
                    num_histos = {}
                    self.histograms[num_key] = num_histos

                    def book_histo(name, *args, **kwargs):
                        # Helper to book a histogram
                        if name not in denom_histos:
                            denom_histos[name] = self.book(os.path.join(
                                region, denom), name, *args, **kwargs)
                        num_histos[name] = self.book(os.path.join(
                            region, denom, numerator), name, *args, **kwargs)

                    mu_binning = array(
                        'd', [10, 12, 15, 20, 30, 50, 100, 200])
                    jet_binning = array(
                        'd', [10, 12, 15, 20, 25, 30, 50, 100, 200])
                    eta_binning = array(
                        'd', [0, 0.25, 0.8, 1.3, 1.8, 2.4])

                    book_histo('muonJetVsLeptonPt', 'Muon Pt',
                               len(jet_binning)-1, jet_binning,
                               len(mu_binning)-1, mu_binning,
                               type=ROOT.TH2F)
                    book_histo('muonJetVsEta', 'Muon Pt',
                               len(jet_binning)-1, jet_binning,
                               len(eta_binning)-1, eta_binning,
                               type=ROOT.TH2F)

                    book_histo('muonPt', 'Muon Pt', 16, 10, 50)
                    book_histo('muonJetPt', 'Muon Jet Pt', 200, 0, 200)
                    book_histo('muonAbsEta', 'Muon Abs Eta', 100, -2.5, 2.5)
                    book_histo('metSignificance', 'MET sig.', 100, 0, 10)
                    book_histo('m1MtToMET', 'Muon 1 MT', 100, 0, 200)


    def process(self):

        def preselection(row):
            if not row.e_m_SS: return False
            #if not row.mu8ele17Pass: return False
            if not row.mPt > 10: return False
            if not row.ePt > 20: return False
            if not row.eMVAIDH2TauWP: return False
            if not row.mAbsEta < 2.4: return False
            if not row.eAbsEta < 2.5: return False
            if not row.mJetBtag < 3.3: return False
            if not row.mPixHits: return False
            if row.eVetoCicTightIso: return False
            if row.muVetoPt5: return False
            if row.bjetCSVVeto: return False
            if row.tauVetoPt20: return False
            if not abs(row.mDZ) < 0.2: return False
            if not abs(row.eDZ) < 0.2: return False
            return True

        def fill(the_histos, row):
            # Get PU weight - fix me
            weight = 1
            the_histos['muonPt'].Fill(row.mPt, weight)
            the_histos['muonJetPt'].Fill(max(row.mJetPt, row.mPt), weight)
            the_histos['muonJetVsLeptonPt'].Fill(max(row.mJetPt, row.mPt), row.mPt, weight)
            the_histos['muonJetVsEta'].Fill(max(row.mJetPt, row.mPt), row.mAbsEta, weight)
            the_histos['muonAbsEta'].Fill(row.mAbsEta, weight)
            the_histos['metSignificance'].Fill(row.metSignificance, weight)

        histos = self.histograms
        for row in self.tree:
            if not preselection(row):
                continue
            region = control_region(row)
            if region is None:
                continue

            def fill_denominator(denominator_tag):
                # This is a QCD or wejets
                fill(histos[(region, denominator_tag)], row)

                if row.mPFIDTight:
                    fill(histos[(region, denominator_tag, 'pfid')], row)

                if row.mRelPFIsoDB < 0.3:
                    fill(histos[(region, denominator_tag, 'iso03')], row)

                if row.mPFIDTight and row.mRelPFIsoDB < 0.3:
                    fill(histos[(region, denominator_tag, 'pfidiso03')], row)

                if row.mPFIDTight and row.mRelPFIsoDB < 0.2:
                    fill(histos[(region, denominator_tag, 'pfidiso02')], row)

                if row.mPFIDTight and row.mRelPFIsoDB < 0.1:
                    fill(histos[(region, denominator_tag, 'pfidiso01')], row)

            fill_denominator('pt10')

            # Barrel/forward
            if row.mAbsEta < 0.8:
                fill_denominator('pt10b')
            elif row.mAbsEta < 1.3:
                fill_denominator('pt10t')
            else:
                fill_denominator('pt10f')

            if row.mPt > 20:
                fill_denominator('pt20')

    def finish(self):
        self.write_histos()
