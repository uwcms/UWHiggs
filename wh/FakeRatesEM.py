'''

Measure fake rates in the E+Mu channel

We measure in QCD (anti-iso mu) and W+jet (iso mu) control regions.

The layout of output is:
    region/denom_tag/var1
    region/denom_tag/var2
    region/denom_tag/num_tag/var1
    region/denom_tag/num_tag/var2

Author: Evan K. Friis, UW

'''

import EMuTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import baseSelections as selections
import os

def control_region(row):
    # Figure out what control region we are in.
    if row.mRelPFIsoDB < 0.15 and row.mMtToMET > 40 and row.eMtToMET < 30:
        return 'wjets'
    elif row.mRelPFIsoDB > 0.3 and row.type1_pfMetEt < 25:
        return 'qcd'
    else:
        return None

class FakeRatesEM(MegaBase):
    tree = 'em/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(FakeRatesEM, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EMuTree.EMuTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):
        for region in ['wjets', 'qcd']:
            for denom in ['pt10','pt20']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos

                for numerator in ['id', 'iso03', 'idiso03',
                                  'idiso02', 'idiso01', 'h2taucuts',
                                  'h2taucuts020', 'h2taucuts025',
                                  'eid13idiso02', 'eid13h2taucuts', 'eid13h2taucuts020',
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

                    book_histo('ePt', 'e Pt', 100, 0, 100)
                    book_histo('eJetPt', 'e Jet Pt', 100, 0, 100)
                    book_histo('eAbsEta', 'e Abs Eta', 100, -2.5, 2.5)
                    #book_histo('metSignificance', 'MET sig.', 100, 0, 10)
                    book_histo('mMtToMET', 'm MT', 100, 0, 200)

    def process(self):

        def preselection(row):
            if not row.e_m_SS: return False
            if not selections.muSelection(row, 'm'):  return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
            if not selections.eSelection(row, 'e'):   return False #applies basic selection (eta, pt > 10, DZ, missingHits, jetBTag, HasConversion and chargedIdTight)
            if not selections.vetos(row):             return False #applies mu bjet e additional tau vetoes
            return True
        #if self.is7TeV:
            #base_selection = 'mu17ele8Pass && ' + base_selection

        def fill(the_histos, row):
            the_histos['ePt'].Fill(row.ePt)
            the_histos['eJetPt'].Fill(max(row.eJetPt, row.ePt))
            the_histos['eAbsEta'].Fill(row.eAbsEta)
            #the_histos['metSignificance'].Fill(row.metSignificance)
            the_histos['mMtToMET'].Fill(row.mMtToMET)

        def fill_region(region,pt_cut):
            fill(histos[(region, pt_cut)], row)

            if row.eRelPFIsoDB < 0.3:
                fill(histos[(region, pt_cut, 'iso03')], row)

            if row.eMVAIDH2TauWP:
                fill(histos[(region, pt_cut, 'id')], row)
                if row.eRelPFIsoDB < 0.3:
                    fill(histos[(region, pt_cut, 'idiso03')], row)

                if row.eRelPFIsoDB < 0.2:
                    fill(histos[(region, pt_cut, 'idiso02')], row)

                if row.eRelPFIsoDB < 0.1:
                    fill(histos[(region, pt_cut, 'idiso01')], row)

                if (row.eRelPFIsoDB < 0.15 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.1:
                    fill(histos[(region, pt_cut, 'h2taucuts')], row)

                if (row.eRelPFIsoDB < 0.20 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.15:
                    fill(histos[(region, pt_cut, 'h2taucuts020')], row)

                if (row.eRelPFIsoDB < 0.25 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.20:
                    fill(histos[(region, pt_cut, 'h2taucuts025')], row)
                    
            if selections.summer_2013_eid(row, 'e'):
                if row.eRelPFIsoDB < 0.2:
                    fill(histos[(region, pt_cut, 'eid13idiso02')], row)

                if (row.eRelPFIsoDB < 0.15 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.1:
                    fill(histos[(region, pt_cut, 'eid13h2taucuts')], row)

                if (row.eRelPFIsoDB < 0.20 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.15:
                    fill(histos[(region, pt_cut, 'eid13h2taucuts020')], row)


        histos = self.histograms
        for row in self.tree:
            if not preselection(row):
                continue
            region = control_region(row)
            if region is None:
                continue
            # This is a QCD or Wjets
            is7TeV = bool('7TeV' in os.environ['jobid'])
            use_iso_trigger = not is7TeV
            mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
            mu8e17 = (row.mu8ele17isoPass and row.ePt >= 20) #if use_iso_trigger else (row.mu17ele8Pass and row.mPt < 20)
            if mu17e8:
                fill_region(region,'pt10')
            if mu8e17:
                fill_region(region,'pt20')

    def finish(self):
        self.write_histos()
