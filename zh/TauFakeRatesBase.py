'''
Base analyzer for hadronic tau fake-rate estimation
'''

from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import os

class TauFakeRatesBase(MegaBase):
    def __init__(self, tree, outfile, wrapper, **kwargs):
        super(TauFakeRatesBase, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = wrapper(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.numerators = ['LooseMVAIso', 'LooseIso', 'MediumMVAIso', 'MediumIso', 'TightIso', 'TightMVAIso']

    def begin(self):
        # Book histograms
        region = 'ztt'
        for denom in ['pt10']:
            denom_key = (region, denom)
            denom_histos = {}
            self.histograms[denom_key] = denom_histos

            for numerator in self.numerators:
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
                book_histo('tauTauInvMass', ';M_{#tau #tau} [GeV/c^{2}];Events/10 [GeV/c^{2}]', 51, 0., 255)

    def process(self):
        # Generic filler function to fill plots after selection
        def fill(the_histos, row, t):
            weight = 1.0
            the_histos['tauPt'].Fill(getattr( row,t+'Pt'), weight)
            the_histos['tauJetPt'].Fill(getattr( row,t+'JetPt'), weight)
            the_histos['tauAbsEta'].Fill(getattr( row,t+'AbsEta'), weight)

        def preselection(self, row):
            if not self.zSelection(row):     return False
            if not row.t1Pt > 10:            return False #in the analysis this value is 20, but the fake rate plot is shown/fitted up to 10
            if not row.t2Pt > 10:            return False
            if not row.t1AbsEta < 2.3:       return False
            if not row.t2AbsEta < 2.3:       return False
            if abs(row.t1DZ) > 0.2:          return False    
            if abs(row.t2DZ) > 0.2:          return False    
            if not bool(row.t1DecayFinding): return False
            if not bool(row.t2DecayFinding): return False
            if not bool(row.t1_t2_SS):       return False
            if not row.t1AntiElectronMedium: return False #in the AN is not specified but seems reasonable
            if not row.t2AntiElectronMedium: return False
            if not row.t1AntiMuonTight:      return False #in the AN is not specified but seems reasonable
            if not row.t2AntiMuonTight:      return False
            return True

        histos = self.histograms

        # Denominator regions (dicts of histograms)
        pt10 = histos[('ztt', 'pt10')]

        # Analyze data.  Select events with a good Z.
        for row in self.tree:
            if not preselection(self, row):
                continue

            # Fill denominator
            fill(pt10, row, 't1')
            fill(pt10, row, 't2')
            pt10['tauTauInvMass'].Fill( row.t1_t2_Mass, 1.)

            for t in ['t1','t2']:
                for num in self.numerators:
                    if bool( getattr(row,t+num) ):
                        fill(histos[('ztt', 'pt10', num)], row, t)


    def finish(self):
        self.write_histos()

