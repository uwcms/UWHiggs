'''
Base analyzer for hadronic tau fake-rate estimation
'''

from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import os

class EMUFakeRatesBase(MegaBase):
    def __init__(self, tree, outfile, wrapper, **kwargs):
        #print 'called: TauFakeRatesBase'
        super(EMUFakeRatesBase, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = wrapper(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.lepton   = ''#self.leptonName()
        self.branchId = ''#self.branchIdentifier()

    def begin(self):
        # Book histograms
        region = 'zlt'
        for denom in ['pt10']:
            denom_key = (region, denom)
            denom_histos = {}
            self.histograms[denom_key] = denom_histos

            for numerator in ['idPassed']:
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

                book_histo(self.lepton+'Pt',     self.lepton+' Pt', 100, 0, 100)
                book_histo(self.lepton+'JetPt',  self.lepton+' Jet Pt', 100, 0, 100)
                book_histo(self.lepton+'AbsEta', self.lepton+' Abs Eta', 100, -2.5, 2.5)
    
    def process(self):
        # Generic filler function to fill plots after selection
        def fill(the_histos, row):
            weight = 1.0
            the_histos[self.lepton+'Pt'].Fill(    getattr( row, self.branchId+'Pt'    ), weight)
            the_histos[self.lepton+'JetPt'].Fill( getattr( row, self.branchId+'JetPt' ), weight)
            the_histos[self.lepton+'AbsEta'].Fill(getattr( row, self.branchId+'AbsEta'), weight)

        def preselection(self, row):
            if not self.zSelection(row):    return False
                #print 'Passed Z Selection'
            if not row.metEt > 20:          return False
                #print 'Passed MET Selection'
            if row.tPt  < 20:               return False
            if row.tAbsEta  > 2.3:          return False
            if abs(row.tDZ) > 0.2:          return False    
                #print 'Passed TAU Selection'
            if not bool(row.tDecayFinding): return False
            if not bool(row.tLooseIso):     return False
                #print 'Passed TAU Isolation'
            ## if not row.tAntiElectronMedium: return False #in the AN is not specified but seems reasonable
            ## if not row.tAntiMuonTight:      return False #in the AN is not specified but seems reasonable
            if not getattr(row,self.branchId+'MtToMET') > 30:           return False
                #print 'Passed MT Isolation'
            #Vetos
            ## if bool(row.muGlbIsoVetoPt10): return False
            ## if bool(row.tauVetoPt20):      return False
            ## if bool(row.eVetoMVAIso):      return False
            return True

        histos = self.histograms
        # Denominator regions (dicts of histograms)
        pt10 = histos[('zlt', 'pt10')]

        # Analyze data.  Select events with a good Z.
        for row in self.tree:
            if not preselection(self, row):
                continue

            # Fill denominator
            #print 'PRESELECTION PASSED!'
            fill(pt10, row)
            if self.lepton_passes_iso(row):
                #print 'ISOLATION PASSED!'
                fill(histos[('zlt', 'pt10', 'idPassed')], row)


    def finish(self):
        self.write_histos()

