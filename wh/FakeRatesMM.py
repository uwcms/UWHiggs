'''

Measure fake rates in dimuon events.

We measure in QCD (anti-iso mu) and W+jet (iso mu) control regions.

The layout of output is:

    region/denom_tag/var1
    region/denom_tag/var2
    region/denom_tag/num_tag/var1
    region/denom_tag/num_tag/var2

Author: Evan K. Friis, UW

'''

from array import array
import MuMuTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import mcCorrectors
import baseSelections as selections
import os
import ROOT

def control_region(row):
    # Figure out what control region we are in.
    if row.m1RelPFIsoDB < 0.15 and row.m1MtToMET > 35 and row.m2MtToMET < 35:
        return 'wjets'
    elif row.m1RelPFIsoDB > 0.3 and row.type1_pfMetEt < 25:
        return 'qcd'
    else:
        return None


class FakeRatesMM(MegaBase):
    tree = 'mm/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(FakeRatesMM, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuMuTree.MuMuTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.lepIds = ['pfidiso02', 'h2taucuts', 'h2taucuts020'] #, 'h2taucuts025']
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')

    def begin(self):
        for region in ['wjets', 'qcd']:#, 'all']:
            for denom in ['pt10', 'pt20']: #, 'pt10b', 'pt10t', 'pt10f']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos
                
                #SPECIAL NTUPLE!
                denom_histos['muonInfo'] = self.book(
                    os.path.join(region, denom),
                    'muonInfo', "muonInfo", 
                    'muonPt:muonJetPt:muonJetCSVBtag:muonPVDXY:weight:'+':'.join(self.lepIds), 
                    type=ROOT.TNtuple)
                
                for numerator in self.lepIds:
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


                    book_histo('muonPt', 'Muon Pt', 200, 0, 200)
                    book_histo('muonJetPt', 'Muon Jet Pt', 200, 0, 200)
                    book_histo('muonAbsEta', 'Muon Abs Eta', 100, -2.5, 2.5)
                    book_histo('muonPtRatio', 'Muon Pt', 100, 0., 1.)
                    book_histo('muonPtDiff', 'Muon Pt', 200, 0., 200.)

                    book_histo('m1MtToMET', 'Muon 1 MT', 100, 0, 200)
                    book_histo('m2MtToMET', 'Muon 2 MT', 100, 0, 200)
                    book_histo('m1m2Mass', 'DiMuon Mass', 100, 0, 200)
                    book_histo("m1JetBtag", "Muon 2 Pt", 100, -10, 3.3)
                    book_histo("m2JetBtag", "Muon 2 Pt", 100, -10, 3.3)

                    book_histo('m2JetptD', "", 200, 0, 1)
                    book_histo('m2Jetmult', "", 50, 0, 50)


    def process(self):

        def preselection(row):
            double_mu_pass =  row.doubleMuPass and \
                row.m1MatchesDoubleMuPaths > 0 and \
                row.m2MatchesDoubleMuPaths > 0
            double_muTrk_pass = row.doubleMuTrkPass and \
                 row.m1MatchesDoubleMuTrkPaths > 0 and \
                 row.m2MatchesDoubleMuTrkPaths > 0
            if not ( double_mu_pass or double_muTrk_pass ): return False

            if not row.m1_m2_SS: return False
            if row.m2Pt > row.m1Pt: return False
            if row.m1_m2_Mass < 20: return False
            if not row.m1Pt > 20: return False
            if not row.m1PFIDTight: return False
            if not selections.muSelection(row, 'm1'):  return False
            if not selections.muSelection(row, 'm2'):  return False
            if not selections.vetos(row):             return False #applies mu bjet e additional tau vetoes
            if not (row.jetVeto40_DR05 >= 1): return False
            return True

        def fill(the_histos, row, fillNtuple=False):
            # Get PU weight - fix me
            weight = 1
            if row.run == 1:
                weight = self.pucorrector(row.nTruePU) * \
                         mcCorrectors.double_muon_trigger(row,'m1','m2')
                        
            the_histos['muonPt'].Fill(row.m2Pt, weight)
            the_histos['muonJetPt'].Fill(max(row.m2JetPt, row.m2Pt), weight)
            the_histos['muonAbsEta'].Fill(row.m2AbsEta, weight)
            the_histos['muonPtRatio'].Fill(row.m2Pt/max(row.m2JetPt, row.m2Pt), weight)
            the_histos['muonPtDiff'].Fill(max(row.m2JetPt, row.m2Pt) - row.m2Pt, weight)

            the_histos['m1MtToMET'].Fill(row.m1MtToMET, weight)
            the_histos['m1m2Mass' ].Fill(row.m1_m2_Mass, weight)
            the_histos['m2MtToMET'].Fill(row.m2MtToMET, weight)
            the_histos["m1JetBtag"].Fill(row.m1JetBtag, weight)
            the_histos["m2JetBtag"].Fill(row.m2JetBtag, weight)

            the_histos['m2JetptD'].Fill(row.m2JetptD, weight)
            the_histos['m2Jetmult'].Fill(row.m2Jetmult, weight)

            if fillNtuple:
                pfidiso02    = float( row.m2PFIDTight and row.m2RelPFIsoDB < 0.2)
                h2taucuts    = float( row.m2PFIDTight and ((row.m2RelPFIsoDB < 0.15 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.1 ))
                h2taucuts020 = float( row.m2PFIDTight and ((row.m2RelPFIsoDB < 0.20 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.15))
                the_histos['muonInfo'].Fill( array("f", [row.m2Pt, row.m2JetPt, max(0, row.m2JetCSVBtag), abs(row.m2PVDXY), weight, pfidiso02, h2taucuts, h2taucuts020] ) )
                
        histos = self.histograms
        for row in self.tree:
            if not preselection(row):
                continue
            region = control_region(row)
            if region is None:
                continue

            def fill_region(region, tag):
                # This is a QCD or Wjets
                fill(histos[(region, tag)], row, True)

                if row.m2PFIDTight:
                    if row.m2RelPFIsoDB < 0.2:
                        fill(histos[(region, tag, 'pfidiso02')], row)

                    if (row.m2RelPFIsoDB < 0.15 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.1:
                        fill(histos[(region, tag, 'h2taucuts')], row)

                    if (row.m2RelPFIsoDB < 0.20 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.15:
                        fill(histos[(region, tag, 'h2taucuts020')], row)
                    
                    #if (row.m2RelPFIsoDB < 0.25 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.20:
                    #    fill(histos[(region, tag, 'h2taucuts025')], row)
                        
            fill_region(region, 'pt10')
            #fill_region('all', 'pt10')

            if row.m2Pt > 20:
                fill_region(region, 'pt20')
                #fill_region('all', 'pt20')

    def finish(self):
        self.write_histos()
