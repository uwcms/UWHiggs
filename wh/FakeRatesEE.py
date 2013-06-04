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
from FinalStateAnalysis.PlotTools.decorators import decorator
#import ROOT.TMath as math
import os
import array
from pprint import pprint
import ROOT
 
class FakeRatesEE(MegaBase):
    tree = 'ee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(FakeRatesEE, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EETree.EETree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):
        for region in ['wjets', 'qcd', 'wjetsNoZmass']:
            for denom in ['pt10', 'pt20']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos

                for numerator in ['mvaid', 'iso03', 'mvaidiso03',
                                  'mvaidiso01','h2taucuts',
                                  'h2taucuts020', 'h2taucuts025',
                                  'mvaidiso02']:
                    num_key = (region, denom, numerator)
                    num_histos = {}
                    self.histograms[num_key] = num_histos

                    def book_histo(name, *args, **kwargs):
                        # Helper to book a histogram
                        if name not in denom_histos:
                            denom_histos[name] = self.book(os.path.join(
                                region, denom), name, *args,**kwargs)
                        num_histos[name] = self.book(os.path.join(
                            region, denom, numerator), name, *args,**kwargs)

                    book_histo('electronPt', 'e Pt', 100, 0, 100)
                    book_histo('electronJetPt', 'e Jet Pt', 100, 0, 100)
                    book_histo('electronAbsEta', 'e Abs Eta', 100, -2.5, 2.5)
                    #book_histo('metSignificance', 'MET sig.', 100, 0, 10)
                    book_histo('e1MtToMET', 'e1 MT', 100, 0, 200)
                    book_histo('e2MtToMET', 'e2 MT', 100, 0, 200)
                    book_histo('e1MtToMET_Z', 'e1 MT', 100, 0, 200)
                    book_histo('e2MtToMET_Z', 'e2 MT', 100, 0, 200)
                    book_histo('e1MtToMET_NoZ', 'e1 MT', 100, 0, 200)
                    book_histo('e2MtToMET_NoZ', 'e2 MT', 100, 0, 200)
                    book_histo('MET', 'MET', 100, 0, 400)
                    book_histo('MET_Z', 'MET, Z Events', 100, 0, 400)
                    book_histo('MET_NoZ', 'MET, No Z Events', 100, 0, 400)
                    book_histo('OSS', 'SS to OS', 2, 0, 2)
                    book_histo('OSS_Z', 'SS to OS, Z events', 2, 0, 2)
                    book_histo('e1e2Mass', 'DiElectron Mass', 100, 0, 200)
                    book_histo('doubleEPrescale', 'prescale', 10, -0.5, 9.5)

                    book_histo('e2JetArea'        , "", 100, 0, 1)
                    book_histo('e2JetEtaEtaMoment', "", 100, 0, 0.2)
                    book_histo('e2JetEtaPhiMoment', "", 100, 0, 0.1)
                    book_histo('e2JetEtaPhiSpread', "", 100, 0, 0.2)
                    book_histo('e2JetPhiPhiMoment', "", 100, 0, 0.1)
                    book_histo('e2JetptD', "", 200, 0, 1)
                    book_histo('e2Jetaxis1', "", 200, 0, 1)
                    book_histo('e2Jetaxis2', "", 200, 0, 1)
                    book_histo('e2Jetmult', "", 50, 0, 50)
                    book_histo('e2JetmultMLPQC', "", 50, 0, 50)
                    book_histo('e2JetmultMLP', "", 50, 0, 50)
                    book_histo('e2JetQGLikelihoodID', "", 200, 0, 1)
                    book_histo('e2JetQGMVAID', "", 200, 0, 1)

                    book_histo('e2JetQGLikelihoodIDvse2JetPt'," ", 100, 0,100, 200, 0, 1, type=ROOT.TH2F)
                    book_histo('e2JetQGMVAIDvse2JetPt'," ", 100, 0,100, 200, 0, 1, type=ROOT.TH2F)

                    book_histo('e2Jetmultvse2JetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)
                    book_histo('e2JetmultMLPvse2JetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)
                    book_histo('e2JetmultMLPQCvse2JetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)
                    
                    book_histo('e2Jetmultvse2JetPt', '', 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('e2JetmultMLPvse2JetPt', '',  200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('e2JetmultMLPQCvse2JetPt', '' , 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('e2JetptDvse2JetPt', '' ,  200, 0, 200, 200, 0, 1,type=ROOT.TH2F)



    def process(self):

        def preselection(row):
            if not row.doubleEPass: return False
            if not (row.e1MatchesDoubleEPath > 0 and \
                row.e2MatchesDoubleEPath > 0): return False 
            if not row.e1Pt > 20: return False
            if not selections.eSelection(row, 'e1'): return False
            if not row.e1MVAIDH2TauWP: return False
            if not selections.eSelection(row, 'e2'): return False
            if not selections.vetos(row): return False
            return True

        def fill(the_histos, row):
            the_histos['electronPt'].Fill(row.e2Pt)
            the_histos['electronJetPt'].Fill(max(row.e2JetPt, row.e2Pt))
            the_histos['electronAbsEta'].Fill(row.e2AbsEta)
            #the_histos['metSignificance'].Fill(row.metSignificance)
            the_histos['e1MtToMET'].Fill(row.e1MtToMET)
            the_histos['e2MtToMET'].Fill(row.e2MtToMET)
            the_histos['MET'].Fill(row.type1_pfMetEt)
            the_histos['OSS'].Fill(int(row.e1_e2_SS))

            the_histos['e2JetArea'].Fill(row.e2JetArea)
            the_histos['e2JetEtaEtaMoment'].Fill( row.e2JetEtaEtaMoment )
            the_histos['e2JetEtaPhiMoment'].Fill( row.e2JetEtaPhiMoment )
            the_histos['e2JetEtaPhiSpread'].Fill( row.e2JetEtaPhiSpread )
            the_histos['e2JetPhiPhiMoment'].Fill( row.e2JetPhiPhiMoment )
            the_histos['e2JetptD'].Fill(row.e2JetptD)
            the_histos['e2Jetaxis1'].Fill(row.e2Jetaxis1)
            the_histos['e2Jetaxis2'].Fill(row.e2Jetaxis2)
            the_histos['e2Jetmult'].Fill(row.e2Jetmult)
            the_histos['e2JetmultMLPQC'].Fill(row.e2JetmultMLPQC)
            the_histos['e2JetmultMLP'].Fill(row.e2JetmultMLP)
            the_histos['e2JetQGLikelihoodID'].Fill(row.e2JetQGLikelihoodID)
            the_histos['e2JetQGMVAID'].Fill(row.e2JetQGMVAID)

            the_histos['e2JetQGLikelihoodIDvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2JetQGLikelihoodID)
            the_histos['e2JetQGMVAIDvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2JetQGMVAID)

            the_histos['e2Jetmultvse2JetptD'].Fill(row.e2JetptD,row.e2Jetmult)
            the_histos['e2JetmultMLPvse2JetptD'].Fill(row.e2JetptD,row.e2JetmultMLP) 
            the_histos['e2JetmultMLPQCvse2JetptD'].Fill(row.e2JetptD,row.e2JetmultMLPQC)

            the_histos['e2Jetmultvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2Jetmult)
            the_histos['e2JetmultMLPvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2JetmultMLP)
            the_histos['e2JetmultMLPQCvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2JetmultMLPQC)
            the_histos['e2JetptDvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2JetptD) 


            if row.e1_e2_Mass > 60 and row.e1_e2_Mass < 120 :
                the_histos['MET_Z'].Fill(row.type1_pfMetEt)
                the_histos['OSS_Z'].Fill(int(row.e1_e2_SS))
                the_histos['e1MtToMET_Z'].Fill(row.e1MtToMET)
                the_histos['e2MtToMET_Z'].Fill(row.e2MtToMET)
            else:
                the_histos['MET_NoZ'].Fill(row.type1_pfMetEt)
                the_histos['e1MtToMET_NoZ'].Fill(row.e1MtToMET)
                the_histos['e2MtToMET_NoZ'].Fill(row.e2MtToMET)
                
            the_histos['e1e2Mass'].Fill(row.e1_e2_Mass)
            the_histos['doubleEPrescale'].Fill(row.doubleEPrescale)

        histos = self.histograms
        #pprint( histos)
        for row in self.tree:
            if not preselection(row):
                continue

            region = selections.control_region_ee(row)
            if region is None:
                continue

            if region == 'zee':
                continue
            # This is a QCD or Wjets

            def make_region_plots(full_region):
                fill(histos[full_region], row)
                if row.e2MVAIDH2TauWP:
                    fill(histos[full_region + ( 'mvaid',)], row)
                if row.e2RelPFIsoDB < 0.3:
                    fill(histos[full_region + ( 'iso03',)], row)
                if row.e2MVAIDH2TauWP and row.e2RelPFIsoDB < 0.3:
                    fill(histos[full_region + ( 'mvaidiso03',)], row)
                if row.e2MVAIDH2TauWP and row.e2RelPFIsoDB < 0.1:
                    fill(histos[full_region + ( 'mvaidiso01',)], row)
                if row.e2MVAIDH2TauWP and row.e2RelPFIsoDB < 0.2:
                    fill(histos[full_region + ( 'mvaidiso02',)], row)
                if row.e2MVAIDH2TauWP and ((row.e2RelPFIsoDB < 0.15 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.1):
                    fill(histos[full_region + ( 'h2taucuts',)], row)
                if row.e2MVAIDH2TauWP and ((row.e2RelPFIsoDB < 0.2 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.15):
                    fill(histos[full_region + ( 'h2taucuts020',)], row)
                if row.e2MVAIDH2TauWP and ((row.e2RelPFIsoDB < 0.25 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.20):
                    fill(histos[full_region + ( 'h2taucuts025',)], row)
                    

            make_region_plots((region, 'pt10'))
            if region == 'wjets' and not (row.e1_e2_Mass > 70 and row.e1_e2_Mass < 110):
                make_region_plots(('wjetsNoZmass', 'pt10'))

            if row.e2Pt > 20:
                make_region_plots((region, 'pt20'))
                if region == 'wjets' and not (row.e1_e2_Mass > 70 and row.e1_e2_Mass < 110):
                    make_region_plots(('wjetsNoZmass', 'pt20'))
            

    def finish(self):
        self.write_histos()

