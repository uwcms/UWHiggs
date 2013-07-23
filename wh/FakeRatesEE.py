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
import mcCorrectors
import os
from array import array
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
        self.pucorrector = mcCorrectors.make_puCorrector('doublee')
        self.defined_eids = selections.electronIds.keys()
        self.iso_points   = ['idiso02', 'h2taucuts', 'h2taucuts020']
        self.lepIds  =  [ '_'.join([i,j])
                          for i in self.defined_eids
                          for j in self.iso_points]

    def begin(self):
        for region in ['wjets', 'qcd', 'wjetsNoZmass', 'qcdNoZmass']:
            for denom in ['pt10', 'pt20']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos

                #SPECIAL NTUPLE!
                denom_histos['electronInfo'] = self.book(
                    os.path.join(region, denom),
                    'electronInfo', "electronInfo", 
                    'electronPt:electronJetPt:weight:'+':'.join(self.lepIds), 
                    type=ROOT.TNtuple)

                for numerator in self.lepIds:
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
                    book_histo('MET', 'MET', 100, 0, 400)
                    book_histo('OSS', 'SS to OS', 2, 0, 2)
                    book_histo('e1e2Mass', 'DiElectron Mass', 100, 0, 200)
                    book_histo('doubleEPrescale', 'prescale', 10, -0.5, 9.5)

                    book_histo('e2JetArea'        , "", 100, 0, 1)
                    book_histo('e2JetptD', "", 200, 0, 1)
                    book_histo('e2Jetmult', "", 50, 0, 50)

                    book_histo('e2Jetmultvse2JetPt', '', 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
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
            if not (row.jetVeto40_DR05 >= 1):             return False
            if not selections.vetos(row): return False
            return True

        def fill(the_histos, row, fillNtuple=False):
            weight = 1
            if row.run == 1:
                weight = self.pucorrector(row.nTruePU)

            the_histos['electronPt'].Fill(row.e2Pt)
            the_histos['electronJetPt'].Fill(max(row.e2JetPt, row.e2Pt))
            the_histos['electronAbsEta'].Fill(row.e2AbsEta)
            #the_histos['metSignificance'].Fill(row.metSignificance)
            the_histos['e1MtToMET'].Fill(row.e1MtToMET)
            the_histos['e2MtToMET'].Fill(row.e2MtToMET)
            the_histos['MET'].Fill(row.type1_pfMetEt)
            the_histos['OSS'].Fill(int(row.e1_e2_SS))

            the_histos['e2JetptD'].Fill(row.e2JetptD)
            the_histos['e2Jetmult'].Fill(row.e2Jetmult)

            the_histos['e2Jetmultvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2Jetmult)
            the_histos['e2JetptDvse2JetPt'].Fill(max(row.e2JetPt, row.e2Pt),row.e2JetptD) 
                
            the_histos['e1e2Mass'].Fill(row.e1_e2_Mass)
            the_histos['doubleEPrescale'].Fill(row.doubleEPrescale)

            if fillNtuple:
                id_iso_vals = [  selections.lepton_id_iso(row, 'e2', label)  for label in self.lepIds]
                #print id_iso_vals
                #print self.lepIds
                id_iso_vals = [float( i ) for i in id_iso_vals ]
                the_histos['electronInfo'].Fill( array("f", [row.e2Pt, row.e2JetPt, weight]+id_iso_vals) )


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
                fill(histos[full_region], row, True)

                for idlabel, idfcn in selections.electronIds.iteritems():                
                    if not idfcn(row, 'e2'): #if it does not pass id skip!
                        continue
                    if row.e2RelPFIsoDB < 0.2:
                        fill(histos[full_region + ( idlabel+'_idiso02',)], row)
                    if (row.e2RelPFIsoDB < 0.15 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.1:
                        fill(histos[full_region + ( idlabel+'_h2taucuts',)], row)
                    if (row.e2RelPFIsoDB < 0.20 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.15:
                        fill(histos[full_region + ( idlabel+'_h2taucuts020',)], row)


            make_region_plots((region, 'pt10'))
            if not (row.e1_e2_Mass > 60 and row.e1_e2_Mass < 120):
                make_region_plots((region+'NoZmass', 'pt10'))

            if row.e2Pt > 20:
                make_region_plots((region, 'pt20'))
                if not (row.e1_e2_Mass > 60 and row.e1_e2_Mass < 120):
                    make_region_plots((region+'NoZmass', 'pt20'))
            

    def finish(self):
        self.write_histos()

