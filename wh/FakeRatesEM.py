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

import ROOT
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
                                  'eid13Looseidiso02', 'eid13Looseh2taucuts', 'eid13Looseh2taucuts020',
                                  'eid13Tightidiso02', 'eid13Tighth2taucuts', 'eid13Tighth2taucuts020',
                                  ]:
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

                    book_histo('ePt', 'e Pt', 100, 0, 100)
                    book_histo('eJetPt', 'e Jet Pt', 100, 0, 100)
                    book_histo('eAbsEta', 'e Abs Eta', 100, -2.5, 2.5)
                    #book_histo('metSignificance', 'MET sig.', 100, 0, 10)
                    book_histo('mMtToMET', 'm MT', 100, 0, 200)
                    book_histo('eJetArea'        , "", 100, 0, 1)
                    book_histo('eJetEtaEtaMoment', "", 100, 0, 0.2)
                    book_histo('eJetEtaPhiMoment', "", 100, 0, 0.1)
                    book_histo('eJetEtaPhiSpread', "", 100, 0, 0.2)
                    book_histo('eJetPhiPhiMoment', "", 100, 0, 0.1)
                    book_histo('eJetptD', "", 200, 0, 1)
                    book_histo('eJetaxis1', "", 200, 0, 1)
                    book_histo('eJetaxis2', "", 200, 0, 1)
                    book_histo('eJetmult', "", 50, 0, 50)
                    book_histo('eJetmultMLPQC', "", 50, 0, 50)
                    book_histo('eJetmultMLP', "", 50, 0, 50)
                    book_histo('eJetQGLikelihoodID', "", 200, 0, 1)
                    book_histo('eJetQGMVAID', "", 200, 0, 1)
                    book_histo('eJetQGLikelihoodIDvseJetPt'," ", 100, 0,100, 200, 0, 1, type=ROOT.TH2F)
                    book_histo('eJetQGMVAIDvseJetPt'," ", 100, 0,100, 200, 0, 1, type=ROOT.TH2F)
                    book_histo('eJetmultvseJetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)
                    book_histo('eJetmultMLPvseJetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)
                    book_histo('eJetmultMLPQCvseJetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)

                    book_histo('eJetmultvseJetPt', '', 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('eJetmultMLPvseJetPt', '',  200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('eJetmultMLPQCvseJetPt', '' , 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('eJetptDvseJetPt', '' ,  200, 0, 200, 200, 0, 1,type=ROOT.TH2F)

    def process(self):

        def preselection(row):
            if not row.e_m_SS: return False
            if not selections.muSelection(row, 'm'):  return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
            if not selections.eSelection(row, 'e'):   return False #applies basic selection (eta, pt > 10, DZ, missingHits, jetBTag, HasConversion and chargedIdTight)
            if not row.eChargeIdTight:                return False
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

            the_histos['eJetArea'].Fill(row.mJetArea)
            the_histos['eJetEtaEtaMoment'].Fill( row.mJetEtaEtaMoment )
            the_histos['eJetEtaPhiMoment'].Fill( row.mJetEtaPhiMoment )
            the_histos['eJetEtaPhiSpread'].Fill( row.mJetEtaPhiSpread )
            the_histos['eJetPhiPhiMoment'].Fill( row.mJetPhiPhiMoment )
            the_histos['eJetptD'].Fill(row.mJetptD)
            the_histos['eJetaxis1'].Fill(row.mJetaxis1)
            the_histos['eJetaxis2'].Fill(row.mJetaxis2)
            the_histos['eJetmult'].Fill(row.mJetmult)
            the_histos['eJetmultMLPQC'].Fill(row.mJetmultMLPQC)
            the_histos['eJetmultMLP'].Fill(row.mJetmultMLP)
            the_histos['eJetQGLikelihoodID'].Fill(row.mJetQGLikelihoodID)
            the_histos['eJetQGMVAID'].Fill(row.mJetQGMVAID)
            the_histos['eJetQGLikelihoodIDvseJetPt'].Fill(max(row.eJetPt, row.ePt),row.eJetQGLikelihoodID)
            the_histos['eJetQGMVAIDvseJetPt'].Fill(max(row.eJetPt, row.ePt),row.eJetQGMVAID)
            the_histos['eJetmultvseJetptD'].Fill(row.eJetptD,row.eJetmult)
            the_histos['eJetmultMLPvseJetptD'].Fill(row.eJetptD,row.eJetmultMLP) 
            the_histos['eJetmultMLPQCvseJetptD'].Fill(row.eJetptD,row.eJetmultMLPQC)

            the_histos['eJetmultvseJetPt'].Fill(max(row.eJetPt, row.ePt),row.eJetmult)
            the_histos['eJetmultMLPvseJetPt'].Fill(max(row.eJetPt, row.ePt),row.eJetmultMLP)
            the_histos['eJetmultMLPQCvseJetPt'].Fill(max(row.eJetPt, row.ePt),row.eJetmultMLPQC)
            the_histos['eJetptDvseJetPt'].Fill(max(row.eJetPt, row.ePt),row.eJetptD) 

 

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
                    fill(histos[(region, pt_cut, 'eid13Looseidiso02')], row)

                if (row.eRelPFIsoDB < 0.15 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.1:
                    fill(histos[(region, pt_cut, 'eid13Looseh2taucuts')], row)

                if (row.eRelPFIsoDB < 0.20 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.15:
                    fill(histos[(region, pt_cut, 'eid13Looseh2taucuts020')], row)
            
            if selections.summer_2013_eid_tight(row, 'e'):
                if row.eRelPFIsoDB < 0.2:
                    fill(histos[(region, pt_cut, 'eid13Tightidiso02')], row)

                if (row.eRelPFIsoDB < 0.15 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.1:
                    fill(histos[(region, pt_cut, 'eid13Tighth2taucuts')], row)

                if (row.eRelPFIsoDB < 0.20 and row.eAbsEta < 1.479) or row.eRelPFIsoDB < 0.15:
                    fill(histos[(region, pt_cut, 'eid13Tighth2taucuts020')], row)
                
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
