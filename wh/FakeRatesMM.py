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

    def begin(self):
        for region in ['wjets', 'qcd', 'all']:
            for denom in ['pt10', 'pt20', 'pt10b', 'pt10t', 'pt10f']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos

                for numerator in ['pfid', 'iso03', 'pfidiso03',
                                  'pfidiso02', 'pfidiso01', 'h2taucuts',
                                  'h2taucuts020', 'h2taucuts025']:
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
                    book_histo('m2JetArea'        , "", 100, 0, 1)
                    book_histo('m2JetEtaEtaMoment', "", 100, 0, 0.2)
                    book_histo('m2JetEtaPhiMoment', "", 100, 0, 0.1)
                    book_histo('m2JetEtaPhiSpread', "", 100, 0, 0.2)
                    book_histo('m2JetPhiPhiMoment', "", 100, 0, 0.1)
                    #book_histo('metSignificance', 'MET sig.', 100, 0, 10)
                    book_histo('m1MtToMET', 'Muon 1 MT', 100, 0, 200)
                    book_histo('m1m2Mass', 'DiMuon Mass', 100, 0, 200)
                    book_histo('m2JetptD', "", 200, 0, 1)
                    book_histo('m2Jetaxis1', "", 200, 0, 1)
                    book_histo('m2Jetaxis2', "", 200, 0, 1)
                    book_histo('m2Jetmult', "", 50, 0, 50)
                    book_histo('m2JetmultMLPQC', "", 50, 0, 50)
                    book_histo('m2JetmultMLP', "", 50, 0, 50)
                    book_histo('m2JetQGLikelihoodID', "", 200, 0, 1)
                    book_histo('m2JetQGMVAID', "", 200, 0, 1)

                    book_histo('m2Jetmultvsm2JetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)
                    book_histo('m2JetmultMLPvsm2JetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)
                    book_histo('m2JetmultMLPQCvsm2JetptD', "",200, 0, 1,50, 0, 50, type=ROOT.TH2F)
                    book_histo('m2JetQGLikelihoodIDvsm2JetPt'," ", 100, 0,100, 200, 0, 1, type=ROOT.TH2F)
                    book_histo('m2JetQGMVAIDvsm2JetPt'," ", 100, 0,100, 200, 0, 1, type=ROOT.TH2F)


                    book_histo('m2Jetmultvsm2JetPt', '', 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('m2JetmultMLPvsm2JetPt', '',  200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('m2JetmultMLPQCvsm2JetPt', '' , 200, 0, 200, 50, 0, 50,type=ROOT.TH2F)
                    book_histo('m2JetptDvsm2JetPt', '' ,  200, 0, 200, 200, 0, 1,type=ROOT.TH2F)





    def process(self):

        def preselection(row):
            if not row.m1_m2_SS: return False
            if not row.doubleMuPass: return False
            if row.m2Pt > row.m1Pt: return False
            if not row.m1Pt > 20: return False
            if not row.m1PFIDTight: return False
            if not row.m2Pt > 10: return False
            if not row.m1AbsEta < 2.4: return False
            if not row.m2AbsEta < 2.4: return False
            if not row.m2JetBtag < 3.3: return False
            if not row.m2PixHits: return False
            if row.eVetoCicTightIso: return False
            if row.muVetoPt5: return False
            if row.bjetCSVVeto: return False
            if row.tauVetoPt20: return False
            if not abs(row.m1DZ) < 0.2: return False
            if not abs(row.m2DZ) < 0.2: return False
            return True

        def fill(the_histos, row):
            # Get PU weight - fix me
            weight = 1
            the_histos['muonPt'].Fill(row.m2Pt, weight)
            the_histos['muonJetPt'].Fill(max(row.m2JetPt, row.m2Pt), weight)
            the_histos['muonJetVsLeptonPt'].Fill(max(row.m2JetPt, row.m2Pt), row.m2Pt, weight)
            the_histos['muonJetVsEta'].Fill(max(row.m2JetPt, row.m2Pt), row.m2AbsEta, weight)
            the_histos['muonAbsEta'].Fill(row.m2AbsEta, weight)
            #the_histos['metSignificance'].Fill(row.metSignificance, weight)
            the_histos['m1MtToMET'].Fill(row.m1MtToMET, weight)
            the_histos['m1m2Mass'].Fill(row.m1_m2_Mass)
            the_histos['m2JetArea'].Fill(row.m2JetArea)
            the_histos['m2JetEtaEtaMoment'].Fill( row.m2JetEtaEtaMoment )
            the_histos['m2JetEtaPhiMoment'].Fill( row.m2JetEtaPhiMoment )
            the_histos['m2JetEtaPhiSpread'].Fill( row.m2JetEtaPhiSpread )
            the_histos['m2JetPhiPhiMoment'].Fill( row.m2JetPhiPhiMoment )
            the_histos['m2JetptD'].Fill(row.m2JetptD)
            the_histos['m2Jetaxis1'].Fill(row.m2Jetaxis1)
            the_histos['m2Jetaxis2'].Fill(row.m2Jetaxis2)
            the_histos['m2Jetmult'].Fill(row.m2Jetmult)
            the_histos['m2JetmultMLPQC'].Fill(row.m2JetmultMLPQC)
            the_histos['m2JetmultMLP'].Fill(row.m2JetmultMLP)
            the_histos['m2JetQGLikelihoodID'].Fill(row.m2JetQGLikelihoodID)
            the_histos['m2JetQGMVAID'].Fill(row.m2JetQGMVAID)
            the_histos['m2Jetmultvsm2JetptD'].Fill(row.m2JetptD,row.m2Jetmult)
            the_histos['m2JetmultMLPvsm2JetptD'].Fill(row.m2JetptD,row.m2JetmultMLP) 
            the_histos['m2JetmultMLPQCvsm2JetptD'].Fill(row.m2JetptD,row.m2JetmultMLPQC)

            the_histos['m2Jetmultvsm2JetPt'].Fill(max(row.m2JetPt, row.m2Pt),row.m2Jetmult)
            the_histos['m2JetmultMLPvsm2JetPt'].Fill(max(row.m2JetPt, row.m2Pt),row.m2JetmultMLP)
            the_histos['m2JetmultMLPQCvsm2JetPt'].Fill(max(row.m2JetPt, row.m2Pt),row.m2JetmultMLPQC)
            the_histos['m2JetptDvsm2JetPt'].Fill(max(row.m2JetPt, row.m2Pt),row.m2JetptD) 

            the_histos['m2JetQGLikelihoodIDvsm2JetPt'].Fill(max(row.m2JetPt, row.m2Pt),row.m2JetQGLikelihoodID)
            the_histos['m2JetQGMVAIDvsm2JetPt'].Fill(max(row.m2JetPt, row.m2Pt),row.m2JetQGMVAID)


        histos = self.histograms
        for row in self.tree:
            if not preselection(row):
                continue
            region = control_region(row)
            if region is None:
                continue

            def fill_region(region, tag):
                # This is a QCD or Wjets
                fill(histos[(region, tag)], row)

                if row.m2RelPFIsoDB < 0.3:
                    fill(histos[(region, tag, 'iso03')], row)

                if row.m2PFIDTight:
                    fill(histos[(region, tag, 'pfid')], row)

                    if row.m2RelPFIsoDB < 0.3:
                        fill(histos[(region, tag, 'pfidiso03')], row)

                    if row.m2RelPFIsoDB < 0.2:
                        fill(histos[(region, tag, 'pfidiso02')], row)

                    if row.m2RelPFIsoDB < 0.1:
                        fill(histos[(region, tag, 'pfidiso01')], row)

                    if (row.m2RelPFIsoDB < 0.15 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.1:
                        fill(histos[(region, tag, 'h2taucuts')], row)

                    if (row.m2RelPFIsoDB < 0.20 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.15:
                        fill(histos[(region, tag, 'h2taucuts020')], row)
                    
                    if (row.m2RelPFIsoDB < 0.25 and row.m2AbsEta < 1.479) or row.m2RelPFIsoDB < 0.20:
                        fill(histos[(region, tag, 'h2taucuts025')], row)
                        
            fill_region(region, 'pt10')
            fill_region('all', 'pt10')

            # Barrel/forward
            if row.m2AbsEta < 0.8:
                fill_region(region,'pt10b')
            elif row.m2AbsEta < 1.3:
                fill_region(region, 'pt10t')
            else:
                fill_region(region, 'pt10f')

            if row.m2Pt > 20:
                fill_region(region, 'pt20')
                fill_region('all', 'pt20')

    def finish(self):
        self.write_histos()
