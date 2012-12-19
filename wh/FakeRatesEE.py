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
#import ROOT.TMath as math
import os
import ROOT
import array

math = ROOT.TMath
def make_corrector_from_th2(filename, path):
    import rootpy.io as io
    import rootpy.plotting
    tfile = io.open(filename)
    hist  = tfile.Get(path).Clone()
    #print hist
    binsx = hist.GetNbinsX()
    binsy = hist.GetNbinsY()
    def refFun(xval,yval):
        #print hist
        xbin = hist.GetXaxis().FindBin(xval)
        xbin = (xbin if xbin <= binsx else binsx ) if xbin >= 1 else 1 #Compute underflow and overflow as first and last bin
        ybin = hist.GetYaxis().FindBin(yval)
        ybin = (ybin if ybin <= binsy else binsy ) if ybin >= 1 else 1 #Compute underflow and overflow as first and last bin
        prob = hist.GetBinContent(xbin,ybin)
        try:
            return prob / (1 - prob)
        except ZeroDivisionError:
            raise ZeroDivisionError(" catched trying to return weight for (%.3f,%.3f) ==> (%i,%i) bin out of (%i,%i). Prob: %.3f. Hist: %s : %s. " % (xval, yval, xbin, ybin, binsx, binsy , prob, filename, path))
    return refFun

frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')

identity         = lambda x,y: 1
charge_flip      = make_corrector_from_th2(frfit_dir+"/chare_flip_prob_map.root", "efficiency_map")         
charge_flip_up   = make_corrector_from_th2(frfit_dir+"/chare_flip_prob_map.root", "efficiency_map_statUp")  
charge_flip_down = make_corrector_from_th2(frfit_dir+"/chare_flip_prob_map.root", "efficiency_map_statDown")

eta_to_theta     = lambda eta: 2*math.ATan( math.Exp(- eta) ) if eta <> 0 else math.Pi() / 2.

def delta_phi(phi1,phi2):
    result = phi1 - phi2;
    while result > math.Pi():
        result -= float(2*math.Pi())
    while result <= -math.Pi():
        result += float(2*math.Pi())
    return result    

def get_lorentz_vector(phi, eta, energy):
    theta = eta_to_theta(eta)
    v1 = ROOT.TLorentzVector(energy*math.Cos(phi)*math.Cos(theta), energy*math.Sin(phi)*math.Cos(theta), energy*math.Sin(theta), energy)
    return v1

def cos_theta_electrons(row):
    v1 = get_versor(row.e1SCPhi, row.e1SCEta)
    v2 = get_versor(row.e2SCPhi, row.e2SCEta)
    return v1.Dot(v2)

def sc_inv_mass(row):
    pt1 = row.e1SCEnergy / math.CosH(row.e1SCEta)
    pt2 = row.e2SCEnergy / math.CosH(row.e2SCEta)
    return math.Sqrt(2*pt1*pt2*( math.CosH(row.e1SCEta - row.e2SCEta) - math.Cos(row.e1SCPhi - row.e2SCPhi) ) )

def control_region(row):
    # Figure out what control region we are in.
    if row.e1_e2_SS and row.e1RelPFIsoDB < 0.15 and row.e1MtToMET > 50 and row.metSignificance > 3:
        return 'wjets'
    elif row.e1_e2_SS and row.e1RelPFIsoDB > 0.3 and row.metSignificance < 3:
        return 'qcd'
    elif row.e1RelPFIsoDB < 0.1 and row.e2RelPFIsoDB < 0.1 and row.e2MVAIDH2TauWP and row.e1MVAIDH2TauWP \
        and not any([ row.muVetoPt5,
                      row.tauVetoPt20,
                      row.eVetoCicTightIso,
                      ]):
        return 'zee'
    else:
        return None
 
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
        for region in ['wjets', 'qcd']:
            for denom in ['pt10', 'pt20']:
                denom_key = (region, denom)
                denom_histos = {}
                self.histograms[denom_key] = denom_histos

                for numerator in ['mvaid', 'iso03', 'mvaidiso03',
                                  'mvaidiso01','h2taucuts']:
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

                    book_histo('e2Pt', 'e Pt', 100, 0, 100)
                    book_histo('e2JetPt', 'e Jet Pt', 100, 0, 100)
                    book_histo('e2AbsEta', 'e Abs Eta', 100, -2.5, 2.5)
                    book_histo('metSignificance', 'MET sig.', 100, 0, 10)
                    book_histo('e1MtToMET', 'm MT', 100, 0, 200)
                    book_histo('e1e2Mass', 'DiElectron Mass', 100, 0, 200)
                    book_histo('doubleEPrescale', 'prescale', 10, -0.5, 9.5)
        # Charge mis-ID measurements
        # Charge mis-ID measurements
        for app in ["","_weight","_weightSysUp","_weightSysDwn"]:
            self.book('charge', "OS"+app+'_ePt', 'e Pt', 100, 0, 100)
            self.book('charge', "OS"+app+'_eAbsEta', 'e Abs Eta', 80, 0, 2.5)
            self.book('charge', "OS"+app+'_SCEnergy', 'electron Super Cluster energy', 500, 0, 1000)
            self.book('charge', "OS"+app+'_SCDPhi', 'electron Super Cluster DeltaPhi', 180, 0, math.Pi())
            self.book('charge', "OS"+app+"_TrkMass", 'OS%s Dielectrons invariant mass; M_{ee} [GeV];counts' % app, 110, 40, 150)
            self.book('charge', "OS"+app+"_SCMass" , 'OS%s Dielectrons Super Cluster invariant mass; M_{ee} [GeV];counts' % app, 110, 40, 150)
        self.book('charge', "SS_SCMass", 'SS Dielectrons Super Cluster invariant mass; M_{ee} [GeV];counts', 110, 40, 150)
        self.book('charge', "SS_TrkMass", 'SS Dielectrons invariant mass; M_{ee} [GeV];counts', 110, 40, 150)
        self.book('charge', "SS"+'_SCEnergy', 'electron Super Cluster energy', 500, 0, 1000)
        self.book('charge', "SS"+'_SCDPhi', 'electron Super Cluster Delta phi', 180, 0, math.Pi())
        self.book('charge', 'SS_ePt', 'e Pt', 100, 0, 100)
        self.book('charge', 'SS_eAbsEta', 'e Abs Eta', 80, 0, 2.5)
            
        #print self.histograms.keys()

    def process(self):

        def preselection(row):
            if not row.doubleEPass: return False
            if not row.e1Pt > 20: return False
            if not selections.eSelection(row, 'e1'): return False
            if not row.e1MVAIDH2TauWP: return False
            if not selections.eSelection(row, 'e2'): return False
            if not selections.vetos(row): return False
            return True

        def fill(the_histos, row):
            the_histos['e2Pt'].Fill(row.e2Pt)
            the_histos['e2JetPt'].Fill(row.e2JetPt)
            the_histos['e2AbsEta'].Fill(row.e2AbsEta)
            the_histos['metSignificance'].Fill(row.metSignificance)
            the_histos['e1MtToMET'].Fill(row.e1MtToMET)
            the_histos['e1e2Mass'].Fill(row.e1_e2_Mass)
            the_histos['doubleEPrescale'].Fill(row.doubleEPrescale)

        histos = self.histograms
        for row in self.tree:
            if not preselection(row):
                continue

            region = control_region(row)
            if region is None:
                continue

            if region == 'zee':
                scmass = sc_inv_mass(row)
                if row.e1_e2_SS:
                    histos['charge/SS_eAbsEta'].Fill(row.e2AbsEta)
                    histos['charge/SS_eAbsEta'].Fill(row.e1AbsEta)
                    histos['charge/SS_ePt'].Fill(row.e2Pt)
                    histos['charge/SS_ePt'].Fill(row.e1Pt)
                    histos['charge/SS_SCEnergy'].Fill(row.e1SCEnergy)
                    histos['charge/SS_SCEnergy'].Fill(row.e2SCEnergy)
                    histos['charge/SS_SCDPhi'].Fill(delta_phi(row.e1SCPhi, row.e2SCPhi) )
                    histos['charge/SS_TrkMass'].Fill(row.e1_e2_Mass)
                    histos['charge/SS_SCMass'].Fill(scmass)
                for app, fcn in zip(["","_weight","_weightSysUp","_weightSysDwn"],[identity, charge_flip, charge_flip_up, charge_flip_down]):
                    weight = ( fcn(row.e1AbsEta,row.e1Pt) + fcn(row.e2AbsEta,row.e2Pt) )
                    #weight = prob / (1 - prob)
                    histos['charge/OS'+app+'_ePt'].Fill(row.e2Pt,         weight ) ## fcn(row.e2AbsEta,row.e2Pt))
                    histos['charge/OS'+app+'_ePt'].Fill(row.e1Pt,         weight ) ## fcn(row.e1AbsEta,row.e1Pt))
                    histos['charge/OS'+app+'_eAbsEta'].Fill(row.e2AbsEta, weight ) ## fcn(row.e2AbsEta,row.e2Pt))
                    histos['charge/OS'+app+'_eAbsEta'].Fill(row.e1AbsEta, weight ) ## fcn(row.e1AbsEta,row.e1Pt))
                    histos['charge/OS'+app+'_SCEnergy'].Fill(row.e1SCEnergy, weight )
                    histos['charge/OS'+app+'_SCEnergy'].Fill(row.e2SCEnergy, weight )
                    histos['charge/OS'+app+'_SCDPhi'].Fill(delta_phi(row.e1SCPhi, row.e2SCPhi), weight )
                    histos['charge/OS'+app+'_TrkMass'].Fill(row.e1_e2_Mass, weight)
                    histos['charge/OS'+app+'_SCMass'].Fill(scmass, weight)
                continue
            # This is a QCD or Wjets
            fill(histos[(region, 'pt10')], row)

            if row.e2MVAIDH2TauWP:
                fill(histos[(region, 'pt10', 'mvaid')], row)

            if row.e2RelPFIsoDB < 0.3:
                fill(histos[(region, 'pt10', 'iso03')], row)

            if row.e2MVAIDH2TauWP and row.e2RelPFIsoDB < 0.3:
                fill(histos[(region, 'pt10', 'mvaidiso03')], row)

            if row.e2MVAIDH2TauWP and row.e2RelPFIsoDB < 0.1:
                fill(histos[(region, 'pt10', 'mvaidiso01')], row)

            if row.e2MVAIDH2TauWP and ((row.e2RelPFIsoDB < 0.15 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.1):
                fill(histos[(region, 'pt10', 'h2taucuts')], row)

            if row.e2Pt > 20:
                fill(histos[(region, 'pt20')], row)
                if row.e2MVAIDH2TauWP:
                    fill(histos[(region, 'pt20', 'mvaid')], row)

                if row.e2RelPFIsoDB < 0.3:
                    fill(histos[(region, 'pt20', 'iso03')], row)

                if row.e2MVAIDH2TauWP and row.e2RelPFIsoDB < 0.3:
                    fill(histos[(region, 'pt20', 'mvaidiso03')], row)

                if row.e2MVAIDH2TauWP and row.e2RelPFIsoDB < 0.1:
                    fill(histos[(region, 'pt20', 'mvaidiso01')], row)

                if row.e2MVAIDH2TauWP and ((row.e2RelPFIsoDB < 0.15 and row.e2AbsEta < 1.479) or row.e2RelPFIsoDB < 0.1):
                    fill(histos[(region, 'pt10', 'h2taucuts')], row)
                

    def finish(self):
        self.write_histos()

