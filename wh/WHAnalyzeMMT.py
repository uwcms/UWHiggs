'''

Analyze MMT events for the WH analysis

'''

from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import fnmatch
import glob
from MuMuTauTree import MuMuTauTree
import os
import FinalStateAnalysis.MetaData.data_views as data_views
import logging
data_views.log.setLevel(logging.INFO)
#import sys
#logging.basicConfig(stream=sys.stderr, level=logging.INFO)
import WHAnalyzerBase
import ROOT
from TwoDimFakeRate import TwoDimFakeRate
import mcCorrectors
import baseSelections as selections

################################################################################
#### Fitted fake rate functions ################################################
################################################################################

# Get fitted fake rate functions
frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')
highpt_mu_fr = build_roofunctor(
    #frfit_dir + '/m_wjets_pt20_pfidiso02_muonJetPt.root',
    frfit_dir + '/m_wjets_pt20_h2taucuts_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
lowpt_mu_fr = build_roofunctor(
    #frfit_dir + '/m_wjets_pt10_pfidiso02_muonJetPt.root',
    frfit_dir + '/m_wjets_pt10_h2taucuts_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
tau_fr = build_roofunctor(
    frfit_dir + '/t_ztt_pt20_mvaloose_tauPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

highpt_mu_qcd_fr = build_roofunctor(
    frfit_dir + '/m_qcd_pt20_h2taucuts_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
lowpt_mu_qcd_fr = build_roofunctor(
    frfit_dir + '/m_qcd_pt10_h2taucuts_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
tau_qcd_fr = build_roofunctor(
    frfit_dir + '/t_ztt_pt20_mvaloose_tauPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

# Get 2D fake rates

fr_data_views = data_views.data_views(
    glob.glob(os.path.join('results', os.environ['jobid'], 'FakeRatesMM', '*.root')),
    glob.glob(os.path.join('inputs', os.environ['jobid'], '*.sum')),
)

def get_view(sample_pattern):
    for sample, sample_info in fr_data_views.iteritems():
        if fnmatch.fnmatch(sample, sample_pattern):
            return sample_info['view']
    raise KeyError("I can't find a view that matches %s, I have: %s" % (
        sample_pattern, " ".join(fr_data_views.keys())))

# FR data, subtracting WZ and ZZ.
mu_fr_ewk_2d = TwoDimFakeRate(
    'wjets/pt10/h2taucuts/muonJetVsLeptonPt', 'wjets/pt10/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

mu_fr_qcd_2d = TwoDimFakeRate(
    'qcd/pt10/h2taucuts/muonJetVsLeptonPt', 'qcd/pt10/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

# eta dependent jet-pt vs pt
mu_fr_ewk_2d_f = TwoDimFakeRate(
    'wjets/pt10f/h2taucuts/muonJetVsLeptonPt', 'wjets/pt10f/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))
mu_fr_qcd_2d_f = TwoDimFakeRate(
    'qcd/pt10f/h2taucuts/muonJetVsLeptonPt', 'qcd/pt10f/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

mu_fr_ewk_2d_t = TwoDimFakeRate(
    'wjets/pt10t/h2taucuts/muonJetVsLeptonPt', 'wjets/pt10t/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))
mu_fr_qcd_2d_t = TwoDimFakeRate(
    'qcd/pt10t/h2taucuts/muonJetVsLeptonPt', 'qcd/pt10t/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

mu_fr_ewk_2d_b = TwoDimFakeRate(
    'wjets/pt10b/h2taucuts/muonJetVsLeptonPt', 'wjets/pt10b/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))
mu_fr_qcd_2d_b = TwoDimFakeRate(
    'qcd/pt10b/h2taucuts/muonJetVsLeptonPt', 'qcd/pt10b/muonJetVsLeptonPt',
    get_view('data'), get_view('WZ*'), get_view('ZZ*'))

if __name__ == "__main__":
    mu_fr_ewk_2d.plot("ewk_2d_frs.png")
    mu_fr_qcd_2d.plot("qcd_2d_frs.png")
    mu_fr_ewk_2d_f.plot("ewk_2df_frs.png")
    mu_fr_qcd_2d_f.plot("qcd_2df_frs.png")
    mu_fr_ewk_2d_b.plot("ewk_2db_frs.png")
    mu_fr_qcd_2d_b.plot("qcd_2db_frs.png")
    mu_fr_ewk_2d_t.plot("ewk_2dt_frs.png")
    mu_fr_qcd_2d_t.plot("qcd_2dt_frs.png")


################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeMMT(WHAnalyzerBase.WHAnalyzerBase):
    tree = 'mmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeMMT, self).__init__(tree, outfile, MuMuTauTree, **kwargs)
        self.hfunc['subMTMass'] = lambda row, weight: (row.m2_t_Mass, weight) if row.m1MtToMET > row.m2MtToMET else (row.m1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')

    def book_histos(self, folder):
        self.book(folder, "weight", "Event weight", 100, 0, 5)
        #self.book(folder, "weight_nopu", "Event weight without PU", 100, 0, 5) #Booked but not filled??
        self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
        self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)
        self.book(folder, "m1Pt", "Muon 1 Pt", 100, 0, 100)
        self.book(folder, "m1JetPt", "Muon 1 Jet Pt", 100, 0, 200)
        self.book(folder, "m2Pt", "Muon 2 Pt", 100, 0, 100)
        self.book(folder, "m2JetBtag", "Muon 2 Pt", 100, -10, 3.3)
        self.book(folder, "m2JetPt", "Muon 2 Jet Pt", 100, 0, 200)
        self.book(folder, "m1AbsEta", "Muon 1 AbsEta", 100, 0, 2.4)
        self.book(folder, "m2AbsEta", "Muon 2 AbsEta", 100, 0, 2.4)
        self.book(folder, "m1_m2_Mass", "Muon 1-2 Mass", 120, 0, 120)
        self.book(folder, "m2_t_Mass", "subleadingMass", 200, 0, 200)
        self.book(folder, "m1_t_Mass", "leadingMass", 200, 0, 200)
        # Rank muons by less MT to MET, for WZ control region
        self.book(folder, "subMTMass", "subMTMass", 200, 0, 200)
        self.book(folder, "m2RelPFIsoDB", "m2Iso", 100, 0, 0.3)
        self.book(folder, "tPt", "Tau Pt", 100, 0, 100)
        self.book(folder, "tAbsEta", "Tau AbsEta", 100, 0, 2.3)
        self.book(folder, "tDecayMode", "Tau AbsEta", 15, -0.5, 14.5)
        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.doubleMuPass:                  return False
        if row.m1Pt < row.m2Pt:                   return False
        if row.m1Pt < 20:                         return False
        if not selections.muSelection(row, 'm1'): return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        if not selections.muSelection(row, 'm2'): return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        if not selections.tauSelection(row, 't'): return False #applies basic selection (eta, pt > 20, DZ)

        if row.m1_m2_Mass < 20:    return False
        if row.LT < 80:            return False

        if not selections.vetos(row): return False #applies mu bjet e additional tau vetoes

        if not row.tAntiElectronLoose: return False
        if row.tCiCTightElecOverlap:   return False

        if not self.trigger_match_m1(row): return False
        if not self.trigger_match_m2(row): return False

        return True

    @staticmethod
    def trigger_match_m1(row):
        return True
        if row.m1DiMuonL3p5PreFiltered8  > 0 or \
           row.m1DiMuonL3PreFiltered7  > 0 or \
           row.m1SingleMu13L3Filtered13  > 0 or \
           row.m1SingleMu13L3Filtered17  > 0 or \
           row.m1DiMuonMu17Mu8DzFiltered0p2  > 0 or \
           row.m1L3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17:
            return True

    @staticmethod
    def trigger_match_m2(row):
        return True
        if row.m2DiMuonL3p5PreFiltered8  > 0 or \
           row.m2DiMuonL3PreFiltered7  > 0 or \
           row.m2SingleMu13L3Filtered13  > 0 or \
           row.m2SingleMu13L3Filtered17  > 0 or \
           row.m2DiMuonMu17Mu8DzFiltered0p2  > 0 or \
           row.m2L3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17:
            return True

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return bool(row.m1_m2_SS)

    def obj1_id(self, row):
        return bool(row.m1PFIDTight) and (
            row.m1RelPFIsoDB < 0.1 or
            (row.m1RelPFIsoDB < 0.15 and row.m1AbsEta < 1.479))

    def obj2_id(self, row):
        return bool(row.m2PFIDTight) and (
            row.m2RelPFIsoDB < 0.1 or
            (row.m2RelPFIsoDB < 0.15 and row.m2AbsEta < 1.479))

    def obj3_id(self, row):
        return bool(row.tLooseMVAIso)

    def anti_wz(self, row):
        return row.tAntiMuonTight and not row.tMuOverlap

    def enhance_wz(self, row):
        # Require the "tau" to be a muon, and require the third muon
        # to have M_Z +- 20
        if row.tAntiMuonTight or not row.tMuOverlap:
            return False
        # Cut on m2 PT > 20
        #if row.m2Pt < 20:
            #return False
        # Make sure any Z is from m1
        m2_good_Z = bool(71 < row.m2_t_Mass < 111)
        return not m2_good_Z

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m1','m2') * \
            mcCorrectors.double_muon_trigger(row,'m1','m2')

    def obj1_weight(self, row):
        return highpt_mu_fr(max(row.m1JetPt, row.m1Pt))
        #return mu_fr_ewk_2d(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        if row.m1AbsEta < 0.8:
            return mu_fr_ewk_2d_b(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        elif row.m1AbsEta < 1.3:
            return mu_fr_ewk_2d_t(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        else:
            return mu_fr_ewk_2d_f(max(row.m1JetPt, row.m1Pt), row.m1Pt)

    def obj2_weight(self, row):
        return lowpt_mu_fr(max(row.m2JetPt, row.m2Pt))
        #return mu_fr_ewk_2d(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        if row.m2AbsEta < 0.8:
            return mu_fr_ewk_2d_b(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        elif row.m2AbsEta < 1.3:
            return mu_fr_ewk_2d_t(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        else:
            return mu_fr_ewk_2d_f(max(row.m2JetPt, row.m2Pt), row.m2Pt)

    def obj3_weight(self, row):
        return tau_fr(row.tPt)

    def obj1_qcd_weight(self, row):
        return highpt_mu_qcd_fr(max(row.m1JetPt, row.m1Pt))
        #return mu_fr_qcd_2d(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        if row.m1AbsEta < 0.8:
            return mu_fr_qcd_2d_b(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        elif row.m1AbsEta < 1.3:
            return mu_fr_qcd_2d_t(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        else:
            return mu_fr_qcd_2d_f(max(row.m1JetPt, row.m1Pt), row.m1Pt)

    def obj2_qcd_weight(self, row):
        return lowpt_mu_qcd_fr(max(row.m2JetPt, row.m2Pt))
        #return mu_fr_qcd_2d(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        if row.m2AbsEta < 0.8:
            return mu_fr_qcd_2d_b(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        elif row.m2AbsEta < 1.3:
            return mu_fr_qcd_2d_t(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        else:
            return mu_fr_qcd_2d_f(max(row.m2JetPt, row.m2Pt), row.m2Pt)

    def obj3_qcd_weight(self, row):
        return tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return not row.m1_t_SS

    def obj1_charge_flip(self, row):
        return 0
