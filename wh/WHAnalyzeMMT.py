'''

Analyze MMT events for the WH analysis

'''

from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import fnmatch
import glob
from MuMuTauTree import MuMuTauTree
import os
import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import FinalStateAnalysis.MetaData.data_views as data_views
import logging
data_views.log.setLevel(logging.INFO)
#import sys
#logging.basicConfig(stream=sys.stderr, level=logging.INFO)
import WHAnalyzerBase
import ROOT
from TwoDimFakeRate import TwoDimFakeRate

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
#### MC-DATA and PU corrections ################################################
################################################################################

# Determine MC-DATA corrections
is7TeV = bool('7TeV' in os.environ['jobid'])
print "Is 7TeV:", is7TeV

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
pu_distributions = glob.glob(os.path.join(
    'inputs', os.environ['jobid'], 'data_DoubleMu*pu.root'))

mc_pu_tag = 'S6' if is7TeV else 'S10'
# Hack to use S6 weights for the HWW 7TeV sample we use in 8TeV
if 'HWW3l' in os.environ.get('megatarget', 'NOTSET') and not is7TeV:
    print "Using S6_600bins PU weights for HWW3l"
    mc_pu_tag = 'S6_600bins'

pu_corrector = PileupWeight.PileupWeight(mc_pu_tag, *pu_distributions)

muon_pog_PFTight_2011 = MuonPOGCorrections.make_muon_pog_PFTight_2011()
muon_pog_PFRelIsoDB02_2011 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB02_2011()

muon_pog_Mu17Mu8_Mu17_2012 = MuonPOGCorrections.make_muon_pog_Mu17Mu8_Mu17_2012()
muon_pog_Mu17Mu8_Mu8_2012 = MuonPOGCorrections.make_muon_pog_Mu17Mu8_Mu8_2012()

# takes etas of muons
muon_pog_Mu17Mu8_2011 = MuonPOGCorrections.muon_pog_Mu17Mu8_eta_eta_2011

# Get object ID and trigger corrector functions
def mc_corrector_2011(row):
    if row.run > 2:
        return 1
    pu = pu_corrector(row.nTruePU)
    #pu = 1
    m1id = muon_pog_PFTight_2011(row.m1Pt, row.m1Eta)
    m2id = muon_pog_PFTight_2011(row.m2Pt, row.m2Eta)
    m1iso = muon_pog_PFRelIsoDB02_2011(row.m1Pt, row.m1Eta)
    m2iso = muon_pog_PFRelIsoDB02_2011(row.m2Pt, row.m2Eta)
    trigger = muon_pog_Mu17Mu8_2011(row.m1Eta, row.m2Eta)
    return pu*m1id*m2id*m1iso*m2iso*trigger

def mc_corrector_2012(row):
    if row.run > 2:
        return 1
    pu = pu_corrector(row.nTruePU)
    m1idiso = H2TauCorrections.correct_mu_idiso_2012(row.m1Pt, row.m1AbsEta)
    m2idiso = H2TauCorrections.correct_mu_idiso_2012(row.m2Pt, row.m2AbsEta)
    m1Trig = muon_pog_Mu17Mu8_Mu17_2012(row.m1Pt, row.m1Eta)
    m2Trig = muon_pog_Mu17Mu8_Mu8_2012(row.m2Pt, row.m2Eta)
    return pu*m1idiso*m2idiso*m1Trig*m2Trig

# Determine which set of corrections to use
mc_corrector = mc_corrector_2011
if not is7TeV:
    mc_corrector = mc_corrector_2012

################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeMMT(WHAnalyzerBase.WHAnalyzerBase):
    tree = 'mmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeMMT, self).__init__(tree, outfile, MuMuTauTree, **kwargs)

    def book_histos(self, folder):
        self.book(folder, "weight", "Event weight", 100, 0, 5)
        self.book(folder, "weight_nopu", "Event weight without PU", 100, 0, 5)
        self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
        self.book(folder, "prescale", "HLT prescale", 26, -5.5, 20.5)
        self.book(folder, "m1Pt", "Muon 1 Pt", 100, 0, 100)
        self.book(folder, "m1JetPt", "Muon 1 Jet Pt", 100, 0, 200)
        self.book(folder, "m2Pt", "Muon 2 Pt", 100, 0, 100)
        self.book(folder, "m2JetBtag", "Muon 2 Pt", 100, -10, 3.3)
        self.book(folder, "m2JetPt", "Muon 2 Jet Pt", 100, 0, 200)
        self.book(folder, "m1AbsEta", "Muon 1 AbsEta", 100, 0, 2.4)
        self.book(folder, "m2AbsEta", "Muon 2 AbsEta", 100, 0, 2.4)
        self.book(folder, "m1m2Mass", "Muon 1-2 Mass", 120, 0, 120)
        self.book(folder, "subMass", "subleadingMass", 200, 0, 200)
        self.book(folder, "leadMass", "leadingMass", 200, 0, 200)
        # Rank muons by less MT to MET, for WZ control region
        self.book(folder, "subMTMass", "subMTMass", 200, 0, 200)
        self.book(folder, "m2Iso", "m2Iso", 100, 0, 0.3)
        self.book(folder, "tPt", "Tau Pt", 100, 0, 100)
        self.book(folder, "tAbsEta", "Tau AbsEta", 100, 0, 2.3)
        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)

    def fill_histos(self, histos, folder, row, weight):
        histos['/'.join(folder + ('nTruePU',))].Fill(row.nTruePU)
        def fill(name, value):
            histos['/'.join(folder + (name,))].Fill(value, weight)
        histos['/'.join(folder + ('weight',))].Fill(weight)

        fill('prescale', row.doubleMuPrescale)
        fill('rho', row.rho)
        fill('nvtx', row.nvtx)
        fill('m1Pt', row.m1Pt)
        fill('m2Pt', row.m2Pt)
        fill('m1JetPt', row.m1JetPt)
        fill('m2JetPt', row.m2JetPt)
        fill('m2JetBtag', row.m2JetBtag)
        fill('m1AbsEta', row.m1AbsEta)
        fill('m2AbsEta', row.m2AbsEta)
        fill('m1m2Mass', row.m1_m2_Mass)
        fill('subMass', row.m2_t_Mass)
        fill('leadMass', row.m1_t_Mass)
        fill('m2Iso', row.m2RelPFIsoDB)
        fill('tPt', row.tPt)
        fill('tAbsEta', row.tAbsEta)
        if row.m1MtToMET > row.m2MtToMET:
            fill('subMTMass', row.m2_t_Mass)
        else:
            fill('subMTMass', row.m1_t_Mass)

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.doubleMuPass:
            return False
        if row.m1Pt < row.m2Pt:
            return False
        if row.m1Pt < 20:
            return False
        if row.m2Pt < 10:
            return False
        if row.tPt < 20:
            return False
        if row.m1AbsEta > 2.4:
            return False
        if row.m2AbsEta > 2.4:
            return False
        if row.tAbsEta > 2.3:
            return False

        if abs(row.m1DZ) > 0.2:
            return False
        if abs(row.m2DZ) > 0.2:
            return False
        if abs(row.tDZ) > 0.2:
            return False

        if row.m1_m2_Mass < 20:
            return False
        if row.LT < 80:
            return False

        if row.muVetoPt5:
            return False
        if row.bjetCSVVeto:
            return False
        if row.tauVetoPt20:
            return False
        if row.eVetoCicTightIso:
            return False

        if not row.m1PixHits:
            return False
        if not row.m2PixHits:
            return False

        ## Fixme use CSV
        if row.m1JetBtag > 3.3:
            return False
        if row.m2JetBtag > 3.3:
            return False

        if not row.tAntiElectronLoose:
            return False
        if row.tCiCTightElecOverlap:
            return False

        if not self.trigger_match_m1(row):
            return False
        if not self.trigger_match_m2(row):
            return False

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
        return mc_corrector(row)

    def obj1_weight(self, row):
        #return highpt_mu_fr(max(row.m1JetPt, row.m1Pt))
        return mu_fr_ewk_2d(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        if row.m1AbsEta < 0.8:
            return mu_fr_ewk_2d_b(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        elif row.m1AbsEta < 1.3:
            return mu_fr_ewk_2d_t(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        else:
            return mu_fr_ewk_2d_f(max(row.m1JetPt, row.m1Pt), row.m1Pt)

    def obj2_weight(self, row):
        #return lowpt_mu_fr(max(row.m2JetPt, row.m2Pt))
        return mu_fr_ewk_2d(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        if row.m2AbsEta < 0.8:
            return mu_fr_ewk_2d_b(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        elif row.m2AbsEta < 1.3:
            return mu_fr_ewk_2d_t(max(row.m2JetPt, row.m2Pt), row.m2Pt)
        else:
            return mu_fr_ewk_2d_f(max(row.m2JetPt, row.m2Pt), row.m2Pt)

    def obj3_weight(self, row):
        return tau_fr(row.tPt)

    def obj1_qcd_weight(self, row):
        #return highpt_mu_qcd_fr(max(row.m1JetPt, row.m1Pt))
        return mu_fr_qcd_2d(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        if row.m1AbsEta < 0.8:
            return mu_fr_qcd_2d_b(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        elif row.m1AbsEta < 1.3:
            return mu_fr_qcd_2d_t(max(row.m1JetPt, row.m1Pt), row.m1Pt)
        else:
            return mu_fr_qcd_2d_f(max(row.m1JetPt, row.m1Pt), row.m1Pt)

    def obj2_qcd_weight(self, row):
        #return lowpt_mu_qcd_fr(max(row.m2JetPt, row.m2Pt))
        return mu_fr_qcd_2d(max(row.m2JetPt, row.m2Pt), row.m2Pt)
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
