'''

Make control plots of Z and ttbar -> emu events.

Author: Evan K. Friis, UW

'''

import EMuTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor

###############################################################################
#### MC-DATA and PU corrections ###############################################
###############################################################################

# Determine MC-DATA corrections
is7TeV = bool('7TeV' in os.environ['jobid'])
print "Is 7TeV:", is7TeV
use_iso_trigger = not is7TeV

mc_pu_tag = 'S6' if is7TeV else 'S10'
if 'HWW3l' in os.environ['megatarget'] and not is7TeV:
    print "Using S6_600bins PU weights for HWW3l"
    mc_pu_tag = 'S6_600bins'
    use_iso_trigger = False

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
pu_distributions = glob.glob(os.path.join(
    'inputs', os.environ['jobid'], 'data_MuEG*pu.root'))
pu_corrector = PileupWeight.PileupWeight(mc_pu_tag, *pu_distributions)

muon_pog_PFTight_2011 = MuonPOGCorrections.make_muon_pog_PFTight_2011()
muon_pog_PFRelIsoDB_2011 = MuonPOGCorrections.\
    make_muon_pog_PFRelIsoDB012_2011()

frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')
lowpt_e_fr = build_roofunctor(
    #frfit_dir + '/e_qcd_pt10_h2taucuts_ePt.root',
    frfit_dir + '/e_wjets_pt10_h2taucuts_ePt.root',
    'fit_efficiency',  # workspace name
    'efficiency'
)


def mc_corrector_2011(row):
    # Get object ID and trigger corrector functions
    if row.run > 2:
        return 1
    pu = pu_corrector(row.nTruePU)
    m1idiso = muon_pog_PFRelIsoDB_2011(row.mPt, row.mAbsEta) * \
        muon_pog_PFTight_2011(row.mPt, row.mAbsEta)
    e2idiso = H2TauCorrections.correct_e_idiso_2011(row.ePt, row.eAbsEta)
    m_trg = H2TauCorrections.correct_mueg_mu_2011(row.mPt, row.mAbsEta)
    e_trg = H2TauCorrections.correct_mueg_e_2011(row.ePt, row.eAbsEta)
    return pu * m1idiso * e2idiso * m_trg * e_trg


def mc_corrector_2012(row):
    if row.run > 2:
        return 1
    pu = pu_corrector(row.nTruePU)
    m1idiso = H2TauCorrections.correct_mu_idiso_2012(row.mPt, row.mAbsEta)
    e2idiso = H2TauCorrections.correct_e_idiso_2012(row.ePt, row.eAbsEta)
    m_trg = H2TauCorrections.correct_mueg_mu_2012(row.mPt, row.mAbsEta)
    e_trg = H2TauCorrections.correct_mueg_e_2012(row.ePt, row.eAbsEta)
    return pu * m1idiso * e2idiso * m_trg * e_trg

# Determine which set of corrections to use
mc_corrector = mc_corrector_2011
if not is7TeV:
    mc_corrector = mc_corrector_2012


class ControlEM(MegaBase):
    tree = 'em/final/Ntuple'

    def __init__(self, tree, outfile, **kwargs):
        super(ControlEM, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EMuTree.EMuTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):
        for folder in ['', '/ss', '/f2', '/f2/w2']:
            self.book('em' + folder, "weight", "Event weight", 100, 0, 5)
            self.book('em' + folder, "weight_nopu", "Event weight without PU",
                      100, 0, 5)

            self.book('em' + folder, "rho", "Fastjet #rho", 100, 0, 25)
            self.book('em' + folder, "nvtx", "Number of vertices",
                      31, -0.5, 30.5)
            self.book('em' + folder, "prescale", "HLT prescale",
                      21, -0.5, 20.5)
            self.book('em' + folder, "group", "HLT group", 21, -0.5, 20.5)

            self.book('em' + folder, "mPt", "Muon 1 Pt", 100, 0, 100)
            self.book('em' + folder, "ePt", "Muon 2 Pt", 100, 0, 100)
            self.book('em' + folder, "mAbsEta", "Muon 1 eta", 100, -2.5, 2.5)
            self.book('em' + folder, "eAbsEta", "Muon 2 eta", 100, -2.5, 2.5)
            self.book('em' + folder, "emMass", "m_{e#mu} (GeV)", 240, 0, 240)

            self.book('em' + folder, 'mPixHits', 'Mu 1 pix hits',
                      10, -0.5, 9.5)

            self.book('em' + folder, 'mJetBtag', 'Mu 1 JetBtag',
                      100, -5.5, 9.5)
            self.book('em' + folder, 'eJetBtag', 'Mu 2 JetBtag',
                      100, -5.5, 9.5)

            # Vetoes
            self.book('em' + folder, 'bjetVeto', 'Number of b-jets',
                      5, -0.5, 4.5)
            self.book('em' + folder, 'bjetCSVVeto', 'Number of b-jets',
                      5, -0.5, 4.5)
            self.book('em' + folder, 'muVetoPt5', 'Number of extra muons',
                      5, -0.5, 4.5)
            self.book('em' + folder, 'tauVetoPt20', 'Number of extra taus',
                      5, -0.5, 4.5)
            self.book('em' + folder, 'eVetoCicTightIso',
                      'Number of extra CiC tight electrons', 5, -0.5, 4.5)

    def correction(self, row):
        return mc_corrector(row)

    def fill_histos(self, row):
        histos = self.histograms
        weight = self.correction(row)

        def fill_folder(x, w):
            histos[x + '/weight'].Fill(w)
            histos[x + '/rho'].Fill(row.rho, w)
            histos[x + '/nvtx'].Fill(row.nvtx, w)
            histos[x + '/prescale'].Fill(row.mu17ele8Prescale, w)
            histos[x + '/group'].Fill(row.mu17ele8Group, w)
            histos[x + '/ePt'].Fill(row.ePt, w)
            histos[x + '/mPt'].Fill(row.mPt, w)
            histos[x + '/eAbsEta'].Fill(row.eAbsEta, w)
            histos[x + '/mAbsEta'].Fill(row.mAbsEta, w)
            histos[x + '/mPixHits'].Fill(row.mPixHits, w)
            histos[x + '/eJetBtag'].Fill(row.eJetBtag, w)
            histos[x + '/mJetBtag'].Fill(row.mJetBtag, w)
            histos[x + '/emMass'].Fill(row.e_m_Mass, w)

            histos[x + '/bjetVeto'].Fill(row.bjetVeto, w)
            histos[x + '/bjetCSVVeto'].Fill(row.bjetCSVVeto, w)
            histos[x + '/muVetoPt5'].Fill(row.muVetoPt5, w)
            histos[x + '/tauVetoPt20'].Fill(row.tauVetoPt20, w)
            histos[x + '/eVetoCicTightIso'].Fill(row.eVetoCicTightIso, w)

        passes_e_id_iso = self.obj2_id(row)
        if row.e_m_SS and passes_e_id_iso:
            fill_folder('em/ss', weight)
        elif passes_e_id_iso:
            fill_folder('em', weight)
        else:
            fill_folder('em/f2', weight)
            fill_folder('em/f2/w2', weight * self.obj2_weight(row))

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not use_iso_trigger:
            if not row.mu17ele8Pass:
                return False
        elif not row.mu17ele8isoPass:
            return False

        if row.mPt < 20:
            return False
        if row.ePt < 10:
            return False
        if row.mAbsEta > 2.4:
            return False
        if row.eAbsEta > 2.5:
            return False
        if abs(row.eDZ) > 0.2:
            return False
        if abs(row.mDZ) > 0.2:
            return False
        if row.muVetoPt5:
            return False
        if row.bjetVeto:
            return False
        if row.tauVetoPt20:
            return False
        if row.eVetoCicTightIso:
            return False
        if not row.mPixHits:
            return False
        if row.eMissingHits:
            return False
        if row.eHasConversion:
            return False
        if not row.eChargeIdTight:
            return False
        return True

    def obj1_id(self, row):
        return bool(row.mPFIDTight) and (
            row.mRelPFIsoDB < 0.1 or
            (row.mRelPFIsoDB < 0.15 and row.mAbsEta < 1.479))

    def obj2_id(self, row):
        return bool(row.eMVAIDH2TauWP) and bool(
            row.eRelPFIsoDB < 0.1 or
            (row.eRelPFIsoDB < 0.15 and row.eAbsEta < 1.479))

    def obj2_weight(self, row):
        fr = lowpt_e_fr(row.ePt)
        return fr / (1. - fr)

    def process(self):
        for row in self.tree:
            if not self.preselection(row):
                continue
            if not self.obj1_id(row):
                continue
            self.fill_histos(row)

    def finish(self):
        self.write_histos()
