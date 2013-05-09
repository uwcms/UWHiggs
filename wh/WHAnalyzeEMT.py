'''

Analyze EMT events for the WH analysis

'''

from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
from EMuTauTree import EMuTauTree
import glob
import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import mcCorrectors
import optimizer as selections
import fakerate_functions as frfits
import ROOT
import math

################################################################################
#### Analysis logic ############################################################
################################################################################
is7TeV = bool('7TeV' in os.environ['jobid'])
use_iso_trigger = not is7TeV

class WHAnalyzeEMT(WHAnalyzerBase):
    tree = 'emt/final/Ntuple'

    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeEMT, self).__init__(tree, outfile, EMuTauTree, **kwargs)
        #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.hfunc['subMass']   = lambda row, weight: (row.e_t_Mass, weight)    if row.ePt < row.mPt else (row.m_t_Mass, weight) 
        self.hfunc['tLeadDR']   = lambda row, weight: (row.m_t_DR,   weight)    if row.ePt < row.mPt else (row.e_t_DR,   weight) 
        self.hfunc['tSubDR']    = lambda row, weight: (row.e_t_DR,   weight)    if row.ePt < row.mPt else (row.m_t_DR,   weight) 
        self.hfunc['pt_ratio' ] = lambda row, weight: (row.ePt/row.mPt, weight) if row.ePt < row.mPt else (row.mPt/row.ePt, weight)
        self.hfunc["e*_t_Mass"] = lambda row, weight: ( frfits.mass_scaler( row.e_t_Mass), weight)
        self.hfunc["e*_m_Mass"] = lambda row, weight: ( frfits.mass_scaler( row.e_m_Mass), weight)
        self.hfunc["subMass*" ] = lambda row, weight: ( frfits.mass_scaler( row.e_t_Mass), weight)    if row.ePt < row.mPt else (row.m_t_Mass, weight)
        self.hfunc["_recoilDaught" ] = lambda row, weight: (math.sqrt(row.recoilDaught) , weight)
        self.hfunc["_recoilWithMet"] = lambda row, weight: (math.sqrt(row.recoilWithMet), weight)

        self.pucorrector = mcCorrectors.make_puCorrector('mueg')

    def book_histos(self, folder):
        self.book(folder, "mPt", "Muon Pt", 100, 0, 100)
        self.book(folder, "ePt", "Electron Pt", 100, 0, 100)
        self.book(folder, "tPt", "Tau Pt", 100, 0, 100)
        self.book(folder, "mAbsEta", "Muon AbsEta", 100, 0, 2.4)
        self.book(folder, "eAbsEta", "Electron AbsEta", 100, 0, 2.5)
        self.book(folder, "tAbsEta", "Tau AbsEta", 100, 0, 2.3)

        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
        self.book(folder, "e_m_Mass", "Electron-Muon Mass", 200, 0, 200)
        self.book(folder, "eChargeIdTight", "Elec charge ID tight", 2, -0.5, 1.5)
        #self.book(folder, "eChargeIdMedium", "Elec charge ID medium", 2, -0.5, 1.5)
        self.book(folder, "e_t_Mass", "Electron-Tau Mass", 200, 0, 200)
        self.book(folder, "subMass", "Subleading Mass", 200, 0, 200)
        #self.book(folder, "bCSVVeto", "BjetCSV", 10, -0.5, 9.5)
        #self.book(folder, "metSig", "MET significance", 100, 0, 15)
        self.book(folder, "tLeadDR", "DR between leading lepton and tau", 100, 0, 5)
        self.book(folder, "tSubDR", "DR between subleading lepton and tau", 100, 0, 5)
        self.book(folder, "LT", "L_T", 100, 0, 300)
        #Charge mis-id special histograms
        if 'c1' in folder:
            self.book(folder, "e*_t_Mass", "Electron-Tau Mass", 200, 0, 200)
            self.book(folder, "subMass*", "Subleading Mass", 200, 0, 200)
            self.book(folder, "e*_m_Mass", "Electron-Muon Mass", 200, 0, 200)
        
        #let's look for osme other possible selections
        self.book(folder, "Mass"          , "mass"          , 800, 0, 800 )
        self.book(folder, "pt_ratio"      , "pt_ratio"      , 100, 0, 1)
        self.book(folder, "tToMETDPhi"    , "tToMETDPhi"    , 100, 0, 4)
        self.book(folder, "type1_pfMetEt"         , "metEt"         , 300, 0, 2000)


    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    def preselection( self, row, cut_flow_trk = None):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not ( abs(row.mGenMotherPdgId) in [15, 24, 23, 21] ): return False
        cut_flow_trk.Fill('obj1 GenMatching') 
        if not ( abs(row.eGenMotherPdgId) in [15, 24, 23, 21] ): return False
        #if -1*row.eCharge*row.eGenPdgId != 11: return False
        cut_flow_trk.Fill('obj2 GenMatching') 
        if row.tGenDecayMode < 0: return False
        cut_flow_trk.Fill('obj3 GenMatching') 

        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        mu8e17 = (row.mu8ele17isoPass and row.ePt >= 20) #if use_iso_trigger else (row.mu17ele8Pass and row.mPt < 20)
        if not (mu17e8 or mu8e17):                return False
        cut_flow_trk.Fill('trigger')

        if not selections.muSelection(row, 'm'):  return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        cut_flow_trk.Fill('obj1 Presel')

        if not selections.eSelection(row, 'e'):   return False #applies basic selection (eta, pt > 10, DZ, missingHits, jetBTag, HasConversion and chargedIdTight)
        if row.eJetBtag > 3.3:                    return False
        cut_flow_trk.Fill('obj2 Presel')

        if not selections.tauSelection(row, 't'): return False #applies basic selection (eta, pt > 20, DZ)
        if row.tMuOverlap:         return False
        if not row.tAntiMuonTight: return False
        cut_flow_trk.Fill('obj3 Presel')

        if row.LT < 80:            return False
        cut_flow_trk.Fill('LT')

        if row.e_m_SS and row.e_t_SS            : return False #remove three SS leptons
        if not selections.vetos(row):             return False #applies mu bjet e additional tau vetoes

        #FIXME: remove after tests
        if not selections.leading_lepton_id_iso(row, 'm'): return False
        cut_flow_trk.Fill('obj1 IDIso')
        if not selections.subleading_lepton_id_iso(row, 'e'): return False
        cut_flow_trk.Fill('obj2 IDIso')
        if not row.tLooseIso3Hits: return False
        cut_flow_trk.Fill('obj3 IDIso')

        if row.muVetoPt5IsoIdVtx: return False
        cut_flow_trk.Fill('mu veto')
        if row.eVetoMVAIsoVtx:    return False
        cut_flow_trk.Fill('e veto')
        if row.tauVetoPt20Loose3HitsVtx: return False
        cut_flow_trk.Fill('tau veto')
        if row.bjetCSVVeto:       return False
        cut_flow_trk.Fill('bjet veto')
        cut_flow_trk.Fill('charge_fakes') #no charge fakes here

        return True

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def sign_cut( row):
        ''' Returns true if muons are SS '''
        return bool(row.e_m_SS)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj1_id( row):
        return selections.leading_lepton_id_iso(row, 'm')

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj2_id( row):
        return selections.subleading_lepton_id_iso(row, 'e')

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj3_id( row):
        return bool(row.tLooseIso3Hits)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def anti_wz( row):
        if row.e_t_Zcompat < 20:
            if not row.tAntiElectronMVA3Medium:
                return False
        elif not row.tAntiElectronMVA3Loose:
            return False
        return True

    def enhance_wz(self, row):
        if row.e_t_Zcompat < 15 and not row.tAntiElectronMVA:
            return True
        return False

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_muon_corrections(row,'m') * \
            mcCorrectors.get_electron_corrections(row,'e') * \
            mcCorrectors.correct_mueg_mu(row.mPt, row.mAbsEta) * \
            mcCorrectors.correct_mueg_e(row.ePt, row.eAbsEta)

    def obj1_weight(self, row):
        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        if mu17e8:
            return frfits.highpt_mu_fr(max(row.mJetPt, row.mPt))
        else:
            return frfits.lowpt_mu_fr(max(row.mJetPt, row.mPt))

    def obj2_weight(self, row):
        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        if mu17e8:
            return frfits.lowpt_e_fr(max(row.eJetPt, row.ePt))
        else:
            return frfits.highpt_e_fr(max(row.eJetPt, row.ePt))

    def obj3_weight(self, row):
        return frfits.tau_fr(row.tPt)

    def obj1_qcd_weight(self, row):
        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        if mu17e8:
            return frfits.highpt_mu_qcd_fr(max(row.mJetPt, row.mPt))
        else:
            return frfits.lowpt_mu_qcd_fr(max(row.mJetPt, row.mPt))

    def obj2_qcd_weight(self, row):
        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        if mu17e8:
            return frfits.lowpt_e_qcd_fr(max(row.eJetPt, row.ePt))
        else:
            return frfits.highpt_e_qcd_fr(max(row.eJetPt, row.ePt))

    def obj3_qcd_weight(self, row):
        return frfits.tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return row.m_t_SS

    ## def obj1_charge_flip(self, row):
    ##     return 0

    def obj2_charge_flip(self, row):
        return frfits.e_charge_flip(row.eAbsEta,row.ePt)

    def obj2_charge_flip_sysup(self, row):
        return frfits.e_charge_flip_up(row.eAbsEta,row.ePt) #lowpt_e_charge_flip
