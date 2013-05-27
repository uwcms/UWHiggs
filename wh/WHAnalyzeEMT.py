'''

Analyze EMT events for the WH analysis

'''

from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
from EMuTauTree import EMuTauTree
import glob
import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits
import ROOT
import math
import optimizer

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
        self.hfunc["subPt"]     = lambda row, weight: (row.e_t_Pt, weight)      if row.ePt < row.mPt else (row.m_t_Pt, weight) 
        self.hfunc['pt_ratio' ] = lambda row, weight: (row.ePt/row.mPt, weight) if row.ePt < row.mPt else (row.mPt/row.ePt, weight)
        self.hfunc["e*_t_Mass"] = lambda row, weight: ( frfits.mass_scaler['h2taucuts']( row.e_t_Mass), weight)
        self.hfunc["e*_m_Mass"] = lambda row, weight: ( frfits.mass_scaler['h2taucuts']( row.e_m_Mass), weight)
        self.hfunc["subMass*" ] = lambda row, weight: ( frfits.mass_scaler['h2taucuts']( row.e_t_Mass), weight) if row.ePt < row.mPt else (row.m_t_Mass, weight)
        self.hfunc["_recoilDaught" ] = lambda row, weight: (math.sqrt(row.recoilDaught) , weight)
        self.hfunc["_recoilWithMet"] = lambda row, weight: (math.sqrt(row.recoilWithMet), weight)

        self.pucorrector = mcCorrectors.make_puCorrector('mueg')

    def book_histos(self, folder):
        for key in optimizer.grid_search:
            prefix = key+'$' if key else ''
            self.book(folder, prefix+"subMass", "Subleading Mass", 200, 0, 200)
            self.book(folder, prefix+"LT", "L_T", 100, 0, 300)
            self.book(folder, prefix+"subPt"   , "SubLeading Pt", 400, 0, 400)
            #Charge mis-id special histograms
            if 'c1' in folder:
                self.book(folder, prefix+"subMass*", "Subleading Mass", 200, 0, 200)

        if len(optimizer.grid_search.keys()) == 1:
            if 'c1' in folder:
                self.book(folder, "e*_t_Mass", "Electron-Tau Mass", 200, 0, 200)
                self.book(folder, "e*_m_Mass", "Electron-Muon Mass", 200, 0, 200)
            self.book(folder, "e_m_Mass", "Electron-Muon Mass", 200, 0, 200)
            self.book(folder, "e_t_Mass", "Electron-Tau Mass", 200, 0, 200)
            self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
            self.book(folder, "mPt", "Muon Pt", 100, 0, 100)
            self.book(folder, "ePt", "Electron Pt", 100, 0, 100)
            self.book(folder, "tPt", "Tau Pt", 100, 0, 100)
            self.book(folder, "mAbsEta", "Muon AbsEta", 100, 0, 2.4)
            self.book(folder, "eAbsEta", "Electron AbsEta", 100, 0, 2.5)
            self.book(folder, "tAbsEta", "Tau AbsEta", 100, 0, 2.3)
            self.book(folder, "eChargeIdTight", "Elec charge ID tight", 2, -0.5, 1.5)
            self.book(folder, "tLeadDR", "DR between leading lepton and tau", 100, 0, 5)
            self.book(folder, "tSubDR", "DR between subleading lepton and tau", 100, 0, 5)
        
            #let's look for osme other possible selections
            self.book(folder, "Mass"          , "mass"          , 800, 0, 800 )
            self.book(folder, "pt_ratio"      , "pt_ratio"      , 100, 0, 1)
            self.book(folder, "tToMETDPhi"    , "tToMETDPhi"    , 100, 0, 4)
            self.book(folder, "type1_pfMetEt"         , "metEt"         , 300, 0, 2000)


    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    def preselection( self, row, cut_flow_trk = None, LT_threshold = 80.):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        mu17e8, mu8e17 = False, False
        if use_iso_trigger:
            mu17e8 = row.mu17ele8isoPass and \
                row.mMatchesMu17Ele8IsoPath > 0 and \
                row.eMatchesMu17Ele8IsoPath > 0 and \
                row.mPt >= 20
            mu8e17 = row.mu8ele17isoPass and \
                row.mMatchesMu8Ele17IsoPath > 0 and \
                row.eMatchesMu8Ele17IsoPath > 0 and \
                row.ePt >= 20
        else:
            mu17e8 = row.mu17ele8Pass and \
                row.mMatchesMu17Ele8Path > 0 and \
                row.eMatchesMu17Ele8Path > 0 and \
                row.mPt >= 20
            mu8e17 = row.mu8ele17Pass and \
                row.mMatchesMu8Ele17Path > 0 and \
                row.eMatchesMu8Ele17Path > 0 and \
                row.ePt >= 20
            
        if not (mu17e8 or mu8e17):                return False
        cut_flow_trk.Fill('trigger')

        if not selections.muSelection(row, 'm'):  return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        cut_flow_trk.Fill('obj1 Presel')

        if not selections.eSelection(row, 'e'):   return False #applies basic selection (eta, pt > 10, DZ, missingHits, jetBTag, HasConversion and chargedIdTight)
        cut_flow_trk.Fill('obj2 Presel')

        if not selections.tauSelection(row, 't'): return False #applies basic selection (eta, pt > 20, DZ)
        if row.tMuOverlap:         return False
        if not row.tAntiMuonTight: return False
        cut_flow_trk.Fill('obj3 Presel')

        if row.LT < LT_threshold:            return False
        cut_flow_trk.Fill('LT')

        if row.e_m_SS and row.e_t_SS: return False #remove three SS leptons
        if not selections.vetos(row): return False #applies mu bjet e additional tau vetoes
        cut_flow_trk.Fill('vetos')
        cut_flow_trk.Fill('charge_fakes') #no charge fakes here

        return True

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def sign_cut( row):
        ''' Returns true if muons are SS '''
        return bool(row.e_m_SS)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj1_id( row, ledleptonId='h2taucuts', subledleptonId='h2taucuts'):
        return selections.lepton_id_iso(row, 'm', ledleptonId) \
            if row.mPt > row.ePt else \
            selections.lepton_id_iso(row, 'm', subledleptonId) 

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj2_id( row, ledleptonId='h2taucuts', subledleptonId='h2taucuts'):
        return selections.lepton_id_iso(row, 'e', ledleptonId) \
            if row.mPt > row.ePt else \
            selections.lepton_id_iso(row, 'e', subledleptonId)

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

    def obj1_weight(self, row, ledleptonId='h2taucuts', subledleptonId='h2taucuts'):
        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        if mu17e8:
            return frfits.highpt_mu_fr[ledleptonId](max(row.mJetPt, row.mPt)) \
                if row.mPt > row.ePt else \
                frfits.highpt_mu_fr[subledleptonId](max(row.mJetPt, row.mPt))
        else:
            return frfits.lowpt_mu_fr[ledleptonId](max(row.mJetPt, row.mPt)) \
                if row.mPt > row.ePt else \
                frfits.lowpt_mu_fr[subledleptonId](max(row.mJetPt, row.mPt)) 

    def obj2_weight(self, row, ledleptonId='h2taucuts', subledleptonId='h2taucuts'):
        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        if mu17e8:
            return frfits.lowpt_e_fr[ledleptonId](max(row.eJetPt, row.ePt)) \
                if row.ePt > row.mPt else \
                frfits.lowpt_e_fr[subledleptonId](max(row.eJetPt, row.ePt)) 
        else:
            return frfits.highpt_e_fr[ledleptonId](max(row.eJetPt, row.ePt)) \
                if row.ePt > row.mPt else \
                frfits.highpt_e_fr[subledleptonId](max(row.eJetPt, row.ePt))

    def obj3_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_fr(row.tPt)

    def obj1_qcd_weight(self, row, ledleptonId='h2taucuts', subledleptonId='h2taucuts'):
        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        if mu17e8:
            return frfits.highpt_mu_qcd_fr[ledleptonId](max(row.mJetPt, row.mPt)) \
                if row.mPt > row.ePt else \
                frfits.highpt_mu_qcd_fr[subledleptonId](max(row.mJetPt, row.mPt))
        else:
            return frfits.lowpt_mu_qcd_fr[ledleptonId](max(row.mJetPt, row.mPt)) \
                if row.mPt > row.ePt else \
                frfits.lowpt_mu_qcd_fr[subledleptonId](max(row.mJetPt, row.mPt)) 

    def obj2_qcd_weight(self, row, ledleptonId='h2taucuts', subledleptonId='h2taucuts'):
        mu17e8 = (row.mu17ele8isoPass and row.mPt >= 20) if use_iso_trigger else (row.mu17ele8Pass and row.mPt >= 20)
        if mu17e8:
            return frfits.lowpt_e_qcd_fr[ledleptonId](max(row.eJetPt, row.ePt)) \
                if row.ePt > row.mPt else \
                frfits.lowpt_e_qcd_fr[subledleptonId](max(row.eJetPt, row.ePt)) 
        else:
            return frfits.highpt_e_qcd_fr[ledleptonId](max(row.eJetPt, row.ePt)) \
                if row.ePt > row.mPt else \
                frfits.highpt_e_qcd_fr[subledleptonId](max(row.eJetPt, row.ePt))

    def obj3_qcd_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return row.m_t_SS

    ## def obj1_charge_flip(self, row):
    ##     return 0

    def obj2_charge_flip(self, row, ledleptonId='h2taucuts', subledleptonId='h2taucuts'):
        return frfits.e_charge_flip[ledleptonId](row.eAbsEta,row.ePt) \
            if row.ePt > row.mPt else \
            frfits.e_charge_flip[subledleptonId](row.eAbsEta,row.ePt)

    def obj2_charge_flip_sysup(self, row, ledleptonId='h2taucuts', subledleptonId='h2taucuts'):
        return frfits.e_charge_flip_up[ledleptonId](row.eAbsEta,row.ePt) \
            if row.ePt > row.mPt else \
            frfits.e_charge_flip_up[subledleptonId](row.eAbsEta,row.ePt) \
