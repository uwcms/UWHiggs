'''

Analyze MMT events for the WH analysis

'''

from MuMuTauTree import MuMuTauTree
import ROOT
import os
from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits
import math
from FinalStateAnalysis.PlotTools.decorators import memo_last
import optimizer
################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeMMT(WHAnalyzerBase):
    tree = 'mmt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        self.channel = 'MMT'
        super(WHAnalyzeMMT, self).__init__(tree, outfile, MuMuTauTree, **kwargs)

	def attr_getter(attribute):
            def f(row, weight):
                return (getattr(row,attribute), weight)
            return f

        def merge_functions(fcn_1, fcn_2):
            def f(row, weight):
                r1, w1 = fcn_1(row, weight)
                r2, w2 = fcn_2(row, weight)
                w = w1 if w1 and w2 else None
                return ((r1, r2), w)
            return f

        lead_iso = self.grid_search['']['leading_iso']
        sublead_iso = self.grid_search['']['subleading_iso']
        
        @memo_last
        def f_par_prob(m1Pt, m1JetPt,
                       m2Pt, m2JetPt,
                       tPt):
            p_m1 = (( frfits.highpt_mu_fr[lead_iso](muonJetPt=max(m1JetPt, m1Pt), muonPt=m1Pt) +\
                      frfits.highpt_mu_qcd_fr[lead_iso](muonJetPt=max(m1JetPt, m1Pt), muonPt=m1Pt) )/2)
            p_m2 = (( frfits.lowpt_mu_fr[sublead_iso](muonJetPt=max(m2JetPt, m2Pt), muonPt=m2Pt) + \
                      frfits.lowpt_mu_qcd_fr[sublead_iso](muonJetPt=max(m2JetPt, m2Pt), muonPt=m2Pt))/2)
            p_t  = frfits.tau_fr(tPt)
            return (p_m1 + p_m2*(1 - p_m1) + p_t*(1 - p_m1)*(1 - p_m2))
            
        def f_prob(row, weight):
            val = f_par_prob(row.m1Pt, row.m1JetPt,
                             row.m2Pt, row.m2JetPt,
                             row.tPt)
            return val, weight

        def log_prob(row, weight):
            prob, weight = f_prob(row, weight)
            return ROOT.TMath.Log10(prob), weight
        
        self.hfunc['faking_prob'] = f_prob
        self.hfunc['log_prob']    = log_prob
        self.hfunc["m2_t_Mass#faking_prob"] = merge_functions( attr_getter('m2_t_Mass'), f_prob  )
        self.hfunc["m2_t_Mass#log_prob"   ] = merge_functions( attr_getter('m2_t_Mass'), log_prob)

        self.hfunc['subMTMass'] = lambda row, weight: (row.m2_t_Mass, weight) if row.m1MtToMET > row.m2MtToMET else (row.m1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.hfunc['pt_ratio' ] = lambda row, weight: (row.m2Pt/row.m1Pt, weight)
        
        self.pucorrector = mcCorrectors.make_puCorrector('doublemu')

    def book_histos(self, folder):
        for key in self.grid_search:
            prefix = key+'$' if key else ''
            self.book(folder, prefix+"m2_t_Mass", "subleadingMass", 200, 0, 200)
            self.book(folder, prefix+"LT" ,  "LT" , 100, 0., 500)
            self.book(folder, prefix+"m2_t_Pt", "subleadingPt", 400, 0, 400)

        if len(self.grid_search.keys()) == 1:

            self.book(folder, "m2_t_Mass#LT" , "subleadingMass", 200, 0, 200, 120, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m2_t_Mass#tPt", "subleadingMass", 200, 0, 200, 200, 0, 200, type=ROOT.TH2F)
            #self.book(folder, "m2_t_Mass#faking_prob" , "subleadingMass", 200, 0, 200, 1100, 0., 1.1, type=ROOT.TH2F)
            #self.book(folder, "m2_t_Mass#log_prob"    , "subleadingMass", 200, 0, 200, 1000, -10,  1, type=ROOT.TH2F)
            #self.book(folder, 'faking_prob'     , "", 1100, 0., 1.1)
            #self.book(folder, 'log_prob'        , "", 1000, -10, 1)

            #Pt
            self.book(folder, "m1Pt#LT" , "subleadingMass", 150, 0, 150, 120, 0, 600, type=ROOT.TH2F)
            self.book(folder, "tPt#LT" , "subleadingMass", 150, 0, 150, 120, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m2Pt#LT" , "subleadingMass", 150, 0, 150, 120, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m1JetPt#LT" , "subleadingMass", 150, 0, 150, 120, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m2JetPt#LT" , "subleadingMass", 150, 0, 150, 120, 0, 600, type=ROOT.TH2F)

            #eta
            self.book(folder, "m1AbsEta#LT" , "subleadingMass", 100, 0, 2.5, 120, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m2AbsEta#LT" , "subleadingMass", 100, 0, 2.5, 120, 0, 600, type=ROOT.TH2F)
            self.book(folder, "tAbsEta#LT" , "subleadingMass", 100, 0, 2.5, 120, 0, 600, type=ROOT.TH2F)

            #DR
            self.book(folder, "m1_t_DR#LT" , "subleadingMass", 100, 0, 10, 120, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m2_t_DR#LT" , "subleadingMass", 100, 0, 10, 120, 0, 600, type=ROOT.TH2F)

            #Jet BTag
            self.book(folder, "m2JetBtag#LT", "Muon 2 Pt", 100, -100, 100, 60, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m1JetBtag#LT", "Muon 2 Pt", 100, -100, 100, 60, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m1JetCSVBtag#LT", "Muon 2 Pt", 120, -5, 1, 60, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m2JetCSVBtag#LT", "Muon 2 Pt", 120, -5, 1, 60, 0, 600, type=ROOT.TH2F)


            #Mt To MET
            self.book(folder, "m1MtToMET#LT", "Muon 2 Pt", 150, 0, 150, 60, 0, 600, type=ROOT.TH2F)
            self.book(folder, "m2MtToMET#LT", "Muon 2 Pt", 150, 0, 150, 60, 0, 600, type=ROOT.TH2F)

            self.book(folder, "m2RelPFIsoDB", "m2Iso", 100, 0, 0.3)
            self.book(folder, "m1_t_Mass", "leadingMass", 200, 0, 200)
            self.book(folder, "m1_m2_Mass", "Muon 1-2 Mass", 120, 0, 120)
            self.book(folder, "m2_t_DR", "m2_t_DR", 100, 0, 5)
            self.book(folder, "m1_t_DR", "m1_t_DR", 100, 0, 5)

            self.book(folder, "m1JetPt", "Muon 1 Jet Pt", 100, 0, 200)
            self.book(folder, "m2JetPt", "Muon 2 Jet Pt", 100, 0, 200)
            # Rank muons by less MT to MET, for WZ control region
            self.book(folder, "weight", "Event weight", 100, 0, 5)
            self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
            self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
            self.book(folder, "m1Pt", "Muon 1 Pt", 100, 0, 100)
            self.book(folder, "m2Pt", "Muon 2 Pt", 100, 0, 100)
            self.book(folder, "m1AbsEta", "Muon 1 AbsEta", 100, 0, 2.4)
            self.book(folder, "m2AbsEta", "Muon 2 AbsEta", 100, 0, 2.4)
            self.book(folder, "tPt", "Tau Pt", 100, 0, 100)
            self.book(folder, "tAbsEta", "Tau AbsEta", 100, 0, 2.3)
            self.book(folder, "subMTMass", "subMTMass", 200, 0, 200)
            #self.book(folder, "tDecayMode", "Tau AbsEta", 15, -0.5, 14.5)
            self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
            self.book(folder, "m1DZ",  "m1DZ", 100, 0., 1)
            self.book(folder, "m2DZ",  "m2DZ", 100, 0., 1)
            self.book(folder, "tDZ" ,  "tDZ" , 100, 0., 1)
            self.book(folder, "doubleMuPrescale", "HLT prescale", 26, -5.5, 20.5)

            #let's look for osme other possible selections
            self.book(folder, "pt_ratio"      , "pt_ratio"      , 100, 0, 1)
            self.book(folder, "tToMETDPhi"    , "tToMETDPhi"    , 100, 0, 4)
            self.book(folder, "Mass"          , "mass"          , 800, 0, 800 )
            self.book(folder, "type1_pfMetEt"         , "metEt"         , 300, 0, 2000)

    def preselection(self, row, cut_flow_trk = None, LT_threshold = 80., taupt_thr = 0.):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        double_mu_pass =  row.doubleMuPass and \
            row.m1MatchesDoubleMuPaths > 0 and \
            row.m2MatchesDoubleMuPaths > 0
        double_muTrk_pass = row.doubleMuTrkPass and \
             row.m1MatchesDoubleMuTrkPaths > 0 and \
             row.m2MatchesDoubleMuTrkPaths > 0
        if not ( double_mu_pass or double_muTrk_pass ): return False
        cut_flow_trk.Fill('trigger')

        if row.m1Pt < row.m2Pt:                   return False
        if row.m1Pt < 20:                         return False
        if not selections.muSelection(row, 'm1'): return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        cut_flow_trk.Fill('obj1 Presel')
        #cut_flow_trk.Fill('pt requirements 1', 'eta requirements 1', 'MissingHits 1', 'HasConversion 1', 'JetBtag 1', 'DZ 1',)


        if not selections.muSelection(row, 'm2'): return False #applies basic selection (eta, pt > 10, DZ, pixHits, jetBTag)
        cut_flow_trk.Fill('obj2 Presel')
        #cut_flow_trk.Fill('pt requirements 2', 'eta requirements 2', 'MissingHits 2', 'HasConversion 2', 'JetBtag 2', 'DZ 2',)

        if not selections.tauSelection(row, 't'): return False #applies basic selection (eta, pt > 20, DZ)
        if not row.tAntiElectronMVA3Loose:        return False
        if row.tPt < taupt_thr: return False
        cut_flow_trk.Fill('obj3 Presel')

        if row.LT < LT_threshold: return False
        cut_flow_trk.Fill('LT')

        if row.m1_m2_SS and row.m1_t_SS: return False #remove three SS leptons
        if row.m1_m2_Mass < 20:          return False
        if not selections.vetos(row):    return False #applies mu bjet e additional tau vetoes
        cut_flow_trk.Fill('vetos')
        #cut_flow_trk.Fill('ChargeIdTight')
        cut_flow_trk.Fill('charge_fakes') #no charge fakes here
        
        #FIXME: ONLY FOR CUT-FLOW PRODUCTION
        #if not row.m1PFIDTight: return False
        #cut_flow_trk.Fill('obj1 ID')
        #if not selections.lepton_id_iso(row, 'm1', 'h2taucuts'): return False
        #cut_flow_trk.Fill('obj1 Iso')
        #if not row.m2PFIDTight: return False
        #cut_flow_trk.Fill('obj2 ID')
        #if not selections.lepton_id_iso(row, 'm2', 'h2taucuts'): return False
        #cut_flow_trk.Fill('obj2 Iso')
        #if not row.tLooseIso3Hits: return False
        #cut_flow_trk.Fill('obj3 IDIso')
        
        return True

    @staticmethod
    def sign_cut(row):
        ''' Returns true if muons are SS '''
        return bool(row.m1_m2_SS)

    @staticmethod
    def obj1_id(row, leadleptonId='h2taucuts', subleadleptonId=None):
        return selections.lepton_id_iso(row, 'm1', leadleptonId)

    @staticmethod
    def obj2_id(row, leadleptonId=None, subleadleptonId='h2taucuts'):
        return selections.lepton_id_iso(row, 'm2', subleadleptonId)

    @staticmethod
    def obj3_id(row, tauId=None, LT_threshold = 80., taupt_thr = 0.):
        if row.LT >= LT_threshold and row.tPt >= taupt_thr:
            return bool(row.tLooseIso3Hits)
        else:
            return bool( getattr(row, tauId) )
    
    @staticmethod
    def anti_wz(row):
        return row.tAntiMuonTight # and not row.tMuOverlap

    def enhance_wz(self, row):
        # Require the "tau" to be a muon, and require the third muon
        # to have M_Z +- 20
        if self.anti_wz(row):
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

    def obj1_weight(self, row, ledleptonId='h2taucuts', subledleptonId=None):
        return frfits.highpt_mu_fr[ledleptonId](muonJetPt=max(row.m1JetPt, row.m1Pt), muonPt=row.m1Pt, muonJetCSVBtag=max(0, row.m1JetCSVBtag)) #muonPVDXY=row.m1PVDXY) #, muonJetBtag=row.m1JetBtag)

    def obj2_weight(self, row, ledleptonId=None, subledleptonId='h2taucuts'):
        return frfits.lowpt_mu_fr[subledleptonId](muonJetPt=max(row.m2JetPt, row.m2Pt), muonPt=row.m2Pt, muonJetCSVBtag=max(0, row.m2JetCSVBtag)) #, muonPVDXY=row.m2PVDXY) #, muonJetBtag=row.m2JetBtag)

    def obj3_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_fr(row.tPt)

    def obj1_qcd_weight(self, row, ledleptonId='h2taucuts', subledleptonId=None):
        return frfits.highpt_mu_qcd_fr[ledleptonId](muonJetPt=max(row.m1JetPt, row.m1Pt), muonPt=row.m1Pt, muonJetCSVBtag=max(0, row.m1JetCSVBtag)) #, muonPVDXY=row.m1PVDXY) #, muonJetBtag=row.m1JetBtag)

    def obj2_qcd_weight(self, row, ledleptonId=None, subledleptonId='h2taucuts'):
        return frfits.lowpt_mu_qcd_fr[subledleptonId](muonJetPt=max(row.m2JetPt, row.m2Pt), muonPt=row.m2Pt, muonJetCSVBtag=max(0, row.m2JetCSVBtag)) #, muonPVDXY=row.m2PVDXY) #, muonJetBtag=row.m2JetBtag)

    def obj3_qcd_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return row.m1_t_SS

    ## def obj1_charge_flip(self, row):
    ##     return 0

    ## def obj2_charge_flip(self, row):
    ##     return 0

