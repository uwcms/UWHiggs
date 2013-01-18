'''

Analyze EET events for the WH analysis

'''
'''

Analyze MMT events for the WH analysis

'''

from EETauTree import EETauTree
import os
import WHAnalyzerBase
import ROOT
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits

################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeEET(WHAnalyzerBase.WHAnalyzerBase):
    tree = 'eet/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeEET, self).__init__(tree, outfile, EETauTree, **kwargs)
        self.hfunc['subMTMass'] = lambda row, weight: (row.m2_t_Mass, weight) if row.m1MtToMET > row.m2MtToMET else (row.m1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
        self.pucorrector = mcCorrectors.make_puCorrector('doublee')

    def book_histos(self, folder):
        self.book(folder, "weight", "Event weight", 100, 0, 5)
        #self.book(folder, "weight_nopu", "Event weight without PU", 100, 0, 5)
        self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
        self.book(folder, "prescale", "HLT prescale", 21, -0.5, 20.5)

        self.book(folder, "e1Pt", "E 1 Pt", 100, 0, 100)
        self.book(folder, "e2Pt", "E 2 Pt", 100, 0, 100)

        self.book(folder, "e1AbsEta", "Muon 1 AbsEta", 100, 0, 2.4)
        self.book(folder, "e2AbsEta", "Muon 2 AbsEta", 100, 0, 2.4)

        self.book(folder, "e1_e2_Mass", "E 1-2 Mass", 120, 0, 120)
        self.book(folder, "e1_t_Mass", "leadingMass", 200, 0, 200)
        self.book(folder, "e2_t_Mass", "subleadingMass", 200, 0, 200)

        self.book(folder, "e2RelPFIsoDB", "e2RelPFIsoDB", 100, 0, 0.3)
        self.book(folder, "tPt", "tPt", 100, 0,100)
        self.book(folder, "tAbsEta", "tAbsEta", 100, 0, 2.3)
        self.book(folder, "metSignificance", "MET significance", 100, 0, 15)
        self.book(folder, "LT", "L_T", 100, 0, 300)

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        #object basic selections
        if not row.doubleEPass:                   return False
        if row.e1Pt < 20:                         return False
        if not selections.eSelection(row, 'e1'):  return False
        if not selections.eSelection(row, 'e2'):  return False
        if not selections.tauSelection(row, 't'): return False

        if row.e1_e2_Mass < 20:       return False
        if row.metSignificance < 2.5: return False
        if row.e1_e2_Mass > 81 \
            and row.e1_e2_Mass < 101: return False
        if row.LT < 100:              return False
        if not selections.vetos(row): return False #applies mu bjet e additional tau vetoes

        if not row.tAntiMuonTight:   return False
        if row.tMuOverlap:           return False
        #'t_ElectronOverlapWP95 < 0.5',
        return True

    def sign_cut(self, row):
        ''' Returns true if muons are SS '''
        return bool(row.e1_e2_SS)

    def obj1_id(self, row):
        return bool(row.e1MVAIDH2TauWP) and bool( row.e1RelPFIsoDB < 0.1 or (row.e1RelPFIsoDB < 0.15 and row.e1AbsEta < 1.479))

    def obj2_id(self, row):
        return bool(row.e2MVAIDH2TauWP) and bool( row.e2RelPFIsoDB < 0.1 or (row.e2RelPFIsoDB < 0.15 and row.e2AbsEta < 1.479))

    def obj3_id(self, row):
        return bool(row.tLooseMVAIso)

    def anti_wz(self, row):
        return row.tAntiElectronMVA and not row.tCiCTightElecOverlap

    def enhance_wz(self, row):
        # Require the "tau" to be a electron, and require the third electron
        # to have M_Z +- 20 
        if self.anti_wz(row):
            return False
        # Make sure any Z is from m1
        e2_good_Z = bool(71 < row.e2_t_Mass < 111)
        return not e2_good_Z

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_electron_corrections(row,'e1','e2')

    def obj1_weight(self, row):
        return highpt_ee_fr(row.e1JetPt)

    def obj2_weight(self, row):
        return lowpt_ee_fr(row.e2JetPt)

    
    def obj1_weight(self, row):
        return frfits.highpt_mu_fr(max(row.m1JetPt, row.m1Pt))

    def obj2_weight(self, row):
        return frfits.lowpt_mu_fr(max(row.m2JetPt, row.m2Pt))

    def obj3_weight(self, row):
        return frfits.tau_fr(row.tPt)

    def obj1_qcd_weight(self, row):
        return frfits.highpt_mu_qcd_fr(max(row.m1JetPt, row.m1Pt))

    def obj2_qcd_weight(self, row):
        return frfits.lowpt_mu_qcd_fr(max(row.m2JetPt, row.m2Pt))

    def obj3_qcd_weight(self, row):
        return frfits.tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return not row.m1_t_SS

    def obj1_charge_flip(self, row):
        return frfits.e_charge_flip(row.e1AbsEta,row.e1Pt)









