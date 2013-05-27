'''

Analyze EET events for the WH analysis

'''

from EETauTree import EETauTree
import os
from WHAnalyzerBase import WHAnalyzerBase, quad, inv_mass
import ROOT
import mcCorrectors
import baseSelections as selections
import fakerate_functions as frfits
import optimizer
import math
import array

#mtr = frfits.mt_likelihood_ratio
################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeEET(WHAnalyzerBase):
    tree = 'eet/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeEET, self).__init__(tree, outfile, EETauTree, **kwargs)
        self.hfunc['subMTMass'] = lambda row, weight: (row.e2_t_Mass, weight) if row.e1MtToMET > row.e2MtToMET else (row.e1_t_Mass, weight) #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in WHAnalyzerBase.fill_histos later
	def make_both_barrel(attribute):
		def f_(row, weight):
			return ( getattr(row, attribute), weight ) if row.e1AbsEta < 1.48 and row.e2AbsEta < 1.48 else (0., 0.)
		return f_

	def make_both_endcap(attribute):
		def f_(row, weight):
			return ( getattr(row, attribute), weight ) if row.e1AbsEta >= 1.48 and row.e2AbsEta >= 1.48 else (0., 0.)
		return f_

	def make_mixed(attribute):
		def f_(row, weight):
			return ( getattr(row, attribute), weight ) \
			    if (row.e1AbsEta >= 1.48 and row.e2AbsEta < 1.48) or \
			    (row.e1AbsEta < 1.48 and row.e2AbsEta >= 1.48) \
		            else (0., 0.)
		return f_
	
		
        self.hfunc['pt_ratio' ] = lambda row, weight: (row.e2Pt/row.e1Pt, weight)
        self.hfunc["e*1_e2_Mass"] = lambda row, weight: ( frfits.mass_scaler['h2taucuts']( row.e1_e2_Mass), weight)
        self.hfunc["e*1_t_Mass" ] = lambda row, weight: ( frfits.mass_scaler['h2taucuts']( row.e1_t_Mass ), weight)
        self.hfunc["e1_e*2_Mass"] = lambda row, weight: ( frfits.mass_scaler['h2taucuts']( row.e1_e2_Mass), weight)
        self.hfunc["e*2_t_Mass" ] = lambda row, weight: ( frfits.mass_scaler['h2taucuts']( row.e2_t_Mass ), weight)
        self.hfunc["logic_cut_met" ] = self.logic_cut_met
	self.hfunc["my_selection_info" ] = self.fill_id_info

	self.hfunc["type1_pfMetEt_barr"] = make_both_barrel("type1_pfMetEt")
	self.hfunc["type1_pfMetEt_endc"] = make_both_endcap("type1_pfMetEt")
	self.hfunc["type1_pfMetEt_mix" ] = make_mixed("type1_pfMetEt")      
	self.hfunc["mva_metEt_barr"] = make_both_barrel("mva_metEt")
	self.hfunc["mva_metEt_endc"] = make_both_endcap("mva_metEt")
	self.hfunc["mva_metEt_mix" ] = make_mixed("mva_metEt")      


        self.hfunc["e2RelPFIsoDB_bar"] = lambda row, weight: ( row.e2RelPFIsoDB, weight ) if row.e2AbsEta < 1.48  else (0., 0.)
        self.hfunc["e1RelPFIsoDB_bar"] = lambda row, weight: ( row.e1RelPFIsoDB, weight ) if row.e1AbsEta < 1.48  else (0., 0.)
        self.hfunc["e2RelPFIsoDB_end"] = lambda row, weight: ( row.e2RelPFIsoDB, weight ) if row.e2AbsEta >= 1.48 else (0., 0.)
        self.hfunc["e1RelPFIsoDB_end"] = lambda row, weight: ( row.e1RelPFIsoDB, weight ) if row.e1AbsEta >= 1.48 else (0., 0.)


        self.hfunc["e1_t_CosThetaStar_barr"] = make_both_barrel("e1_t_CosThetaStar")
        self.hfunc["e1_t_CosThetaStar_endc"] = make_both_endcap("e1_t_CosThetaStar")
        self.hfunc["e1_t_CosThetaStar_mix" ] = make_mixed("e1_t_CosThetaStar")      
	self.hfunc["e1_e2_Mass_barr"] = make_both_barrel("e1_e2_Mass")
	self.hfunc["e1_e2_Mass_endc"] = make_both_endcap("e1_e2_Mass")
	self.hfunc["e1_e2_Mass_mix" ] = make_mixed("e1_e2_Mass")
        self.hfunc["electron_rejection_study" ] = self.electron_rejection_study
	self.hfunc["tau_id_study" ] = self.tau_id_study

        self.pucorrector = mcCorrectors.make_puCorrector('doublee')

    

    @staticmethod
    def logic_cut_met( row, weight):
        z_Mass_distance = row.e1_e2_Mass - frfits.mass_scaler['h2taucuts'](91.2) \
            if row.e1_e2_SS else \
            row.e1_e2_Mass - 91.2
        z_Mass_distance = abs(z_Mass_distance)
        if row.e1AbsEta < 1.48 and row.e2AbsEta < 1.48: #both barrel
            return False \
		    if row.mva_metEt < 25 and \
		    z_Mass_distance < 10 \
		    else True
        elif row.e1AbsEta < 1.48 or row.e2AbsEta < 1.48: #at least one in barrel
            return False \
		    if row.mva_metEt < 40 and \
		    z_Mass_distance < 10 \
		    else True
        else: #both in endcap
            return False \
		    if row.mva_metEt < 40 and \
		    z_Mass_distance < 20 \
	        else True
        return True

    @staticmethod
    def electron_rejection_study( row, weight):
        if row.tAntiElectronMVA3VTight:
            return (4.5, weight)
        elif row.tAntiElectronMVA3Tight:
            return (3.5, weight)
        elif row.tAntiElectronMVA3Medium:
            return (2.5, weight)
        elif row.tAntiElectronMVA3Loose:
            return (1.5, weight)
        return (0.5, weight)

    @staticmethod
    def tau_id_study( row, weight):
        if row.tTightIso3Hits:
            return (2.5, weight)
        elif row.tMediumIso3Hits:
            return (1.5, weight)
        return (0.5, weight)

    @staticmethod
    def fill_id_info(row, weight):
	return array.array("f", [row.e1_e2_Mass, row.e1AbsEta, row.e2AbsEta, row.type1_pfMetEt, row.e2_t_CosThetaStar, row.e1_t_CosThetaStar, row.e1_e2_CosThetaStar, weight] ), None

    def book_histos(self, folder):
	#PLOTS TO FILL IN ANY CASE
	for key in optimizer.grid_search:
	    prefix = key+'$' if key else ''
	    self.book(folder, prefix+"LT", "L_T", 100, 0, 300)
            self.book(folder, prefix+"e2_t_Mass", "subleadingMass", 200, 0, 200)
            self.book(folder, prefix+"e2_t_Pt",   "subleadingPt", 400, 0, 400)
            #Charge mis-id special histograms
            if 'c2' in folder:
                self.book(folder, prefix+"e*2_t_Mass", "subleadingMass with misid sclaing correction", 200, 0, 200)

        if len(optimizer.grid_search.keys()) == 1:
            if 'c1' in folder:
                self.book(folder, "e*1_e2_Mass", "E 1-2 Mass with misid sclaing correction", 120, 0, 120)
                self.book(folder, "e*1_t_Mass", "leadingMass with misid sclaing correction", 200, 0, 200)
            elif 'c2' in folder:
                self.book(folder, "e1_e*2_Mass", "E 1-2 Mass with misid sclaing correction", 120, 0, 120)
            self.book(folder, "e2RelPFIsoDB", "e2RelPFIsoDB", 30, 0, 0.3)
            self.book(folder, "e1RelPFIsoDB", "e1RelPFIsoDB", 30, 0, 0.3)
            self.book(folder, "e1_e2_Mass", "E 1-2 Mass", 120, 0, 120)
            self.book(folder, "e1_t_Mass", "leadingMass", 200, 0, 200)
            self.book(folder, "weight", "Event weight", 100, 0, 5)
            self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
            self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
            self.book(folder, "e1Pt", "E 1 Pt", 100, 0, 100)
            self.book(folder, "e2Pt", "E 2 Pt", 100, 0, 100)
            self.book(folder, "e1AbsEta", "Muon 1 AbsEta", 100, 0, 2.5)
            self.book(folder, "e2AbsEta", "Muon 2 AbsEta", 100, 0, 2.5)
	    self.book(folder, "tPt", "tPt", 100, 0,100)
	    self.book(folder, "tAbsEta", "tAbsEta", 100, 0, 2.3)

            self.book(folder, "e2RelPFIsoDB_bar", "e2RelPFIsoDB", 30, 0, 0.3)
            self.book(folder, "e1RelPFIsoDB_bar", "e1RelPFIsoDB", 30, 0, 0.3)
            self.book(folder, "e2RelPFIsoDB_end", "e2RelPFIsoDB", 30, 0, 0.3)
            self.book(folder, "e1RelPFIsoDB_end", "e1RelPFIsoDB", 30, 0, 0.3)
            self.book(folder, "tau_id_study", "tau_id_study", 3, 0, 3)
            self.book(folder, "electron_rejection_study", "electron_rejection_study", 5, 0, 5)

            #let's look for osme other possible selections
            self.book(folder, "Mass"          , "Mass"      , 100, 0, 1)
            self.book(folder, "pt_ratio"      , "pt_ratio"      , 100, 0, 1)
            self.book(folder, "tToMETDPhi"    , "tToMETDPhi"    , 100, 0, 4)
            self.book(folder, "e1_e2_Pt"       , "lepRecoil"     , 600, 0, 8000)
            self.book(folder, "e1_e2_DR"       , "e1_e2_DR"      , 500, 0, 10)
            self.book(folder, "e1_e2_CosThetaStar", "e1_e2_CosThetaStar", 110, 0., 1.1)
            self.book(folder, "e1_t_CosThetaStar" , "e1_t_CosThetaStar" , 110, 0., 1.1)
            self.book(folder, "e2_t_CosThetaStar" , "e2_t_CosThetaStar" , 110, 0., 1.1)
            self.book(folder, "type1_pfMetEt"         , "metEt"         , 300, 0, 300)
            self.book(folder, "mva_metEt"         , "metEt"         , 300, 0, 300)

	    #split into both e in barr/endcap/mixed
	    self.book(folder, "type1_pfMetEt_barr", "metEt", 300, 0, 300)
	    self.book(folder, "type1_pfMetEt_endc", "metEt", 300, 0, 300)
	    self.book(folder, "type1_pfMetEt_mix" , "metEt", 300, 0, 300)

	    self.book(folder, "mva_metEt_barr", "metEt", 300, 0, 300)
	    self.book(folder, "mva_metEt_endc", "metEt", 300, 0, 300)
	    self.book(folder, "mva_metEt_mix" , "metEt", 300, 0, 300)

            self.book(folder, "e1_t_CosThetaStar_barr", "e1_t_CosThetaStar" , 110, 0., 1.1)
            self.book(folder, "e1_t_CosThetaStar_endc", "e1_t_CosThetaStar" , 110, 0., 1.1)
            self.book(folder, "e1_t_CosThetaStar_mix" , "e1_t_CosThetaStar" , 110, 0., 1.1)

	    self.book(folder, "e1_e2_Mass_barr", "E 1-2 Mass", 120, 0, 120)
	    self.book(folder, "e1_e2_Mass_endc", "E 1-2 Mass", 120, 0, 120)
	    self.book(folder, "e1_e2_Mass_mix" , "E 1-2 Mass", 120, 0, 120)
            
    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    def preselection(self, row, cut_flow_trk = None, LT_threshold = 80.):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.doubleEPass:                   return False
	if not (row.e1MatchesDoubleEPath > 0 and \
		row.e2MatchesDoubleEPath > 0): return False 
        cut_flow_trk.Fill('trigger')

        if not selections.eSelection(row, 'e1'):  return False
        cut_flow_trk.Fill('obj1 Presel')

        if not selections.eSelection(row, 'e2'):  return False
        cut_flow_trk.Fill('obj2 Presel')

        if not selections.tauSelection(row, 't'): return False
        if not row.tAntiMuonLoose:   return False
        cut_flow_trk.Fill('obj3 Presel')

        if row.LT < LT_threshold: return False
        cut_flow_trk.Fill('LT')

        if row.e1_e2_SS and row.e1_t_SS         : return False #remove three SS leptons
        if row.e1_e2_Mass < 20:                   return False
        if not selections.vetos(row):             return False #applies mu bjet e additional tau vetoes
        cut_flow_trk.Fill('vetos')

        #REMOVE CHARGE FAKES!
        if not self.logic_cut_met(row, 1.): return False
        cut_flow_trk.Fill('charge_fakes')  

        return True

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def sign_cut(row):
        ''' Returns true if muons are SS '''
        return bool(row.e1_e2_SS)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj1_id(row, leadleptonId='h2taucuts', subleadleptonId=None):
        return selections.lepton_id_iso(row, 'e1', leadleptonId)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj2_id(row, leadleptonId=None, subleadleptonId='h2taucuts'):
        return selections.lepton_id_iso(row, 'e2', subleadleptonId)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def obj3_id(row):
        return bool(row.tLooseIso3Hits)

    #There is no call to self, so just promote it to statucmethod, to allow usage by other dedicated analyzers
    @staticmethod
    def anti_wz(row):
        if row.e2_t_Zcompat < 10 or row.e1_t_Zcompat < 10 :
            if not row.tAntiElectronMVA3Tight:
                return False
        if row.e2_t_Zcompat < 20 or row.e1_t_Zcompat < 20 :
            if not row.tAntiElectronMVA3Medium:
                return False
        return row.tAntiElectronMVA3Loose

    def enhance_wz(self, row):
        # Require the "tau" to be a electron, and require the third electron
        # to have M_Z +- 20 
        if self.anti_wz(row):
            return False
        # Make sure any Z is from e1
        e2_good_Z = bool(71 < row.e2_t_Mass < 111)
        return not e2_good_Z

    def event_weight(self, row):
        if row.run > 2:
            return 1.
        return self.pucorrector(row.nTruePU) * \
            mcCorrectors.get_electron_corrections(row,'e1','e2')

   
    def obj1_weight(self, row, leadleptonId='h2taucuts', subleadleptonId=None):
        return frfits.highpt_ee_fr[leadleptonId](max(row.e1JetPt, row.e1Pt))

    def obj2_weight(self, row, leadleptonId=None, subleadleptonId='h2taucuts'):
	return frfits.lowpt_ee_fr[subleadleptonId](max(row.e1JetPt, row.e1Pt))

    def obj3_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_fr(row.tPt)

    def obj1_qcd_weight(self, row, leadleptonId='h2taucuts', subleadleptonId=None):
        return frfits.highpt_ee_qcd_fr[leadleptonId](max(row.e1JetPt, row.e1Pt))

    def obj2_qcd_weight(self, row, leadleptonId=None, subleadleptonId='h2taucuts'):
        return frfits.lowpt_ee_qcd_fr[subleadleptonId](max(row.e2JetPt, row.e2Pt))

    def obj3_qcd_weight(self, row, notUsed1=None, notUsed2=None):
        return frfits.tau_qcd_fr(row.tPt)

    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return not row.e1_t_SS

    def obj1_charge_flip(self, row, leadleptonId='h2taucuts', subleadleptonId=None):
        return frfits.e_charge_flip[leadleptonId](row.e1AbsEta,row.e1Pt) #highpt_e_charge_flip

    def obj2_charge_flip(self, row, leadleptonId=None, subleadleptonId='h2taucuts'):
        return frfits.e_charge_flip[subleadleptonId](row.e2AbsEta,row.e2Pt) #lowpt_e_charge_flip

    def obj1_charge_flip_sysup(self, row, leadleptonId='h2taucuts', subleadleptonId=None):
        return frfits.e_charge_flip_up[leadleptonId](row.e1AbsEta,row.e1Pt) #highpt_e_charge_flip

    def obj2_charge_flip_sysup(self, row, leadleptonId=None, subleadleptonId='h2taucuts'):
        return frfits.e_charge_flip_up[subleadleptonId](row.e2AbsEta,row.e2Pt) #lowpt_e_charge_flip







