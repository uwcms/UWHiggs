'''

Test: LFV H->MuTau Analysis Skeleton

Authors: M.C, Aaron Levine, Evan K. Friis, UW

'''

import MuTauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import ROOT
import math

###### Because I need to add a bunch more branches to the ntuple...
from math import sqrt, pi

def deltaPhi(phi1, phi2):
  PHI = abs(phi1-phi2)
  if PHI<=pi:
      return PHI
  else:
      return 2*pi-PHI


def fullMT(met,mupt,taupt, metphi, muphi, tauphi):
	mux=mupt*math.cos(muphi)
        muy=mupt*math.sin(muphi)
        metx=met*math.cos(metphi)
        mety=met*math.sin(metphi)
        taux=taupt*math.cos(tauphi)
        tauy=taupt*math.sin(tauphi)
	full_et=met+mupt+taupt   # for muon and tau I am approximating pt~et (M<<P)
	full_x=metx+mux+taux
        full_y=mety+muy+tauy
	full_mt_2 = full_et*full_et-full_x*full_x-full_y*full_y
	full_mt=0
	if (full_mt_2>0):
		full_mt= math.sqrt(full_mt_2)
	return full_mt

#def fullMass(met,mupt,taupt, metphi, muphi, mupz, tauphi,taupz):
#        mux=mupt*math.cos(muphi)
#        muy=mupt*math.sin(muphi)
#        metx=met*math.cos(metphi)
#        mety=met*math.sin(metphi)
#        taux=taupt*math.cos(tauphi)
#        tauy=taupt*math.sin(tauphi)
#	metz=taupz # this is just an assumption!
#	full_x=metx+mux+taux
#        full_y=mety+muy+tauy
#        full_z=metz+muz+tauz
#	full_p_2=full_x*full_x+full_y*full_y+full_z*full_z
#


################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################

# Determine MC-DATA corrections
is7TeV = bool('7TeV' in os.environ['jobid'])
print "Is 7TeV:", is7TeV

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
pu_distributions = glob.glob(os.path.join(
	'inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))
pu_corrector = PileupWeight.PileupWeight(
    'S6' if is7TeV else 'S10', *pu_distributions)

muon_pog_PFTight_2012 = MuonPOGCorrections.make_muon_pog_PFTight_2012ABCD()
muon_pog_PFRelIsoDB012_2012 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB012_2012ABCD()
muon_pog_IsoMu24eta2p1_2012 = MuonPOGCorrections.make_muon_pog_IsoMu24eta2p1_2012()



def mc_corrector_2012(row):
    if row.run > 2:
        return 1
    pu = pu_corrector(row.nTruePU)
    m1id = muon_pog_PFTight_2012(row.mPt, row.mEta)
    m1iso = muon_pog_PFRelIsoDB012_2012(row.mPt, row.mEta)
    m_trg =  muon_pog_IsoMu24eta2p1_2012(row.mPt, row.mAbsEta) 
    return pu*m1id*m1iso*m_trg

mc_corrector = mc_corrector_2012

class AnalyzeMuTauTightPFMET(MegaBase):
    tree = 'mt/final/Ntuple'

    def __init__(self, tree, outfile, **kwargs):
        super(AnalyzeMuTauTightPFMET, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuTauTree.MuTauTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):

        names=["gg","antiisomuongg","antiisotaugg","ssgg", "ssantiisotaugg", "ssantiisomuongg","ggBaselineTightTauIso", "ggBaseline","antiisomuonggBaselineTightTauIso","antiisomuonggBaseline", "checkZs", 
"checkZsTauIDTightIso", "checkZsTauIDLooseIso", "checkZsTauID"
]

	for x in range(0,14):

            self.book(names[x], "weight", "Event weight", 100, 0, 5)
            self.book(names[x], "weight_nopu", "Event weight without PU", 100, 0, 5)
    
            self.book(names[x], "rho", "Fastjet #rho", 100, 0, 25)
            self.book(names[x], "nvtx", "Number of vertices", 31, -0.5, 30.5)
            self.book(names[x], "prescale", "HLT prescale", 21, -0.5, 20.5)
    
            self.book(names[x], "mPt", "Muon  Pt", 500,0,500)
            self.book(names[x], "mEta", "Muon  eta", 200, -2.5, 2.5)
            self.book(names[x], "mMtToMVAMET", "Muon MT (MVA)", 500,0,500)
            self.book(names[x], "mMtToPfMet_Ty1", "Muon MT (PF Ty1)", 500,0,500)
            self.book(names[x], "mCharge", "Muon Charge", 5, -2, 2)
            self.book(names[x], "mRelPFIsoDB", "Muon Isolation", 1000, 0, 1)
   
            self.book(names[x], "tPt", "Tau  Pt", 500,0,500)
            self.book(names[x], "tJetPt", "Tau  Jet Pt", 500,0,500)
            self.book(names[x], "tEta", "Tau  eta", 200, -2.5, 2.5)
            self.book(names[x], "tMtToMVAMET", "Tau MT (MVA)", 500,0,500)
            self.book(names[x], "tMtToPfMet_Ty1", "Tau MT (PF Ty1)", 500,0,500)
            self.book(names[x], "tCharge", "Tau  Charge", 5, -2, 2)
            self.book(names[x], "tLooseIso", "Tau Isolation Loose", 3, 0, 2)
            self.book(names[x], "tTightIso3Hits", "Tau Isolation Tight", 3, 0, 2)
            self.book(names[x], "tDecayMode", "Tau Decay Mode", 21, -0.5, 20.5)
 
            self.book(names[x], 'mPixHits', 'Mu 1 pix hits', 10, -0.5, 9.5)
            self.book(names[x], 'mJetBtag', 'Mu 1 JetBtag', 100, -5.5, 9.5)
    
    	    self.book(names[x], "LT", "ht", 500, 0, 500)
   	    self.book(names[x],"fullMT_mva","fullMT_mva",500,0,500); 
            self.book(names[x],"fullMT_type1","fullMT_type1",500,0,500);

            self.book(names[x], "type1_pfMetEt", "Type1 MET", 500,0,500)
            self.book(names[x], "mva_metEt", "MVA MET", 500,0,500)
   
            self.book(names[x], "pfMetPhi","Type1 MET Phi",100,-3.15,3.15)
            self.book(names[x], "mva_metPhi","MVA MET Phi",100,-3.15,3.15)
 
            self.book(names[x], "m_t_Mass", "Muon + Tau Mass", 500,0,500)
            self.book(names[x], "m_t_Pt", "Muon + Tau Pt", 500,0,500)
            self.book(names[x], "m_t_DR", "Muon + Tau DR", 100, 0, 10)
            self.book(names[x], "m_t_DPhi", "Muon + Tau DPhi", 100, 0, 4)
            self.book(names[x], "m_t_SS", "Muon + Tau SS", 5, -2, 2)
            self.book(names[x], "m_t_ToMETDPhi_Ty1", "Muon Tau DPhi to MET", 100, 0, 4)
    
            # Vetoes
            self.book(names[x], 'bjetVeto', 'Number of b-jets', 5, -0.5, 4.5)
            self.book(names[x], 'bjetCSVVeto', 'Number of b-jets', 5, -0.5, 4.5)
            self.book(names[x], 'muVetoPt5IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
            self.book(names[x], 'muVetoPt15IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
            self.book(names[x], 'tauVetoPt20', 'Number of extra taus', 5, -0.5, 4.5)
            self.book(names[x], 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
    
            self.book(names[x], 'jetVeto20', 'Number of extra jets', 5, -0.5, 4.5)
            self.book(names[x], 'jetVeto30', 'Number of extra jets', 5, -0.5, 4.5)
    
            self.book(names[x], "mPhiMtPhi", "", 100, 0,4)
            self.book(names[x], "mPhiMETPhiMVA", "", 100, 0,4)
            self.book(names[x], "tPhiMETPhiMVA", "", 100, 0,4)
            self.book(names[x], "mPhiMETPhiType1", "", 100, 0,4)
            self.book(names[x], "tPhiMETPhiType1", "", 100, 0,4)
    

    def correction(self, row):
        return mc_corrector(row)

    def fill_histos(self, row,name='gg'):
        histos = self.histograms
        weight = self.correction(row)
        histos[name+'/weight'].Fill(weight)
        histos[name+'/weight_nopu'].Fill(self.correction(row))
        histos[name+'/rho'].Fill(row.rho, weight)
        histos[name+'/nvtx'].Fill(row.nvtx, weight)
        histos[name+'/prescale'].Fill(row.doubleMuPrescale, weight)

        histos[name+'/mPt'].Fill(row.mPt, weight)
        histos[name+'/mEta'].Fill(row.mEta, weight)
	histos[name+'/mMtToMVAMET'].Fill(row.mMtToMVAMET,weight)
        histos[name+'/mMtToPfMet_Ty1'].Fill(row.mMtToPfMet_Ty1,weight)
        histos[name+'/mCharge'].Fill(row.mCharge, weight)
        histos[name+'/mRelPFIsoDB'].Fill(row.mRelPFIsoDB, weight)

        histos[name+'/tPt'].Fill(row.tPt, weight)
        histos[name+'/tJetPt'].Fill(row.tJetPt, weight)
        histos[name+'/tEta'].Fill(row.tEta, weight)
        histos[name+'/tMtToMVAMET'].Fill(row.tMtToMVAMET,weight)
        histos[name+'/tMtToPfMet_Ty1'].Fill(row.tMtToPfMet_Ty1,weight)
        histos[name+'/tCharge'].Fill(row.tCharge, weight)
        histos[name+'/tLooseIso'].Fill(row.tLooseIso, weight)
        histos[name+'/tTightIso3Hits'].Fill(row.tTightIso3Hits, weight)
        histos[name+'/tDecayMode'].Fill(row.tDecayMode, weight)

	histos[name+'/LT'].Fill(row.LT,weight)
        histos[name+'/fullMT_mva'].Fill(fullMT(row.mva_metEt,row.mPt,row.tPt,row.mva_metPhi, row.mPhi, row.tPhi),weight)
        histos[name+'/fullMT_type1'].Fill(fullMT(row.type1_pfMetEt,row.mPt,row.tPt,row.type1_pfMetPhi, row.mPhi, row.tPhi),weight)

	histos[name+'/type1_pfMetEt'].Fill(row.type1_pfMetEt,weight)
	histos[name+'/mva_metEt'].Fill(row.mva_metEt,weight)

        histos[name+'/pfMetPhi'].Fill(row.pfMetPhi,weight)
        histos[name+'/mva_metPhi'].Fill(row.mva_metPhi,weight)

        histos[name+'/m_t_Mass'].Fill(row.m_t_Mass,weight)
        histos[name+'/m_t_Pt'].Fill(row.m_t_Pt,weight)
        histos[name+'/m_t_DR'].Fill(row.m_t_DR,weight)
        histos[name+'/m_t_DPhi'].Fill(row.m_t_DPhi,weight)
        histos[name+'/m_t_SS'].Fill(row.m_t_SS,weight)
	histos[name+'/m_t_ToMETDPhi_Ty1'].Fill(row.m_t_ToMETDPhi_Ty1,weight)

        histos[name+'/mPixHits'].Fill(row.mPixHits, weight)
        histos[name+'/mJetBtag'].Fill(row.mJetBtag, weight)

        histos[name+'/bjetVeto'].Fill(row.bjetVeto, weight)
        histos[name+'/bjetCSVVeto'].Fill(row.bjetCSVVeto, weight)
        histos[name+'/muVetoPt5IsoIdVtx'].Fill(row.muVetoPt5IsoIdVtx, weight)
        histos[name+'/muVetoPt15IsoIdVtx'].Fill(row.muVetoPt15IsoIdVtx, weight)
        histos[name+'/tauVetoPt20'].Fill(row.tauVetoPt20, weight)
        histos[name+'/eVetoCicTightIso'].Fill(row.eVetoCicTightIso, weight)
        histos[name+'/jetVeto20'].Fill(row.jetVeto20, weight)
        histos[name+'/jetVeto30'].Fill(row.jetVeto30, weight)

        histos[name+'/mPhiMtPhi'].Fill(deltaPhi(row.mPhi,row.tPhi),weight)
        histos[name+'/mPhiMETPhiMVA'].Fill(deltaPhi(row.mPhi,row.mva_metPhi),weight)
        histos[name+'/tPhiMETPhiMVA'].Fill(deltaPhi(row.tPhi,row.mva_metPhi),weight)
        histos[name+'/mPhiMETPhiType1'].Fill(deltaPhi(row.mPhi,row.type1_pfMetPhi),weight)
        histos[name+'/tPhiMETPhiType1'].Fill(deltaPhi(row.tPhi,row.type1_pfMetPhi),weight)







    def baselineCheck(self, row):
        if not row.isoMu24eta2p1Pass:
            return False
        if row.mPt < 26:
            return False
        if abs(row.mEta) >= 2.1:
            return False
        if row.tPt<20 :
            return False
        if abs(row.tEta)>=2.3 :
            return False
        if row.mCharge*row.tCharge!=-1:
            return False
        if not bool(row.mPFIDTight):
	    return False 
        if not bool(row.mRelPFIsoDB <0.12):
	    return False 
        if not bool(abs(row.mDZ) < 0.2):
	    return False
	if  not row.tLooseIso:
	    return False
	if not row.tAntiElectronLoose:
		return False
	if not row.tAntiMuonTight:
		return False
	if not row.tDecayFinding:
		return False
	if not bool (row.muVetoPt5IsoIdVtx<1):
		return False
	if not bool (row.eVetoCicTightIso<1):
		return False 
        return True


    def presel(self, row):
	if not row.isoMu24eta2p1Pass:
            return False
        if row.mPt < 26:
            return False
        if abs(row.mEta) >= 2.1:
            return False
        return True

    def kinematics(self, row):
        if row.mPt<30 :
            return False
        if row.tPt<20 :
            return False
        if abs(row.tEta)>=2.3 :
            return False
	if row.mMtToPfMet_Ty1<30:
	    return False
        return True

    def ggtightKIN(self,row):
       if row.mPt < 40:
           return False
       if row.tPt < 40:
           return False
       if row.mMtToPfMet_Ty1<30:
            return False
       if row.tMtToPfMet_Ty1>20:
            return False
       return True	

    def ggtightANGLE(self,row):
       if (deltaPhi(row.tPhi,row.pfMetPhi)>0.3):
          return False
       if (deltaPhi(row.mPhi,row.tPhi)<2.):    
          return False
       if (deltaPhi(row.mPhi,row.pfMetPhi)<2):
          return False
       return True

    def ggtight(self,row):
       if row.mPt < 26:
           return False
       if row.tPt < 40:
           return False
       if row.mMtToPfMet_Ty1<30:
            return False
       if row.tMtToPfMet_Ty1>20:
            return False
       if (deltaPhi(row.tPhi,row.pfMetPhi)>0.3):
          return False
       if (deltaPhi(row.mPhi,row.tPhi)<2.5):    
          return False
       if (deltaPhi(row.mPhi,row.pfMetPhi)<2.5):
          return False
       return True

    def njets(self,row):
	return bool(row.jetVeto30<1)

    def vbf(self,row):
#	if(row.vbfNJets<2):
#	    return False
	if(abs(row.vbfDeta)<3):
	    return False
        return True

    def oppositesign(self,row):
	if row.mCharge*row.tCharge!=-1:
            return False
	return True

    def obj1_id(self, row):
        return bool(row.mPFIDTight)  and bool(abs(row.mDZ) < 0.2) 

    def obj2_id(self, row):
	return  row.tAntiElectronLoose and row.tAntiMuonTight and row.tDecayFinding

    def vetos(self,row):
	return bool (row.muVetoPt5IsoIdVtx<1) and bool (row.eVetoCicTightIso<1) 

    def obj1_iso(self, row):
        return bool(row.mRelPFIsoDB <0.12)

    def obj2_iso(self, row):
	return  row.tLooseIso

    def obj2_iso_tight(self, row):
        return  row.tTightIso3Hits

    def obj1_antiiso(self, row):
        return bool(row.mRelPFIsoDB >0.2) 

    def obj2_antiiso(self, row):
        return  not row.tTightIso3Hits

    def selectZs(self, row):
	return bool(row.muVetoPt15IsoIdVtx==1)

    def process(self):
        event =0
        sel=False
        for row in self.tree:
            if event!=row.evt:
                event=row.evt
                sel = False
            if sel==True:
                continue

	    if not self.presel(row):
		continue
            if not self.obj1_id(row):
                continue
            if not self.kinematics(row):
                continue
            if self.obj1_iso(row) and self.selectZs(row):
                self.fill_histos(row,'checkZs')
		if self.obj2_id(row):
                 self.fill_histos(row,'checkZsTauID')
		if self.obj2_id(row) and self.obj2_iso(row):
		 self.fill_histos(row,'checkZsTauIDLooseIso')
		if self.obj2_id(row) and self.obj2_iso_tight(row):
		 self.fill_histos(row,'checkZsTauIDTightIso')
            if not self.vetos(row):
                continue
	    sel = True
            if self.obj1_iso(row) and self.oppositesign(row) and  row.tAntiElectronLoose and row.tAntiMuonTight:
                        self.fill_histos(row,'ggBaselineNOTauIsoRelTauID')
            if not self.obj2_id(row):
                continue

	    if self.obj1_iso(row) and self.obj2_iso(row) and self.oppositesign(row):
                        self.fill_histos(row,'ggBaseline')
	    if self.obj1_antiiso(row) and self.obj2_iso(row) and  self.oppositesign(row):
                        self.fill_histos(row,'antiisomuonggBaseline')

            if self.obj1_iso(row) and self.obj2_iso_tight(row) and self.oppositesign(row):
                        self.fill_histos(row,'ggBaselineTightTauIso')
            if self.obj1_antiiso(row) and self.obj2_iso_tight(row) and  self.oppositesign(row):
                        self.fill_histos(row,'antiisomuonggBaselineTightTauIso')

            if self.obj1_iso(row) and self.oppositesign(row):
                        self.fill_histos(row,'ggBaselineNOTauIso')

	    if self.ggtight(row):
		if self.obj1_iso(row) and self.obj2_iso_tight(row) and self.oppositesign(row): 
	                self.fill_histos(row,'gg')
                if self.obj1_iso(row) and self.oppositesign(row):
                        self.fill_histos(row,'notauisogg')
                if self.obj1_iso(row) and self.obj2_iso_tight(row) and not self.oppositesign(row):
                        self.fill_histos(row,'ssgg')
                if self.obj1_antiiso(row) and self.obj2_iso_tight(row) and  self.oppositesign(row):
                        self.fill_histos(row,'antiisomuongg')
                if self.obj1_iso(row) and self.obj2_antiiso(row) and  self.oppositesign(row):
                        self.fill_histos(row,'antiisotaugg')
		if self.obj1_antiiso(row) and self.obj2_iso_tight(row) and not self.oppositesign(row):
                        self.fill_histos(row,'ssantiisomuongg')
                if self.obj1_iso(row) and not self.obj2_iso(row) and not self.oppositesign(row):
                        self.fill_histos(row,'ssantiisotaugg')


    def finish(self):
        self.write_histos()