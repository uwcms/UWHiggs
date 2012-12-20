import ROOT
from ROOT import TVector3 # for delta R

from MOOSEY.OrderedDict import OrderedDict
from MOOSEY.cuts import CutflowDecision
from math import fabs,pi


#define the selection of a good event
def passes_mumu_trigger_data(HLT_Fired,HLT_Prescale,run):
    trigger_fired = HLT_Fired
    prescale = HLT_Prescale

    return (trigger_fired == 1 and prescale == 1)

#define the selection of a good event
def passes_mumu_trigger_mc(HLT_Fired,HLT_Prescale,evt_frac):
    trigger_fired = HLT_Fired
    prescale = HLT_Prescale
    
    return (trigger_fired == 1 and prescale == 1)

def passes_ee_trigger_data(HLT_Fired,HLT_Prescale,HLT_Index,run):
    trigger_fired = 0
    trigger_prescale = 0
    index = -1
    if run >= 160431 and run <= 167913:
        index = HLT_Index[83]
        trigger_fired = HLT_Fired[HLT_Index[83]]
        prescale = HLT_Prescale[HLT_Index[83]]
    elif run >= 170249 and run <= 180252:
        index = HLT_Index[86]
        trigger_fired = HLT_Fired[HLT_Index[86]]
        prescale = HLT_Prescale[HLT_Index[86]]

    return (trigger_fired == 1 and prescale == 1 and index > 0)

def passes_ee_trigger_mc(HLT_Fired,HLT_Prescale,HLT_Index,evt_frac):
    trigger_fired = 0
    trigger_prescale = 0
    index = -1
    #if evt_frac < 0.236:
    #    index = HLT_Index[83]
    #    trigger_fired = HLT_Fired[HLT_Index[83]]
    #    prescale = HLT_Prescale[HLT_Index[83]]
    #else:
    index = HLT_Index[86]
    trigger_fired = HLT_Fired[HLT_Index[86]]
    prescale = HLT_Prescale[HLT_Index[86]]

    return (trigger_fired == 1 and prescale == 1 and index > 0)
    
def good_vtx(nGoodVtx):    
    return nGoodVtx != 0

def no_scraping(IsTracksGood):
    return IsTracksGood != 0

mumu_good_event_data_reqs = OrderedDict([['trigger',
                                          [passes_mumu_trigger_data,
                                           ['mu17mu8Pass','mu17mu8Prescale','run'],
                                           0]],
                                         ['goodvtx',
                                          [good_vtx,
                                           ['nvtx'],
                                           0]]])

mumu_good_event_mc_reqs = OrderedDict([['trigger',
                                        [passes_mumu_trigger_mc,
                                         ['mu17mu8Pass','mu17mu8Prescale','eventFraction'],
                                         0]],
                                       ['goodvtx',
                                        [good_vtx,
                                         ['nvtx'],
                                         0]]])

ee_good_event_data_reqs = OrderedDict([['trigger',
                                        [passes_ee_trigger_data,
                                         ['HLT','HLTprescale','HLTIndex','run'],
                                         0]],
                                       ['goodvtx',
                                        [good_vtx,
                                         ['nGoodVtx'],
                                         0]],
                                       ['noscraping',
                                        [no_scraping,
                                         ['IsTracksGood'],
                                         0]]])
ee_good_event_mc_reqs = OrderedDict([['trigger',
                                      [passes_ee_trigger_mc,
                                       ['HLT','HLTprescale','HLTIndex','eventFraction'],
                                       0]],
                                     ['goodvtx',
                                      [good_vtx,
                                       ['nGoodVtx'],
                                       0]],
                                     ['noscraping',
                                      [no_scraping,
                                       ['IsTracksGood'],
                                       0]]])

#single muon selection
def mu_pt(mu):
    return (mu.Pt() > 10)

def mu_eta(mu):
    return (abs(mu.Eta()) < 2.4)

#validated the use of IDHZG2011/2012
#just needs simple converter
def mu_id(the_id):
    bit = (the_id == 1.0)
    return bit

def mu_iso(mu,pfchg,pfneut,pfpho,ea,rho):    
    rhoprime = max(rho,0)
    eff_iso = max( pfneut + pfpho - ea*rhoprime, 0.0 )
    iso = (pfchg + eff_iso)/mu.Pt()
    return (iso < 0.12)

def mu_trg_match_data(trigMatches,run):
    hasMatch = True    
    #if run >= 160431 and run <= 163869:
    #    hasMatch = (trigMatches[13] > 0)
    #elif run >= 165088 and run <= 178380:
    #    hasMatch = (trigMatches[14] > 0)
    #elif run >= 178420 and run <= 180252:
    #    hasMatch = (trigMatches[15] > 0)
    return hasMatch

def mu_trg_match_mc(trigMatches,evt_frac):
    hasMatch = True
    #if evt_frac < 0.046:
    #    hasMatch = (trigMatches[13] > 0)
    #elif evt_frac < 0.83:
    #    hasMatch = (trigMatches[14] > 0)
    #else:
    #    hasMatch = (trigMatches[15] > 0)
    return hasMatch

muon_selection_mc_reqs = OrderedDict([['mupt1',
                                       [mu_pt,['ell1'],0]],
                                      ['mueta1',
                                       [mu_eta,['ell1'],0]],
                                      ['mupt2',
                                       [mu_pt,['ell2'],0]],
                                      ['mueta2',
                                       [mu_eta,['ell2'],0]],
                                      ['muid1',
                                       [mu_id,['m1IDHZG2012'],0]],
                                      ['muid2',
                                       [mu_id,['m2IDHZG2012'],0]],
                                      ['muiso1',
                                       [mu_iso,['ell1','m1PFChargedIso','m1PFNeutralIso',
                                                'm1PFPhotonIso','m1EffectiveArea2012','m1RhoHZG2012'],0]],
                                      ['muiso2',
                                       [mu_iso,['ell2','m2PFChargedIso','m2PFNeutralIso',
                                                'm2PFPhotonIso','m2EffectiveArea2012','m2RhoHZG2012'],0]],
                                      ['trigger_match',
                                       [mu_trg_match_mc,['muTrgFixed','eventFraction'],0]]])

muon_selection_data_reqs = OrderedDict([['mupt1',
                                         [mu_pt,['ell1'],0]],
                                        ['mueta1',
                                         [mu_eta,['ell1'],0]],
                                        ['mupt2',
                                         [mu_pt,['ell2'],0]],
                                        ['mueta2',
                                         [mu_eta,['ell2'],0]],
                                        ['muid1',
                                         [mu_id,['m1IDHZG2012'],0]],
                                        ['muid2',
                                         [mu_id,['m2IDHZG2012'],0]],
                                        ['muiso1',
                                         [mu_iso,['ell1','m1PFChargedIso','m1PFNeutralIso',
                                                  'm1PFPhotonIso','m1EffectiveArea2012','m1RhoHZG2012'],0]],
                                        ['muiso2',
                                         [mu_iso,['ell2','m2PFChargedIso','m2PFNeutralIso',
                                                  'm2PFPhotonIso','m2EffectiveArea2012','m2RhoHZG2012'],0]],
                                        ['trigger_match',
                                         [mu_trg_match_data,['metEt','run'],0]]])


def ele_pt(thept):    
    return thept > 10

def ele_eta(eleEta):
    return ((fabs(eleEta) < 1.4442 or (fabs(eleEta) > 1.566 and fabs(eleEta) < 2.5)))

#WP85 from 2011 simple cutbased
eidcuts = {'misshits':[0,0],
           'dist':[0.02,0.02],
           'dcot':[0.02,0.02],
           'sihih':[0.01,0.031],
           'dphi':[0.039,0.028],
           'deta':[0.005,0.007]}

def ele_ID(eleSCEta,eleMissHits,eleDist,eleDcot,eleSihih,eleDphivtx,eleDetavtx):
    etabin = int(fabs(eleSCEta) > 1.5)    

    result1 = (eleMissHits == eidcuts['misshits'][etabin] and
               ( fabs(eleDist) > eidcuts['dist'][etabin] or
                 fabs(eleDcot) > eidcuts['dcot'][etabin] ) and
               eleSihih < eidcuts['sihih'][etabin] and
               fabs(eleDphivtx) < eidcuts['dphi'][etabin] and
               fabs(eleDetavtx) < eidcuts['deta'][etabin])
    
    return (result1)

eisocuts = [0.053,0.042]
def ele_iso(elePt,scEta,trk,ecal,hcal,rho):

    iso = -1
    eta1 = -1
    if(fabs(scEta) < 1.4442):
        eta1 = 0
        iso = (trk + max(ecal-1.0,0.0) + hcal - rho*pi*0.3*0.3)/elePt        
    else:
        eta1 = 1
        iso = (trk + ecal + hcal - rho*pi*0.3*0.3)/elePt

    result1 = iso < eisocuts[eta1]    
    return (result1)

def ele_trg_match_data(trigMatch,run):
    hasMatch = False
    if run >= 160404 and run <= 167913:
        hasMatch = (trigMatch[18] > 0)
    elif run >= 170249 and run <= 180252:
        hasMatch = (trigMatch[21] > 0)
    return hasMatch

def ele_trg_match_mc(trigMatch,evt_frac):
    hasMatch = False
    #if evt_frac < 0.236:
    #    hasMatch = (trigMatch[18] > 0)
    #else:
    hasMatch = (trigMatch[21] > 0)
    return hasMatch

def d0_pv(d0):
    return abs(d0) < 0.02

def dz_pv(dz):
    return abs(dz) < 0.05

electron_selection_data_reqs = OrderedDict([['ept',
                                             [ele_pt,['eleCorPt'],0]],
                                            ['eeta',
                                             [ele_eta,['eleSCEta'],0]],
                                            ['eID',
                                             [ele_ID,['eleSCEta','eleConvMissinghit','eleConvDist',
                                                      'eleConvDcot','eleSigmaIEtaIEta',
                                                      'eledPhiAtVtx','eledEtaAtVtx'],0]],
                                            ['ed0',
                                             [d0_pv,['elePVD0'],0]],
                                            ['edz',
                                             [dz_pv,['elePVDz'],0]],
                                            ['eiso',
                                             [ele_iso,['eleCorPt','eleSCEta','eleIsoTrkDR03',
                                                       'eleIsoEcalDR03','eleIsoHcalSolidDR03','rho'],0]],
                                            ['trigger_match',
                                             [ele_trg_match_data,['eleTrgFixed','run'],0]]])

electron_selection_mc_reqs = OrderedDict([['ept',
                                           [ele_pt,['eleCorPt'],0]],
                                          ['eeta',
                                           [ele_eta,['eleSCEta'],0]],
                                          ['eID',
                                           [ele_ID,['eleSCEta','eleConvMissinghit','eleConvDist',
                                                    'eleConvDcot','eleSigmaIEtaIEta',
                                                    'eledPhiAtVtx','eledEtaAtVtx'],0]],
                                          ['ed0',
                                           [d0_pv,['elePVD0'],0]],
                                          ['edz',
                                           [dz_pv,['elePVDz'],0]],
                                          ['eiso',
                                           [ele_iso,['eleCorPt','eleSCEta','eleIsoTrkDR03',
                                                     'eleIsoEcalDR03','eleIsoHcalSolidDR03','rho'],0]],
                                          ['trigger_match',
                                           [ele_trg_match_mc,['eleTrgFixed','eventFraction'],0]]])

def z_oneleg20(ell1, ell2):
    return (ell1.Pt() > 20 or ell2.Pt() > 20)

def z_ss(z_ss):
    return z_ss == 0

def z_mass(z):
    return z.M() > 50

mumu_selection_reqs = OrderedDict([['z_mass',
                                    [z_mass,['Z'],0]],
                                   ['z_ss',
                                    [z_ss,['m1_m2_SS'],0]],
                                   ['z_oneleg20',
                                    [z_oneleg20,['ell1','ell2'],0]]
                                   ])
ee_selection_reqs = OrderedDict([['z_mass',
                                  [z_mass,['Z'],0]],
                                 ['z_ss',
                                  [z_ss,['e1_e2_SS'],0]],
                                 ['z_oneleg20',
                                  [z_oneleg20,['ell1','ell2'],0]]
                                 ])

#photon kinematic selection
def photon_et(pho,Zg):    
    return (pho.Pt() > 15 and pho.Pt()/Zg.M() > 15.0/110.0)

def photon_eta(sc_eta):
    return (abs(sc_eta) < 1.4442 or (abs(sc_eta) > 1.566  and abs(sc_eta) < 2.5 ))

photon_kine_reqs = OrderedDict([['phoet',[photon_et,['gam','Zg'],0]],
                                ['phoeta',[photon_eta,['gSCEta'],0]]])

def photon_CBID_MEDIUM(bit):
    return bit == 1.0

def photon_HoE(hoe):
    return (hoe < 0.05)

def photon_pixveto(hasPixelSeed):
    return (hasPixelSeed == 0)

trk_iso_etslope  = [0.0010,0.0010]
trk_iso_rhoslope = [0.0167,0.0320]
def photon_trkiso(trkIso,phEt,sc_eta,rho):
    eta_bin = int(fabs(sc_eta) > 1.566)
    return (trkIso - trk_iso_etslope[eta_bin]*phEt - trk_iso_rhoslope[eta_bin]*rho < 2.0)

ecal_iso_etslope  = [0.0060,0.0060]
ecal_iso_rhoslope = [0.1830,0.0900]
def photon_ecaliso(ecalIso,phEt,sc_eta,rho):
    eta_bin = int(fabs(sc_eta) > 1.566)
    return (ecalIso - ecal_iso_etslope[eta_bin]*phEt - ecal_iso_rhoslope[eta_bin]*rho < 4.2)

hcal_iso_etslope  = [0.0025,0.0025]
hcal_iso_rhoslope = [0.0620,0.1800]
def photon_hcaliso(hcalIso,phEt,sc_eta,rho):
    eta_bin = int(fabs(sc_eta) > 1.566)
    return (hcalIso - hcal_iso_etslope[eta_bin]*phEt - hcal_iso_rhoslope[eta_bin]*rho < 2.2)

def photon_sihih(sihih,sipip,sc_eta):    
    return (((fabs(sc_eta) < 1.4442 and sihih < 0.011 and sihih > 0.001 and sipip > 0.001) or
            (fabs(sc_eta) < 2.5 and fabs(sc_eta) > 1.566 and sihih < 0.030)))

def photon_like_jet(phEt,sc_eta,rho,hasPixelSeed,hoe,sihih,sipip,
                    trkIso,ecalIso,hcalIso):
    eta_bin = int(fabs(sc_eta) > 1.566)
    result = photon_HoE(hoe) #start from common H-over-E cut

    #plj sigma ieta ieta cut
    result *= ( ( not eta_bin and sihih < 0.014) or
                ( eta_bin and sihih < 0.035) )

    #plj isolation selection
    result *= (trkIso  < min(5*(3.5 + trk_iso_ptslope[eta_bin]*phEt  + trk_iso_rhoslope[eta_bin]*rho), phEt*0.2))
    result *= (ecalIso < min(5*(3.5 + ecal_iso_ptslope[eta_bin]*phEt + ecal_iso_rhoslope[eta_bin]*rho),phEt*0.2))
    result *= (hcalIso < min(5*(3.5 + hcal_iso_ptslope[eta_bin]*phEt + hcal_iso_rhoslope[eta_bin]*rho),phEt*0.2))

    #plj photon veto
    return (result and not( photon_sihih(sihih,sipip,sc_eta)        and
                            photon_trkiso(trkIso,phEt,sc_eta,rho)   and
                            photon_ecaliso(ecalIso,phEt,sc_eta,rho) and
                            photon_hcaliso(ecalIso,phEt,sc_eta,rho)     )) 
                           
    
    

charged_leptons = [11,13,15]
weak_bosons = [23,24]
quarks = [1,2,3,4,5,6]
photon = 22
gluon = 21
proton = 2212
def isnotgenmatchifsr(motherPID,gMotherPID):
    if( abs(motherPID) < 22 ): return False
    if( abs(motherPID) == proton or abs(gMotherPID) == proton ): return False
    if( abs(motherPID) in charged_leptons and abs(gMotherPID) in charged_leptons ): return False #secondary radiation from lepton
    if( abs(motherPID) in charged_leptons and abs(gMotherPID) in weak_bosons ): return False #primary radiation from lepton
    if( motherPID == photon and abs(gMotherPID) in charged_leptons ): return False #photon had some SMC interaction?
    if( motherPID == gluon or abs(motherPID) in quarks ): return False #ISR photon or gluon fragmentation
    return True

#runs over MC particles
def event_has_ifsr(PID, gMotherPID):
    if( abs(PID) == photon and abs(gMotherPID) == 23 ): return True
    if( abs(PID) == photon and (abs(gMotherPID) in quarks or abs(gMotherPID) == gluon) ): return True
    return False

photon_not_ifsr = OrderedDict([['ifsr_veto',[isnotgenmatchifsr,['phoGenMomPID','phoGenGMomPID'],0]]])
event_has_ifsr = OrderedDict([['has_ifsr',[event_has_ifsr,['mcPID','mcGMomPID'],0]]])


photon_idiso_reqs = OrderedDict([['phoCBID_MEDIUM',[photon_CBID_MEDIUM,['gCBID_MEDIUM'],0]]])
#                                         ['phohasPixelSeed',[photon_pixveto,['phohasPixelSeed'],0]],
#                                         ['phoTrkIso',[photon_trkiso,['phoTrkIsoHollowDR04','phoEt','phoSCEta','rho'],0]],
#                                         ['phoEcalIso',[photon_ecaliso,['phoEcalIsoDR04','phoEt','phoSCEta','rho'],0]],
#                                         ['phoHcalIso',[photon_hcaliso,['phoHcalIsoDR04','phoEt','phoSCEta','rho'],0]]])

photon_like_jet_idiso = OrderedDict([['photonLikeJet',[photon_like_jet,['phoCorEt','phoSCEta','rho','phohasPixelSeed',
                                                                        'phoHoverE','phoSigmaIEtaIEta','phoSigmaIPhiIPhi',
                                                                        'phoTrkIsoHollowDR04','phoEcalIsoDR04',
                                                                        'phoHcalIsoDR04']]]])

photon_sihih_req = OrderedDict([['phosihih',[photon_sihih,['phoSigmaIEtaIEta','phoSigmaIPhiIPhi','phoSCEta'],0]]])
photon_mc_sihih_req = OrderedDict([['phosihih',[photon_sihih,['phoCorSigmaIEtaIEta','phoSigmaIPhiIPhi','phoSCEta'],0]]])

photon_idiso_reqs = OrderedDict( photon_idiso_reqs.items() ) #+
#photon_sihih_req.items() )
photon_mc_idiso_reqs = OrderedDict( photon_idiso_reqs.items() )# +
#                                    photon_mc_sihih_req.items() )

def ell_gamma_dr(ell1,ell2,pho):
    mindr = min(pho.DeltaR(ell1),
                pho.DeltaR(ell2))
    return mindr > 0.4


#trigger and eventsel
muon_triggers_data = CutflowDecision(mumu_good_event_data_reqs)
muon_triggers_mc = CutflowDecision(mumu_good_event_mc_reqs)
electron_triggers_data = CutflowDecision(ee_good_event_data_reqs)
electron_triggers_mc = CutflowDecision(ee_good_event_mc_reqs)
event_has_ifsr = CutflowDecision(event_has_ifsr)

#single object IDs
mu_cuts_data = CutflowDecision(muon_selection_data_reqs)
mu_cuts_mc   = CutflowDecision(muon_selection_mc_reqs)
e_cuts_data  = CutflowDecision(electron_selection_data_reqs)
e_cuts_mc    = CutflowDecision(electron_selection_mc_reqs)
photon_cuts_data  = CutflowDecision(photon_kine_reqs.items() +
                                    photon_idiso_reqs.items())
photon_cuts_plj   = CutflowDecision(photon_kine_reqs.items() +
                                    photon_like_jet_idiso.items())
photon_cuts_mc  = CutflowDecision(photon_kine_reqs.items() +
                                  photon_mc_idiso_reqs.items())
photon_cuts_noifsr = CutflowDecision(photon_kine_reqs.items() +
                                     photon_mc_idiso_reqs.items() +
                                     photon_not_ifsr.items())

#z candidate cuts
mumu_cuts = CutflowDecision(mumu_selection_reqs)
ee_cuts = CutflowDecision(ee_selection_reqs)

#zgamma candidate cuts
#llg_cuts = CutflowDecision(zgamma_selection_reqs)

#mumug gamma total selection
#mumug_cuts = CompositeCutflow()
#mumug_cuts.addCutflow('trigger',muon_triggers)
#mumug_cuts.addCutflow('mu1',mu_cuts)
#mumug_cuts.addCutflow('mu2',mu_cuts)
#mumug_cuts.addCutflow('pho',photon_cuts)
#mumug_cuts.addCutflow('z',mumu_cuts)
#mumug_cuts.addCut('mindr',zgamma_selection_reqs['mindr'])

#eeg candidate cuts
#eeg_cuts = CompositeCutflow()
#eeg_cuts.addCutflow('trigger',electron_triggers)
#eeg_cuts.addCutflow('e1',e_cuts)
#eeg_cuts.addCutflow('e2',e_cuts)
#eeg_cuts.addCutflow('pho',photon_cuts)
#eeg_cuts.addCutflow('z',ee_cuts)
#eeg_cuts.addCut('mindr',zgamma_selection_reqs['mindr'])
