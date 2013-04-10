#!/usr/bin/env python

import ROOT
from ROOT import TFile, TTree, gDirectory, TH1F

import sys

if len(sys.argv) != 2:
    sys.exit("this program accepts one argument (the file name!)")

print "Opening %s"%(sys.argv[1])

pwd = gDirectory.GetPath()
file = TFile.Open(sys.argv[1])
gDirectory.cd(pwd)

eventCount = file.Get("ee").Get("eventCount")
eeNtuple  = file.Get("ee").Get("final").Get("Ntuple")
eegNtuple = file.Get("eeg").Get("final").Get("Ntuple")

print "Initial: %i"%(eventCount.GetEntries())

def ecal_fiducial(eta):
    absEta = abs(eta)
    return (absEta < 1.4442 or (absEta > 1.566 and absEta < 2.5))

def trigger_req(event,i):
    return (event.doubleETightPass[i] == 1 and event.doubleETightPrescale[i] == 1)

def vtx_req(event,i):
    return (event.nvtx[i] > 0)

# [barrel, endcap], LOOSE WP
dEtaCut  = [0.007,0.009]
dPhiCut  = [0.15,0.10]
sihihCut = [0.010,0.030]
HoECut   = [0.12,0.10]
d0Cut    = [0.02, 0.02]
dZCut    = [0.2, 0.2]
ooemoopCut  = [0.05,0.05]
hasConvCut  = [False,False]
missHitsCut = [1,1]

def e_id_fun(event,i):    
    idxe1 = int(abs(event.e1Eta[i]) >= 1.566)
    idxe2 = int(abs(event.e2Eta[i]) >= 1.566)
    
    return \
     (abs(event.e1deltaEtaSuperClusterTrackAtVtx[i]) < dEtaCut[idxe1] and abs(event.e2deltaEtaSuperClusterTrackAtVtx[i]) < dEtaCut[idxe2] and
      abs(event.e1deltaPhiSuperClusterTrackAtVtx[i]) < dPhiCut[idxe1] and abs(event.e2deltaPhiSuperClusterTrackAtVtx[i]) < dPhiCut[idxe2] and
      event.e1SigmaIEtaIEta[i] < sihihCut[idxe1] and event.e2SigmaIEtaIEta[i] < sihihCut[idxe2] and
      event.e1HadronicOverEM[i] < HoECut[idxe1] and event.e2HadronicOverEM[i] < HoECut[idxe2] and
      abs(event.e1PVDXY[i]) < d0Cut[idxe1] and abs(event.e2PVDXY[i]) < d0Cut[idxe2] and
      abs(event.e1PVDZ[i]) < dZCut[idxe1] and abs(event.e2PVDZ[i]) < dZCut[idxe2] and
      ( abs(1.0 - event.e1eSuperClusterOverP[i])/event.e1ecalEnergy[i] < ooemoopCut[idxe1] and
        abs(1.0 - event.e2eSuperClusterOverP[i])/event.e2ecalEnergy[i] < ooemoopCut[idxe2]    ) and
      event.e1HasMatchedConversion[i] == int(hasConvCut[idxe1]) and event.e2HasMatchedConversion[i] == int(hasConvCut[idxe2]) and
      int(event.e1MissingHits[i]) <= missHitsCut[idxe1] and  int(event.e2MissingHits[i]) <= missHitsCut[idxe2])

def e_id(event,i):
    return (event.e1Pt[i] > 10. and event.e2Pt[i] > 10. and
            ecal_fiducial(event.e1SCEta[i]) and ecal_fiducial(event.e2SCEta[i]) and
            e_id_fun(event,i) and
            event.e1NearMuonVeto[i] == 0.0 and event.e2NearMuonVeto[i] == 0.0)

def e_iso(event,i):    

    e1Iso = ( (event.e1PFChargedIso[i] +
               max(event.e1PFNeutralIso[i] + event.e1PFPhotonIso[i]
                   -event.e1EffectiveArea2012Data[i]*
                   max(event.e1RhoHZG2012[i],0.0), 
                   0.0))/event.e1Pt[i] )
    e2Iso = ( (event.e2PFChargedIso[i] +
               max(event.e2PFNeutralIso[i] + event.e2PFPhotonIso[i]
                   -event.e2EffectiveArea2012Data[i]*
                   max(event.e2RhoHZG2012[i],0.0), 
                   0.0))/event.e2Pt[i] )

    return (e1Iso < 0.4 and e2Iso < 0.4)


#the pieces of the electron MVA id
sihihCutMVA = [0.014,0.035]
HoECutMVA   = [0.15,0.10]
missHitsMVA = [0,0]
def e_mva_preselection(event,i):
    idxe1 = int(abs(event.e1SCEta[i]) >= 1.479)
    idxe2 = int(abs(event.e2SCEta[i]) >= 1.479)
    
    ecalIso1 = event.e1EcalIsoDR03[i]
    ecalIso2 = event.e2EcalIsoDR03[i]

    return ( event.e1Pt[i] > 10 and
             event.e2Pt[i] > 10 and
             event.e1SigmaIEtaIEta[i] < sihihCutMVA[idxe1] and
             event.e2SigmaIEtaIEta[i] < sihihCutMVA[idxe2] and
             event.e1HadronicOverEM[i] < HoECutMVA[idxe1] and
             event.e2HadronicOverEM[i] < HoECutMVA[idxe2] and
             event.e1TrkIsoDR03[i]/event.e1Pt[i] < 0.2 and
             event.e2TrkIsoDR03[i]/event.e2Pt[i] < 0.2 and
             ecalIso1/event.e1Pt[i] < 0.2 and
             ecalIso2/event.e2Pt[i] < 0.2 and
             event.e1HcalIsoDR03[i]/event.e1Pt[i] < 0.2 and
             event.e2HcalIsoDR03[i]/event.e2Pt[i] < 0.2 and
             int(event.e1MissingHits[i]) <= missHitsMVA[idxe1] and
             int(event.e2MissingHits[i]) <= missHitsMVA[idxe2] )

mvacuts = [0.5,0.12,0.6]
def e_mva_idiso(event,i):    
    e1idx = int(abs(event.e1SCEta[i]) > 0.8) + int(abs(event.e1SCEta[i]) > 1.479)
    e2idx = int(abs(event.e2SCEta[i]) > 0.8) + int(abs(event.e2SCEta[i]) > 1.479)
    e1mva = event.e1MVATrig[i]
    e2mva = event.e2MVATrig[i]

    return ( event.e1Pt > 10 and
             event.e2Pt > 10 and
             e1mva > mvacuts[e1idx] and
             e2mva > mvacuts[e2idx] and
             e_iso(event,i)             )

## this is for the ISISO MVA
#dEtaCutMVA  = [0.007,0.009]
#dPhiCutMVA  = [0.15,0.10]
#sihihCutMVA = [0.010,0.030]
#HoECutMVA   = [0.12,0.10]
#def e_mva_preselection(event,i):
#    idxe1 = int(abs(event.e1SCEta[i]) >= 1.479)
#    idxe2 = int(abs(event.e2SCEta[i]) >= 1.479)
#
#    ecalIso1 = event.e1EcalIsoDR03[i]
#    ecalIso2 = event.e2EcalIsoDR03[i]
#
#    if idxe1 == 0: ecalIso1 = max(ecalIso1 - 1.0, 0.0)
#    if idxe2 == 0: ecalIso2 = max(ecalIso2 - 1.0, 0.0)
#    
#    return ( event.e1Pt[i] > 10. and event.e2Pt[i] > 10. and
#             #abs(event.e1SCEta[i]) < 2.5 and
#             #abs(event.e2SCEta[i]) < 2.5 and
#             abs(event.e1deltaEtaSuperClusterTrackAtVtx[i]) < dEtaCutMVA[idxe1] and
#             abs(event.e2deltaEtaSuperClusterTrackAtVtx[i]) < dEtaCutMVA[idxe2] and
#             abs(event.e1deltaPhiSuperClusterTrackAtVtx[i]) < dPhiCutMVA[idxe1] and
#             abs(event.e2deltaPhiSuperClusterTrackAtVtx[i]) < dPhiCutMVA[idxe2] and
#             event.e1SigmaIEtaIEta[i] < sihihCutMVA[idxe1] and
#             event.e2SigmaIEtaIEta[i] < sihihCutMVA[idxe2] and
#             event.e1HadronicOverEM[i] < HoECutMVA[idxe1] and
#             event.e2HadronicOverEM[i] < HoECutMVA[idxe2] and
#             event.e1TrkIsoDR03[i]/event.e1Pt[i] < 0.2 and
#             event.e2TrkIsoDR03[i]/event.e2Pt[i] < 0.2 and
#             ecalIso1/event.e1Pt[i] < 0.2 and
#             ecalIso2/event.e2Pt[i] < 0.2 and
#             event.e1HcalIsoDR03[i]/event.e1Pt[i] < 0.2 and
#             event.e2HcalIsoDR03[i]/event.e2Pt[i] < 0.2 and
#             int(event.e1MissingHits[i]) <= 1 and
#             int(event.e2MissingHits[i]) <= 1 and
#             abs(event.e1PVDZ[i]) < 0.1 and
#             abs(event.e2PVDZ[i]) < 0.1 ) #and
             #event.e1NearMuonVeto[i] == 0.0 and
             #event.e2NearMuonVeto[i] == 0.0) 

## id iso mva
#mvacuts = [-0.82,-0.98] # pt < 20, pt > 20 
#def e_mva_idiso(event,i):
#    e1idx = int(event.e1Pt[i] > 20)
#    e2idx = int(event.e2Pt[i] > 20)
#    e1mva = event.e1MVATrigIDISO[i]
#    e2mva = event.e2MVATrigIDISO[i]
#
#    return ( e1mva > mvacuts[e1idx] and
#             e2mva > mvacuts[e2idx]      )

from ROOT import TLorentzVector

lep1,lep2 = TLorentzVector(),TLorentzVector()
pho = TLorentzVector()

def z_id(event,i):
    lep1.SetPtEtaPhiM(event.e1PtCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e1EtaCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e1PhiCorrReg_Summer12_DR53X_HCP2012[i],
                      0.000511)
    lep2.SetPtEtaPhiM(event.e2PtCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e2EtaCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e2PhiCorrReg_Summer12_DR53X_HCP2012[i],
                      0.000511)
    
    return ( (event.e1Pt[i] > 20 or
              event.e2Pt[i] > 20) and
             (lep1+lep2).M() > 50 and
             event.charge[i] == 0)

def eleVeto(event,i):
    eVeto = event.gConvSafeElectronVeto[i]

    return eVeto == True

def HoverE(event,i):
    singleTowerHoE = event.gSingleTowerHadronicOverEm[i]
    return singleTowerHoE<0.05

def sihih(event,i):
    sihih = event.gSigmaIEtaIEta[i]
    ascEta = abs(event.gSCEta[i])

    if( ascEta < 1.5 ):
        return (sihih < 0.011)        
    else:
        return (sihih < 0.033)
        
def phoIso(event,i):
    result = []
    ascEta = abs(event.gSCEta[i])

    rho = event.gRho[i]
    pfChgIso = event.gPFChargedIso[i]
    pfChgEA  = event.gEffectiveAreaCHad[i]
    
    chgIso = max(pfChgIso - rho*pfChgEA,0.0)

    pfNeutIso= event.gPFNeutralIso[i]
    pfNeutEA = event.gEffectiveAreaNHad[i]    

    neutIso = max(pfNeutIso - 0.04*event.gPt[i] - rho*pfNeutEA,0.0)

    pfPhoIso = event.gPFPhotonIso[i]
    pfPhoEA  = event.gEffectiveAreaPho[i]

    phoIso = max(pfPhoIso - 0.005*event.gPt[i] - rho*pfPhoEA,0.0)

    if( ascEta < 1.5 ):        
        result.append(chgIso < 1.5)
        result.append(neutIso < 1.0)
        result.append(phoIso < 0.7)
    else:       
        result.append(chgIso < 1.2)
        result.append(neutIso < 1.5)
        result.append(phoIso < 1.0)
        
    return (result.count(True)==3)


def good_photon(event,i):
    lep1.SetPtEtaPhiM(event.e1PtCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e1EtaCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e1PhiCorrReg_Summer12_DR53X_HCP2012[i],
                      0.000511)
    lep2.SetPtEtaPhiM(event.e2PtCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e2EtaCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e2PhiCorrReg_Summer12_DR53X_HCP2012[i],
                      0.000511)
    pho.SetPtEtaPhiM(event.gPt[i],
                     event.gEta[i],
                     event.gPhi[i],
                     0.0)
    pt_over_m = event.gPt[i]/(lep1+lep2+pho).M()
    ascEta = abs(event.gSCEta[i])
    
    return ( event.gPt[i] > 15.0 and
             ecal_fiducial(ascEta) and
             pt_over_m > 15.0/110.0 and             
             event.gCBID_MEDIUM[i] == 1.0)

def pho_fiducial(event,i):
    pt_over_m = event.gPt[i]/event.Mass[i]
    ascEta = abs(event.gSCEta[i])
    
    return ( event.gPt[i] > 15.0 and
             pt_over_m > 15.0/110.0 and
             ecal_fiducial(ascEta) )

def photon_dr(event,i):
    return min(event.e1_g_DR[i],event.e2_g_DR[i]) > 0.4

def zg_mass_low(event,i):
    lep1.SetPtEtaPhiM(event.e1PtCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e1EtaCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e1PhiCorrReg_Summer12_DR53X_HCP2012[i],
                      0.000511)
    lep2.SetPtEtaPhiM(event.e2PtCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e2EtaCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e2PhiCorrReg_Summer12_DR53X_HCP2012[i],
                      0.000511)
    pho.SetPtEtaPhiM(event.gPt[i],
                     event.gEta[i],
                     event.gPhi[i],
                     0.0)
    return (lep1+lep2+pho).M() > 115.0 

def zg_mass_high(event,i):
    lep1.SetPtEtaPhiM(event.e1PtCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e1EtaCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e1PhiCorrReg_Summer12_DR53X_HCP2012[i],
                      0.000511)
    lep2.SetPtEtaPhiM(event.e2PtCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e2EtaCorrReg_Summer12_DR53X_HCP2012[i],
                      event.e2PhiCorrReg_Summer12_DR53X_HCP2012[i],
                      0.000511)
    pho.SetPtEtaPhiM(event.gPt[i],
                     event.gEta[i],
                     event.gPhi[i],
                     0.0)
    return (lep1+lep2+pho).M() < 180.0

def electron_id_debug(event,i):
    if( trigger_req(event,i) and
        vtx_req(event,i) and
        e_mva_preselection(event,i) and
        e_mva_idiso(event,i) and
        z_id(event,i) ):
        print "ELECTRON :: run: %i  evt: %i  pt: %.4f  eReg: %.6f  scEta: %.6f hoe: %.4f" \
              "  sieie: %.6f  dEtaIn: %.6f  dPhiIn: %.6f" \
              "  nhitsmiss: %i dZ: %.4f  trkIso: %.4f  ecalIso: %.4f  hcalIso: %.4f" \
              "  pfChHad: %.6f  pfNeut: %.6f  pfPho: %.6f  ea: %.6f  rho: %.6f" \
              "  IDISOmva: %.4f" \
              %(event.run[i], event.evt[i], event.e1Pt[i],
                event.e1ECorrReg_Summer12_DR53X_HCP2012[i],event.e1SCEta[i],
                event.e1HadronicOverEM[i], event.e1SigmaIEtaIEta[i],
                event.e1deltaEtaSuperClusterTrackAtVtx[i],
                event.e1deltaPhiSuperClusterTrackAtVtx[i],
                event.e1MissingHits[i],event.e1PVDZ[i], event.e1TrkIsoDR03[i],
                event.e1EcalIsoDR03[i],event.e1HcalIsoDR03[i],
                event.e1PFChargedIso[i], 
                event.e1PFNeutralIso[i],
                event.e1PFPhotonIso[i],
                event.e1EffectiveArea2012Data[i],
                event.e1RhoHZG2012[i],        
                event.e1MVATrig[i])
        print "ELECTRON :: run: %i  evt: %i  pt: %.4f  eReg: %.6f  scEta: %.6f hoe: %.4f" \
              "  sieie: %.6f  dEtaIn: %.6f  dPhiIn: %.6f" \
              "  nhitsmiss: %i dZ: %.4f  trkIso: %.4f  ecalIso: %.4f  hcalIso: %.4f" \
              "  pfChHad: %.6f  pfNeut: %.6f  pfPho: %.6f  ea: %.6f  rho: %.6f" \
              "  IDISOmva: %.4f" \
              %(event.run[i], event.evt[i], event.e2Pt[i],
                event.e2ECorrReg_Summer12_DR53X_HCP2012[i],event.e2SCEta[i],
                event.e2HadronicOverEM[i], event.e2SigmaIEtaIEta[i],
                event.e2deltaEtaSuperClusterTrackAtVtx[i],
                event.e2deltaPhiSuperClusterTrackAtVtx[i],
                event.e2MissingHits[i],event.e2PVDZ[i], event.e2TrkIsoDR03[i],
                event.e2EcalIsoDR03[i],event.e2HcalIsoDR03[i],
                event.e2PFChargedIso[i],
                event.e2PFNeutralIso[i],
                event.e2PFPhotonIso[i],
                event.e2EffectiveArea2012Data[i],
                event.e2RhoHZG2012[i],
                event.e2MVATrig[i])
        

def photon_id_debug(event,i):
    if( trigger_req(event,i) and
        vtx_req(event,i) and
        e_mva_preselection(event,i) and
        e_mva_idiso(event,i) and
        z_id(event,i) and
        #good_photon(event,i) and
        #photon_dr(event,i) and
        #zg_mass_low(event,i) and
        #zg_mass_high(event,i) and # ):
        event.evt[i] == 17380 ):
        lep1.SetPtEtaPhiM(event.e1PtCorrReg_Summer12_DR53X_HCP2012[i],
                          event.e1EtaCorrReg_Summer12_DR53X_HCP2012[i],
                          event.e1PhiCorrReg_Summer12_DR53X_HCP2012[i],
                          0.000511)
        lep2.SetPtEtaPhiM(event.e2PtCorrReg_Summer12_DR53X_HCP2012[i],
                          event.e2EtaCorrReg_Summer12_DR53X_HCP2012[i],
                          event.e2PhiCorrReg_Summer12_DR53X_HCP2012[i],
                          0.000511)
        pho.SetPtEtaPhiM(event.gPt[i],
                         event.gEta[i],
                         event.gPhi[i],
                         0.0)
        print "PHOTON :: run %i  evt: %i  pt: %.4f  scEta: %0.6f  hoe: %f" \
              "  sieie: %f  pfCh: %.6f  pfNe: %.6f  pfGa: %.6f  rho: %f  EACh: %.3f   EANeut: %.3f   EAPho: %.3f" \
              "  dR1: %.6f  dR2: %.6f  M_Z: %.6f M_Zg: %.6f" \
              %(event.run[i], event.evt[i], event.gPt[i], event.gSCEta[i],
                event.gSingleTowerHadronicDepth1OverEm[i] +
                event.gSingleTowerHadronicDepth2OverEm[i] ,
                event.gSigmaIEtaIEta[i],
                event.gPFChargedIso[i],
                event.gPFNeutralIso[i],
                event.gPFPhotonIso[i],
                event.gRho[i],
                event.gEffectiveAreaCHad[i],
                event.gEffectiveAreaNHad[i],
                event.gEffectiveAreaPho[i],
                event.e1_g_DR[i],
                event.e2_g_DR[i],
                (lep1+lep2).M(),
                (lep1+lep2+pho).M(),
                )        

def process_tuple(tuple,cut_list,counts,printer=None):    
    for event in tuple:
        one_passes = False
        counts_evt = [0 for cut in cut_list]
    
        for i in range(event.N_PATFinalState):

            if printer is not None:
                printer(event,i)
            
            cut_bits = [cut(event,i) for cut in cut_list]
            one_passes = one_passes or (cut_bits.count(True) == len(cut_list))
            
            passed_last = True
            kbit = 0            
            while passed_last and kbit < len(cut_bits):
                counts_evt[kbit] += 1*cut_bits[kbit]            
                passed_last = cut_bits[kbit]
                kbit += 1

        for i,count in enumerate(counts_evt):
            counts[i] += 1*(count > 0)

        counts[len(cut_list)] += int(one_passes)

cut_list_ee = [trigger_req, #HLT
               vtx_req, #PV selection
               e_mva_preselection, #10 GeV && ID
               e_mva_idiso, #ISO
               z_id #Z ID
               ]
counts_ee = [0 for cut in cut_list_ee] + [0]

cut_list_eeg = list(cut_list_ee)
cut_list_eeg += [good_photon,
                 good_photon,
                 good_photon,
                 good_photon,
                 good_photon, #good photon
                 photon_dr, #delta r lepton-photon
                 zg_mass_low,
                 zg_mass_high
                 ]
counts_eeg = [0 for cut in cut_list_eeg] + [0]

process_tuple(eeNtuple,cut_list_ee,counts_ee)#,electron_id_debug)

print 'HLT     : %i'%(counts_ee[0])
print 'VTX     : %i'%(counts_ee[1])
print 'Elec preMVA : %i'%(counts_ee[2])
print 'Elec MVA ID : %i'%(counts_ee[3])
print 'Z Sel   : %i'%(counts_ee[4])
print
print 'Total ee: %i'%(counts_ee[5])
print

process_tuple(eegNtuple,cut_list_eeg,counts_eeg,printer=photon_id_debug)
    
print "Fiducial Cuts   : %i"%(counts_eeg[len(cut_list_ee)])
print "Electron Veto   : %i"%(counts_eeg[len(cut_list_ee)+1])
print "ST HoE          : %i"%(counts_eeg[len(cut_list_ee)+2])
print "SIHIH           : %i"%(counts_eeg[len(cut_list_ee)+3])
print "PF Iso          : %i"%(counts_eeg[len(cut_list_ee)+4])
print "DR(l,g) > 0.4    : %i"%(counts_eeg[len(cut_list_ee)+5])
print "ZG Mass > 115    : %i"%(counts_eeg[len(cut_list_ee)+6])
print "ZG Mass < 180    : %i"%(counts_eeg[len(cut_list_ee)+7])
print
print "Total eeg        : %i"%(counts_eeg[len(cut_list_eeg)])
print

file.Close()
