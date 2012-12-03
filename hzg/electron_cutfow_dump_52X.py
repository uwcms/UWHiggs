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

def ecal_fiducial(absEta):
    return (absEta < 1.4442 or (absEta > 1.566 and absEta < 2.5))

def trigger_req(event,i):
    return (event.doubleEPass[i] == 1 and event.doubleEPrescale[i] == 1)

def vtx_req(event,i):
    return (event.pvIsValid[i] == 1 and event.pvIsFake[i] == 0)

def e_id(event,i):
    return (event.e1Pt[i] > 10. and event.e2Pt[i] > 10. and
            ecal_fiducial(event.e1AbsEta[i]) and ecal_fiducial(event.e2AbsEta[i]) and
            event.e1CBID_LOOSE[i] == 1.0 and event.e2CBID_LOOSE[i] == 1.0)

def e_iso(event,i):    

    e1Iso = event.e1RelPFIsoDB[i]
    e2Iso = event.e2RelPFIsoDB[i]
    
    """
    e1Iso = ( (event.e1PFChargedIso[i] +
               max(event.e1PFNeutralIso[i] + event.e1PFPhotonIso[i]
                   -event.e1EffectiveArea2012Data[i]*
                   max(event.gRho[i],0.0), #hacky hack of hackiness
                   0.0))/event.e1Pt[i] )
    e2Iso = ( (event.e2PFChargedIso[i] +
               max(event.e2PFNeutralIso[i] + event.e2PFPhotonIso[i]
                   -event.e2EffectiveArea2012Data[i]*
                   max(event.gRho[i],0.0), #hacky hack of hackiness
                   0.0))/event.e2Pt[i] )
    """

    return (e1Iso < 0.4 and e2Iso < 0.4)

def z_id(event,i):
    return ( (event.e1Pt[i] > 20 or event.e2Pt[i] > 20) and
             event.e1_e2_Mass[i] > 50 and
             event.e1_e2_SS[i] == 0)



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
    pfChgIso = event.gPFChargedIsov2[i]
    pfChgEA  = event.gEffectiveAreaCHad[i]
    
    chgIso = max(pfChgIso - rho*pfChgEA,0.0)

    pfNeutIso= event.gPFNeutralIsov2[i]
    pfNeutEA = event.gEffectiveAreaNHad[i]    

    neutIso = max(pfNeutIso - 0.04*event.gPt[i] - rho*pfNeutEA,0.0)

    pfPhoIso = event.gPFPhotonIsov2[i]
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
    pt_over_m = event.gPt[i]/event.Mass[i]
    ascEta = abs(event.gSCEta[i])
    
    return ( pt_over_m > 15.0/110.0 and
             ecal_fiducial(ascEta) and
             event.gCBID_MEDIUM[i] == 1.0)

def pho_fiducial(event,i):
    pt_over_m = event.gPt[i]/event.Mass[i]
    ascEta = abs(event.gSCEta[i])
    
    return ( pt_over_m > 15.0/110.0 and
             ecal_fiducial(ascEta) )

def photon_dr(event,i):
    return min(event.e1_g_DR[i],event.e2_g_DR[i]) > 0.4

def zg_mass_low(event,i):
    return event.Mass[i] > 115.0 

def zg_mass_high(event,i):
    return event.Mass[i] < 180.0


def photon_id_debug(event,i):
    """if( trigger_req(event,i) and
        vtx_req(event,i) and
        mu_id(event,i) and
        mu_iso(event,i) and
        z_id(event,i) and
        good_photon(event,i) and
        photon_dr(event,i) and
        zg_mass_low(event,i) and
        zg_mass_high(event,i) ):
        """
    if (True):
        print "ALLPHOTON :: run %i  evt: %i  pt:%.4f  scEta: %0.6f  hoe: %f" \
              "  sieie: %f  pfCh: %.6f  pfNe: %.6f  pfGa: %.6f  rho: %f  EACh: %.3f   EANeut: %.3f   EAPho: %.3f" \
              %(event.run[i], event.evt[i], event.gPt[i], event.gSCEta[i],
                event.gSingleTowerHadronicDepth1OverEm[i] +
                event.gSingleTowerHadronicDepth2OverEm[i] ,
                event.gSigmaIEtaIEta[i],
                event.gPFChargedIsov2[i],
                event.gPFNeutralIsov2[i],
                event.gPFPhotonIsov2[i],
                event.gRho[i],
                event.gEffectiveAreaCHad[i],
                event.gEffectiveAreaNHad[i],
                event.gEffectiveAreaPho[i])        

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
               e_id, #10 GeV && ID
               e_iso, #ISO
               z_id #Z ID
               ]
counts_ee = [0 for cut in cut_list_ee] + [0]

cut_list_eeg = list(cut_list_ee)
cut_list_eeg += [pho_fiducial,
                 eleVeto,
                 HoverE,
                 sihih,
                 phoIso, #good photon
                 photon_dr, #delta r lepton-photon
                 zg_mass_low,
                 zg_mass_high
                 ]
counts_eeg = [0 for cut in cut_list_eeg] + [0]

process_tuple(eeNtuple,cut_list_ee,counts_ee)

print 'HLT     : %i'%(counts_ee[0])
print 'VTX     : %i'%(counts_ee[1])
print 'Elec ID : %i'%(counts_ee[2])
print 'Elec Iso: %i'%(counts_ee[3])
print 'Z Sel   : %i'%(counts_ee[4])
print
print 'Total ee: %i'%(counts_ee[5])
print

process_tuple(eegNtuple,cut_list_eeg,counts_eeg)#,printer=photon_id_debug)
    
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
