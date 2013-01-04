#!/usr/bin/env python

import ROOT
from ROOT import TFile, TTree, gDirectory, TH1F

import sys,cPickle

if len(sys.argv) != 2:
    sys.exit("this program accepts one argument (the file name!)")

print "Opening %s"%(sys.argv[1])

pwd = gDirectory.GetPath()
file = TFile.Open(sys.argv[1])
gDirectory.cd(pwd)

eventCount = file.Get("mm").Get("eventCount")
mmNtuple  = file.Get("mm").Get("final").Get("Ntuple")
mmgNtuple = file.Get("mmg").Get("final").Get("Ntuple")

print "Initial: %i"%(eventCount.GetEntries())

def trigger_req(event,i):
    return (event.mu17mu8Pass[i] == 1 and event.mu17mu8Prescale[i] == 1)

def vtx_req(event,i):
    return (event.pvIsValid[i] == 1 and event.pvIsFake[i] == 0)

def mu_id(event,i):
    return (event.m1Pt[i] > 10. and event.m2Pt[i] > 10. and
            event.m1AbsEta[i] < 2.4 and event.m2AbsEta[i] < 2.4 and
            event.m1IDHZG2012[i] == 1.0 and event.m2IDHZG2012[i] == 1.0)

def mu_iso(event,i):    
    
    
    m1Iso = ( (event.m1PFChargedIso[i] +
               max(event.m1PFNeutralIso[i] + event.m1PFPhotonIso[i]
                   -0.5*event.m1PFPUChargedIso[i],                   
                   0.0))/event.m1Pt[i] )
    m2Iso = ( (event.m2PFChargedIso[i] +
               max(event.m2PFNeutralIso[i] + event.m2PFPhotonIso[i]
                   -0.5*event.m2PFPUChargedIso[i],
                   0.0))/event.m2Pt[i] )

    return (m1Iso < 0.12 and m2Iso < 0.12)

def z_id(event,i):
    return ( (event.m1Pt[i] > 20 or event.m2Pt[i] > 20) and
             event.m1_m2_Mass[i] > 50 and
             event.m1_m2_SS[i] == 0)



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
    pt_over_m = event.gPt[i]/event.Mass[i]
    ascEta = abs(event.gSCEta[i])
    
    return ( event.gPt[i] > 15.0 and
             (ascEta < 1.4442 or (ascEta > 1.566 and ascEta < 2.5)) and
             pt_over_m > 15.0/110.0 and             
             event.gCBID_MEDIUM[i] == 1.0 
             )

def pho_fiducial(event,i):
    pt_over_m = event.gPt[i]/event.Mass[i]
    ascEta = abs(event.gSCEta[i])
    
    return ( pt_over_m > 15.0/110.0 and
             (ascEta < 1.4442 or (ascEta > 1.566 and ascEta < 2.5)) )

def photon_dr(event,i):
    return min(event.m1_g_DR[i],event.m2_g_DR[i]) > 0.4

def zg_mass_low(event,i):
    return event.Mass[i] > 115.0 

def zg_mass_high(event,i):
    return event.Mass[i] < 180.0


def photon_id_debug(event,i):
    if( trigger_req(event,i) and
        vtx_req(event,i) and
        mu_id(event,i) and
        mu_iso(event,i) and
        z_id(event,i) and
        good_photon(event,i) and
        photon_dr(event,i) and
        zg_mass_low(event,i) and
        zg_mass_high(event,i)
        ):        
    #if (True):
        print "PHOTON :: entry: %i run %i  evt: %i  pt:%.4f  scEta: %0.6f  eVeto: %i  hoe: %f" \
              "  sieie: %f  pfCh: %.6f  pfNe: %.6f  pfGa: %.6f  rho: %f  EACh: %.3f   EANeut: %.3f   EAPho: %.3f" \
              %(i, event.run[i], event.evt[i], event.gPt[i], event.gSCEta[i],
                event.gConvSafeElectronVeto[i],
                event.gSingleTowerHadronicDepth1OverEm[i] +
                event.gSingleTowerHadronicDepth2OverEm[i] ,
                event.gSigmaIEtaIEta[i],
                event.gPFChargedIso[i],
                event.gPFNeutralIso[i],
                event.gPFPhotonIso[i],
                event.gRho[i],
                event.gEffectiveAreaCHad[i],
                event.gEffectiveAreaNHad[i],
                event.gEffectiveAreaPho[i])        

mm_events = {}
def mm_dump(event):
    key = (i,event.run[i],event.lumi[i],event.evt[i]) #index on entry,run,lumi,event
    entry = {}

mmg_events = {}
def mmg_dump(event):
    key = (event.run[0],event.lumi[0],event.evt[0])
    entries = []

    for i in range(event.N_PATFinalState):        
        entry = {
            #trigger, vertex
            'trigger':event.mu17mu8Pass[i],'ngoodvtx':event.nvtx[i],
            #rhos
            'muonRho':event.m1RhoHZG2012[i],'photonRho':event.gRho[i],
            # masses charge delta-r
            'zgMass':event.Mass[i],'zMass':event.m1_m2_Mass[i],
            'sameSign':event.m1_m2_SS[i],
            'm1_g_DR':event.m1_g_DR[i],'m2_g_DR':event.m2_g_DR[i],
            #muon1
            'm1Pt':event.m1Pt[i],'m1Eta':event.m1Eta[i],'m1Phi':event.m1Phi[i],
            'm1ID':event.m1IDHZG2012[i],
            'm1PFChargedIso':event.m1PFChargedIso[i],'m1PFNeutralIso':event.m1PFNeutralIso[i],
            'm1PFPhotonIso':event.m1PFPhotonIso[i],            
            #muon 2
            'm2Pt':event.m2Pt[i],'m2Eta':event.m2Eta[i],'m2Phi':event.m2Phi[i],
            'm2ID':event.m2IDHZG2012[i],
            'm2PFChargedIso':event.m2PFChargedIso[i],'m2PFNeutralIso':event.m2PFNeutralIso[i],
            'm2PFPhotonIso':event.m2PFPhotonIso[i],
            #photon
            'gPt':event.gPt[i],'gEta':event.gEta[i],'gSCEta':event.gSCEta[i],
            'gPhi':event.gPhi[i],'gConvSafeElectronVeto':event.gConvSafeElectronVeto[i],
            'gSingleTowerHadronicOverEm':event.gSingleTowerHadronicOverEm[i],
            'gSigmaIEtaIEta':event.gSigmaIEtaIEta[i],'gPFChargedIso':event.gPFChargedIso[i],
            'gPFNeutralIso':event.gPFNeutralIso[i],
            'gEffectiveAreaCHad':event.gEffectiveAreaCHad[i],
            'gEffectiveAreaNHad':event.gEffectiveAreaNHad[i],
            'gEffectiveAreaPho':event.gEffectiveAreaPho[i]            
            }
        entries.append(entry)
    mmg_events[key] = entries
        

def process_tuple(tuple,cut_list,counts,dumper=None):    
    for event in tuple:
        one_passes = False
        counts_evt = [0 for cut in cut_list]

        
    
        for i in range(event.N_PATFinalState):            

            if dumper is not None:
                dumper(event,i)
            
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

cut_list_mm = [trigger_req, #HLT
               vtx_req, #PV selection
               mu_id, #10 GeV && ID
               mu_iso, #ISO
               z_id #Z ID
               ]
counts_mm = [0 for cut in cut_list_mm] + [0]

cut_list_mmg = list(cut_list_mm)
cut_list_mmg += [good_photon,
                 eleVeto,
                 HoverE,
                 sihih,
                 phoIso, #good photon
                 photon_dr, #delta r lepton-photon
                 zg_mass_low,
                 zg_mass_high
                 ]
counts_mmg = [0 for cut in cut_list_mmg] + [0]

process_tuple(mmNtuple,cut_list_mm,counts_mm)

print 'HLT     : %i'%(counts_mm[0])
print 'VTX     : %i'%(counts_mm[1])
print 'Muon ID : %i'%(counts_mm[2])
print 'Muon Iso: %i'%(counts_mm[3])
print 'Z Sel   : %i'%(counts_mm[4])
print
print 'Total MM: %i'%(counts_mm[5])
print

process_tuple(mmgNtuple,cut_list_mmg,counts_mmg)#,printer=photon_id_debug)
    
print "Fiducial Cuts   : %i"%(counts_mmg[len(cut_list_mm)])
print "Electron Veto   : %i"%(counts_mmg[len(cut_list_mm)+1])
print "ST HoE          : %i"%(counts_mmg[len(cut_list_mm)+2])
print "SIHIH           : %i"%(counts_mmg[len(cut_list_mm)+3])
print "PF Iso          : %i"%(counts_mmg[len(cut_list_mm)+4])
print "DR(l,g) > 0.4    : %i"%(counts_mmg[len(cut_list_mm)+5])
print "ZG Mass > 115    : %i"%(counts_mmg[len(cut_list_mm)+6])
print "ZG Mass < 180    : %i"%(counts_mmg[len(cut_list_mm)+7])
print
print "Total mmg        : %i"%(counts_mmg[len(cut_list_mmg)])
print

cPickle.dump( mmg_events, open( "mmg_event_dump.pkl", "wb" ) )

file.Close()
