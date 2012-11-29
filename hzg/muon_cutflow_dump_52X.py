#!/usr/bin/env python

import ROOT
from ROOT import TFile, TTree, gDirectory, TH1F

import sys

def makeIsoString(prefix):
    result  = '(%sPFChargedIso + '%(prefix)
    result += 'max(%sPFNeutralIso + %sPFPhotonIso - '
    '%sEffectiveArea2012*max(%sRhoHZG2012,0.0), 0.0))'(prefix,prefix,prefix,prefix)
    result += '/%sPt < 0.12'          
    return result

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
    return (event.doubleMuPass[i] == 1 and event.doubleMuPrescale[i] == 1)

def vtx_req(event,i):
    return (event.pvIsValid[i] == 1 and event.pvIsFake[i] == 0)

def mu_id(event,i):
    return (event.m1Pt[i] > 10 and event.m2Pt[i] > 10 and
            event.m1AbsEta[i] < 2.4 and event.m2AbsEta[i] < 2.4 and
            event.m1IDHZG2012[i] == 1 and event.m2IDHZG2012[i] == 1)

def mu_iso(event,i):
    m1Iso = ( (event.m1PFChargedIso[i] +
               max(event.m1PFNeutralIso[i] + event.m1PFPhotonIso[i]
                   -event.m1EffectiveArea2012[i]*max(event.m1RhoHZG2012[i],0.0),
                   0.0))/event.m1Pt[i] )
    m2Iso = ( (event.m2PFChargedIso[i] +
               max(event.m2PFNeutralIso[i] + event.m2PFPhotonIso[i]
                   -event.m2EffectiveArea2012[i]*max(event.m2RhoHZG2012[i],0.0),
                   0.0))/event.m2Pt[i] )

    return (m1Iso < 0.12 and m2Iso < 0.12)

def z_id(event,i):
    return ( (event.m1Pt[i] > 20 or event.m2Pt[i] > 20) and
             event.m1_m2_Mass[i] > 50 )

def good_photon(event,i):
    pt_over_m = event.gPt[i]/event.Mass[i]
    scEta = event.gSCEta[i]
    
    return ( pt_over_m > 15.0/110.0 and
             (abs(scEta) < 1.4442 or (abs(scEta) > 1.560 and abs(scEta) < 2.5)) and
             event.gCBID_MEDIUM[i] == 1 )

def photon_dr(event,i):
    return min(event.m1_g_DR[i],event.m2_g_DR[i]) > 0.4

def zg_mass_low(event,i):
    return event.Mass[i] > 115.0 

def zg_mass_high(event,i):
    return event.Mass[i] < 180.0

cut_list_mm = [trigger_req, #HLT
               vtx_req, #PV selection
               mu_id, #10 GeV && ID
               mu_iso, #ISO
               z_id #Z ID
               ]
counts_mm = [0 for cut in cut_list_mm] + [0]

cut_list_mmg = list(cut_list_mm)
cut_list_mmg += [good_photon, #good photon
                 photon_dr, #delta r lepton-photon
                 zg_mass_low,
                 zg_mass_high
                 ]
counts_mmg = [0 for cut in cut_list_mmg] + [0]

for event in mmNtuple:
    one_passes = False
    counts_evt = [0 for cut in cut_list_mm]
    
    for i in range(event.N_PATFinalState):
        
        cut_bits = [cut(event,i) for cut in cut_list_mm]
        one_passes = one_passes or (cut_bits.count(True) == len(cut_list_mm))

        passed_last = True
        kbit = 0
        
        while passed_last and kbit < len(cut_bits):
            counts_evt[kbit] += 1*cut_bits[kbit]            
            passed_last = cut_bits[kbit]
            kbit += 1

    for i,count in enumerate(counts_evt):
        counts_mm[i] += 1*(count > 0)
        
    counts_mm[len(cut_list_mm)] += int(one_passes)

print 'HLT     : %i'%(counts_mm[0])
print 'VTX     : %i'%(counts_mm[1])
print 'Muon ID : %i'%(counts_mm[2])
print 'Muon Iso: %i'%(counts_mm[3])
print 'Z Sel   : %i'%(counts_mm[4])
print 'Total MM: %i'%(counts_mm[5])
print


for event in mmgNtuple:
    one_passes = False
    counts_evt = [0 for cut in cut_list_mmg]
    
    for i in range(event.N_PATFinalState):
        
        cut_bits = [cut(event,i) for cut in cut_list_mmg]
        one_passes = one_passes or (cut_bits.count(True) == len(cut_list_mmg))

        passed_last = True
        kbit = 0
        
        while passed_last and kbit < len(cut_bits):
            counts_evt[kbit] += 1*cut_bits[kbit]            
            passed_last = cut_bits[kbit]
            kbit += 1

    for i,count in enumerate(counts_evt):
        counts_mmg[i] += 1*(count > 0)
        
    counts_mmg[len(cut_list_mmg)] += int(one_passes)

print ">=1 Good Photon : %i"%(counts_mmg[len(cut_list_mm)])
print "DR(l,g) > 0.4   : %i"%(counts_mmg[len(cut_list_mm)+1])
print "ZG Mass > 115   : %i"%(counts_mmg[len(cut_list_mm)+2])
print "ZG Mass < 180   : %i"%(counts_mmg[len(cut_list_mm)+3])
print "Total mmg       : %i"%(counts_mmg[len(cut_list_mmg)]) 

file.Close()
