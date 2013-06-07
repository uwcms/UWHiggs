from FinalStateAnalysis.PlotTools.decorators import memo
from FinalStateAnalysis.Utilities.struct import struct

@memo
def getVar(name, var):
    return name+var

#OBJECT SELECTION
def muSelection(row, name):
    if getattr( row, getVar(name,'Pt')) < 10:       return False
    if getattr( row, getVar(name,'AbsEta')) > 2.4:  return False
    if not getattr( row, getVar(name,'PixHits')):   return False
    if getattr( row, getVar(name,'JetBtag')) > 3.3: return False
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2: return False
    return True

def eSelection(row, name):
    if getattr( row, getVar(name,'Pt')) < 10:           return False
    if getattr( row, getVar(name,'AbsEta')) > 2.5:      return False
    if getattr( row, getVar(name,'MissingHits')):       return False
    if getattr( row, getVar(name,'HasConversion')):     return False
    if not getattr( row, getVar(name,'ChargeIdTight')): return False
    if getattr( row, getVar(name,'JetBtag')) > 3.3:     return False
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2:     return False
    return True
    
def tauSelection(row, name):
    if getattr( row, getVar(name,'Pt')) < 20:          return False
    if getattr( row, getVar(name,'AbsEta')) > 2.3:     return False
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2:    return False
    return True


#VETOS
def vetos(row):
    if row.muVetoPt5IsoIdVtx: return False
        #if row.bjetCSVVeto:       return False
    if row.eVetoMVAIsoVtx:    return False
    if row.tauVetoPt20Loose3HitsVtx: return False
    return True

#LEPTON ID-ISO
def summer_2013_eid(row, name):
    mva_output = getattr(row, getVar(name, 'MVATrigNoIP'))
    pT    = getattr(row, getVar(name, 'Pt'))
    abseta= getattr(row, getVar(name, 'AbsEta'))
    if pT < 20    and abseta < 0.8:
        return ( mva_output > -0.5375 )
    elif pT < 20  and 0.8 < abseta < 1.479:
        return ( mva_output > -0.375 )
    elif pT < 20  and abseta > 1.479:
        return ( mva_output > -0.025 )
    elif pT > 20  and abseta < 0.8:
        return ( mva_output > 0.325 )
    elif pT > 20  and 0.8 < abseta < 1.479:
        return ( mva_output > 0.775 )
    elif pT > 20  and abseta > 1.479:
        return ( mva_output > 0.775 )

def summer_2013_eid_tight(row, name):
    mva_output = getattr(row, getVar(name, 'MVATrigNoIP'))
    pT    = getattr(row, getVar(name, 'Pt'))
    abseta= getattr(row, getVar(name, 'AbsEta'))
    if pT < 20 and abseta < 0.8:
        return ( mva_output > -0.35 )
    elif pT < 20 and 0.8 < abseta < 1.479:
        return ( mva_output > 0.0 )
    elif pT < 20 and abseta > 1.479:
        return ( mva_output > 0.025 )
    elif pT > 20 and abseta < 0.8:
        return ( mva_output > 0.7 )
    elif pT > 20 and 0.8 < abseta < 1.479:
        return ( mva_output > 0.9 )
    elif pT > 20 and abseta > 1.479:
        return ( mva_output > 0.8375 )

def lepton_id_iso(row, name, label):
    'One function to rule them all'
    LEPTON_ID = False
    isolabel  = label.split('eid13')[-1]
    if name[0] == 'e':
        if not label.startswith('eid13'):
            LEPTON_ID = bool(getattr(row, getVar(name, 'MVAIDH2TauWP')))
        elif label.startswith('eid13Loose'):
            LEPTON_ID = summer_2013_eid(row, name)
        elif label.startswith('eid13Tight'):
            LEPTON_ID = summer_2013_eid_tight(row, name)
    else:
        LEPTON_ID = getattr(row, getVar(name, 'PFIDTight'))
    if not LEPTON_ID:
        return False
    RelPFIsoDB   = getattr(row, getVar(name, 'RelPFIsoDB'))
    AbsEta       = getattr(row, getVar(name, 'AbsEta'))
    if isolabel == 'h2taucuts':
        return bool( RelPFIsoDB < 0.1 or (RelPFIsoDB < 0.15 and AbsEta < 1.479))
    if isolabel == 'h2taucuts020':
        return bool( RelPFIsoDB < 0.15 or (RelPFIsoDB < 0.20 and AbsEta < 1.479))
    if isolabel == 'idiso02':
        return bool( RelPFIsoDB < 0.20 )

def control_region_ee(row):
    '''Figure out what control region we are in. Shared among two codes, to avoid mismatching copied here'''
    if  row.e1_e2_SS and lepton_id_iso(row, 'e1', 'h2taucuts') and row.e1MtToMET > 30: # and row.e2MtToMET < 30:# and row.e2MtToMET < 30: #and row.metEt > 30: #row.metSignificance > 3:
        return 'wjets'
    elif row.e1_e2_SS and row.e1RelPFIsoDB > 0.3 and row.type1_pfMetEt < 25: #and row.metSignificance < 3: #
        return 'qcd'
    elif lepton_id_iso(row,'e1', 'h2taucuts') and lepton_id_iso(row,'e2', 'h2taucuts') \
        and not any([ row.muVetoPt5IsoIdVtx,
                      row.tauVetoPt20Loose3HitsVtx,
                      row.eVetoMVAIsoVtx,
                      ]):
        return 'zee'
    else:
        return None


