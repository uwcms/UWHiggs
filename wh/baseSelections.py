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
def h2taucuts(row, name):
    LEPTON_ID    = getattr(row, getVar(name, 'MVAIDH2TauWP')) \
        if name[0] == 'e' else \
        getattr(row, getVar(name, 'PFIDTight'))
    RelPFIsoDB   = getattr(row, getVar(name, 'RelPFIsoDB'))
    AbsEta       = getattr(row, getVar(name, 'AbsEta'))
    return bool(LEPTON_ID) and bool( RelPFIsoDB < 0.1 or (RelPFIsoDB < 0.15 and AbsEta < 1.479))

def h2taucuts020(row, name):
    LEPTON_ID    = getattr(row, getVar(name, 'MVAIDH2TauWP')) \
        if name[0] == 'e' else \
        getattr(row, getVar(name, 'PFIDTight'))
    RelPFIsoDB   = getattr(row, getVar(name, 'RelPFIsoDB'))
    AbsEta       = getattr(row, getVar(name, 'AbsEta'))
    return bool(LEPTON_ID) and bool( RelPFIsoDB < 0.15 or (RelPFIsoDB < 0.20 and AbsEta < 1.479))

def idiso02(row, name):
    LEPTON_ID    = getattr(row, getVar(name, 'MVAIDH2TauWP')) \
        if name[0] == 'e' else \
        getattr(row, getVar(name, 'PFIDTight'))
    RelPFIsoDB   = getattr(row, getVar(name, 'RelPFIsoDB'))
    AbsEta       = getattr(row, getVar(name, 'AbsEta'))
    return bool(LEPTON_ID) and bool( RelPFIsoDB < 0.20 )

lepton_ids = {
    'h2taucuts'    : h2taucuts,
    'h2taucuts020' : h2taucuts020,
    'idiso02'      : idiso02,
    }


def control_region_ee(row):
    '''Figure out what control region we are in. Shared among two codes, to avoid mismatching copied here'''
    if  row.e1_e2_SS and h2taucuts(row, 'e1') and row.e1MtToMET > 30: # and row.e2MtToMET < 30:# and row.e2MtToMET < 30: #and row.metEt > 30: #row.metSignificance > 3:
        return 'wjets'
    elif row.e1_e2_SS and row.e1RelPFIsoDB > 0.3 and row.type1_pfMetEt < 25: #and row.metSignificance < 3: #
        return 'qcd'
    elif h2taucuts(row,'e1') and h2taucuts(row,'e2') \
        and not any([ row.muVetoPt5IsoIdVtx,
                      row.tauVetoPt20Loose3HitsVtx,
                      row.eVetoMVAIsoVtx,
                      ]):
        return 'zee'
    else:
        return None


