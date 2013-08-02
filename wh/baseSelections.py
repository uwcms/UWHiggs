from FinalStateAnalysis.PlotTools.decorators import memo
from FinalStateAnalysis.Utilities.struct import struct
from electronids import electronIds

@memo
def getVar(name, var):
    return name+var

@memo
def splitEid(label):
    return label.split('_')[-1], label.split('_')[0] 

#OBJECT SELECTION
def muSelection(row, name):
    if getattr( row, getVar(name,'Pt')) < 10:       return False
    if getattr( row, getVar(name,'AbsEta')) > 2.4:  return False
    if not getattr( row, getVar(name,'PixHits')):   return False
    if getattr( row, getVar(name,'JetCSVBtag')) > 0.8: return False
    #if getattr( row, getVar(name,'JetBtag')) > 3.3: return False #was 3.3 
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2: return False
    return True

def eSelection(row, name):
    if getattr( row, getVar(name,'Pt')) < 10:           return False
    if getattr( row, getVar(name,'AbsEta')) > 2.5:      return False
    if getattr( row, getVar(name,'MissingHits')):       return False
    if getattr( row, getVar(name,'HasConversion')):     return False
    if not getattr( row, getVar(name,'ChargeIdTight')): return False
    if getattr( row, getVar(name,'JetCSVBtag')) > 0.8: return False
    #if getattr( row, getVar(name,'JetBtag')) > 3.3:     return False
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
    if row.bjetCSVVetoZHLike: return False
        #if row.bjetCSVVeto:       return False
    if row.eVetoMVAIsoVtx:    return False
    if row.tauVetoPt20Loose3HitsVtx: return False
    return True

def lepton_id_iso(row, name, label): #label in the format eidtype_isotype
    'One function to rule them all'
    LEPTON_ID = False
    isolabel, eidlabel = splitEid(label) #memoizes to be faster!
    if name[0] == 'e':
        LEPTON_ID = electronIds[eidlabel](row, name)
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
    if  row.e1_e2_SS and lepton_id_iso(row, 'e1', 'eid12Medium_h2taucuts') and row.e1MtToMET > 30: 
        return 'wjets'
    elif row.e1_e2_SS and row.e1RelPFIsoDB > 0.3 and row.type1_pfMetEt < 25: #and row.metSignificance < 3: #
        return 'qcd'
    elif lepton_id_iso(row,'e1', 'eid12Medium_h2taucuts') and lepton_id_iso(row,'e2', 'eid12Medium_h2taucuts') \
        and not any([ row.muVetoPt5IsoIdVtx,
                      row.tauVetoPt20Loose3HitsVtx,
                      row.eVetoMVAIsoVtx,
                      ]):
        return 'zee'
    else:
        return None


