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
    if getattr( row, getVar(name,'Pt')) < 30:       return False
    if getattr( row, getVar(name,'AbsEta')) > 2.1:  return False
    if not getattr( row, getVar(name,'PixHits')):   return False
    if getattr( row, getVar(name,'JetCSVBtag')) > 0.8: return False
    #if getattr( row, getVar(name,'JetBtag')) > 3.3: return False #was 3.3 
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2: return False
    return True

def eSelection(row, name):
    eAbsEta = getattr( row, getVar(name,'AbsEta'))
    if getattr( row, getVar(name,'Pt')) < 30:           return False #was 20
    if eAbsEta > 2.3:      return False
    if getattr( row, getVar(name,'MissingHits')):       return False
    if getattr( row, getVar(name,'HasConversion')):     return False
    if eAbsEta > 1.4442 and eAbsEta < 1.566: return False
#    if not getattr( row, getVar(name,'ChargeIdTight')): return False
    if not getattr( row, getVar(name,'ChargeIdLoose')): return False
    if getattr( row, getVar(name,'JetCSVBtag')) > 0.8:  return False
    ###if getattr( row, getVar(name,'JetBtag')) > 3.3:     return False
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2:     return False
    return True
    
def tauSelection(row, name):
    if getattr( row, getVar(name,'Pt')) < 30:          return False
    if getattr( row, getVar(name,'AbsEta')) > 2.3:     return False
    if abs(getattr( row, getVar(name,'DZ'))) > 0.2:    return False
    return True


#VETOS
def vetos(row):
    if row.muVetoPt5IsoIdVtx: return False
    if row.eVetoMVAIsoVtx:    return False
    if row.eVetoCicTightIso:   return False # change it to loose
    if row.tauVetoPt20:        return False
    
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
    if isolabel == 'etauiso012' or isolabel == 'mutauiso012': 
        return bool( RelPFIsoDB < 0.12 ) 
    if isolabel == 'etauiso01' or isolabel == 'mutauiso01': 
        return bool( RelPFIsoDB < 0.1 ) 
    ##put the new iso
    if isolabel == 'mvaLoose' :
        if AbsEta < 0.8:
            return bool( RelPFIsoDB < 0.426 ) 
        if AbsEta > 0.8 and AbsEta < 1.479 :
            return bool(RelPFIsoDB < 0.481)
        if AbsEta > 1.479 and AbsEta < 2.5:
            return bool(RelPFIsoDB < 0.390)
        return False
    if isolabel == 'mvaTight':
        if AbsEta < 0.8:
            return bool( RelPFIsoDB < 0.105 ) 
        if AbsEta > 0.8 and AbsEta < 1.479 :
            return bool(RelPFIsoDB < 0.178)
        if AbsEta > 1.479 and AbsEta < 2.5:
            return bool(RelPFIsoDB < 0.150)
        return False
        
        

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


