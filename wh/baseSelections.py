from FinalStateAnalysis.PlotTools.decorators import memo

@memo
def getVar(name, var):
    return name+var

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

def vetos(row):
    if row.muVetoPt5:         return False
    if row.bjetCSVVeto:       return False
    if row.tauVetoPt20:       return False
    if row.eVetoCicTightIso:  return False
    return True

def control_region_ee(row):
    '''Figure out what control region we are in. Shared among two codes, to avoid mismatching copied here'''
    if  row.e1_e2_SS and row.e1RelPFIsoDB < 0.15 and row.e1MtToMET > 30 and row.e2MtToMET < 30:# and row.e2MtToMET < 30: #and row.metEt > 30: #row.metSignificance > 3:
        return 'wjets'
    elif row.e1_e2_SS and row.e1RelPFIsoDB > 0.3 and row.metEt < 25: #and row.metSignificance < 3: #
        return 'qcd'
    elif row.e1RelPFIsoDB < 0.1 and row.e2RelPFIsoDB < 0.1 and row.e2MVAIDH2TauWP and row.e1MVAIDH2TauWP \
        and not any([ row.muVetoPt5,
                      row.tauVetoPt20,
                      row.eVetoCicTightIso,
                      ]):
        return 'zee'
    else:
        return None


