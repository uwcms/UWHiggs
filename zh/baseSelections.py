'''
Common selection used in ZH analysis
'''

import os
from FinalStateAnalysis.PlotTools.decorators import memo

@memo
def getVar(name, var):
    return name+var

# Determine MC-DATA corrections
is7TeV = bool('7TeV' in os.environ['jobid'])

iso_base_name = 'RelPFIsoDB' #'RelPFIsoDBZhLike'
def objIsolation(row,muon):
    return getattr(row, getVar(muon, iso_base_name) )

def muIsoTight(row,muon):
    return objIsolation(row,muon) < 0.15

def muIsoLoose(row,muon):
    return objIsolation(row,muon) < 0.25

def elIsoTight(row,el):
    return objIsolation(row,el) < 0.10

def elIsoLoose(row,el):
    return objIsolation(row,el) < 0.25


def MuTriggerMatching(row):
    '''
    Applies trigger matching according to the run period
    '''
    if is7TeV:
        return row.m1MatchesDoubleMu2011Paths > 0 and row.m2MatchesDoubleMu2011Paths > 0
    else:
        return ( row.m1MatchesDoubleMu2011Paths > 0 or row.m1MatchesMu17TrkMu8Path > 0 ) and ( row.m2MatchesDoubleMu2011Paths > 0 or row.m2MatchesMu17TrkMu8Path > 0 )
 
def ElTriggerMatching(row):
    '''
    Applies trigger matching
    '''
    return row.e1MatchesDoubleEPath > 0 and row.e2MatchesDoubleEPath > 0

def Vetos(row):
    '''
    applies b-tag, muon, electron and tau veto
    '''
    if bool(row.bjetCSVVeto):      return False
    #if bool(row.bjetCSVVetoZHLikeNoJetId): return False
    if bool(row.muGlbIsoVetoPt10): return False
    if bool(row.tauHpsVetoPt20):   return False
    if bool(row.eVetoMVAIso):      return False
    return True

def overlap(row,*args):
    return any( map( lambda x: x < 0.1, [getattr(row,'%s_%s_DR' % (l1,l2) ) for l1 in args for l2 in args if l1 <> l2 and hasattr(row,'%s_%s_DR' % (l1,l2) )] ) )

def eleID(row, name):
    if getattr(row, getVar(name, 'Pt') )   < 10: return False
    if getattr(row, getVar(name, 'AbsEta') ) < 0.8    and getattr(row, getVar(name, 'MVANonTrig')) > 0.5: return True
    if getattr(row, getVar(name, 'AbsEta') ) >= 0.8   and getattr(row, getVar(name, 'AbsEta')) < 1.479 and getattr(row, getVar(name, 'MVANonTrig')) > 0.12: return True
    if getattr(row, getVar(name, 'AbsEta') ) >= 1.479 and getattr(row, getVar(name, 'MVANonTrig')) > 0.6: return True
    return False

def ZMuMuSelectionNoVetos(row):
    '''
    Z Selection as AN
    '''
    #Z Selection
    if not (row.doubleMuPass or row.doubleMuTrkPass):  return False
    if row.m1Pt < row.m2Pt:                            return False
    if row.m1Pt < 20:                                  return False
    if row.m2Pt < 10:                                  return False
    if row.m1AbsEta > 2.4:                             return False
    if row.m2AbsEta > 2.4:                             return False
    if abs(row.m1DZ) > 0.1:                            return False
    if abs(row.m2DZ) > 0.1:                            return False
    #if not bool(row.m1PFIDTight):                      return False
    #if not muIsoLoose(row,'m1'):                       return False
    #if not bool(row.m2PFIDTight):                      return False
    #if not muIsoLoose(row,'m2'):                       return False
    #if bool(row.m1_m2_SS):                             return False
    if row.m1_m2_Mass < 60 or row.m1_m2_Mass > 120 :   return False
    return True
#return MuTriggerMatching(row)

def ZMuMuSelection(row):
    return ZMuMuSelectionNoVetos(row) and Vetos(row)

def ZEESelectionNoVetos(row):
    '''
    Z Selection as AN
    '''
    if not row.doubleEPass:                          return False
    if row.e1Pt < row.e2Pt:                          return False
    if row.e1Pt < 20:                                return False
    if row.e2Pt < 10:                                return False
    if row.e1AbsEta > 2.5:                           return False
    if row.e2AbsEta > 2.5:                           return False
    if abs(row.e1DZ) > 0.1:                          return False
    if abs(row.e2DZ) > 0.1:                          return False
    #if not eleID(row, 'e1'):                         return False
    #if not elIsoLoose(row, 'e1'):                    return False
    #if not eleID(row, 'e2'):                         return False
    #if not elIsoLoose(row, 'e1'):                    return False
    #if bool(row.e1_e2_SS):                           return False
    if row.e1_e2_Mass < 60 or row.e1_e2_Mass > 120 : return False
    return True
#return ElTriggerMatching(row)

def ZEESelection(row):
    return ZEESelectionNoVetos(row) and Vetos(row)

def signalMuonSelection(row,muId):
    '''
    Basic selection for signal muons (the ones coming from Higgs). No Isolation applied
    '''
    if getattr(row, getVar(muId,'Pt') ) < 10:              return False
    if getattr(row, getVar(muId,'AbsEta') ) > 2.4:         return False
    if abs(getattr(row, getVar(muId,'DZ') )) > 0.1:        return False
        #if not bool(getattr(row, '%sPFIDTight' % muId) ): return False
    return True

def signalTauSelection(row, tauId, ptThr = 20):
    '''
    Basic selection for signal hadronic (the ones coming from Higgs). No Isolation is applied, but DecayMode is
    '''
    if not bool( getattr( row, getVar(tauId, 'DecayFinding') ) ):      return False
    if getattr( row, getVar(tauId, 'Pt') )  < ptThr:                   return False
    if getattr( row, getVar(tauId, 'AbsEta') )  > 2.3:                 return False
    if abs(getattr( row, getVar(tauId, 'DZ') ) ) > 0.1:                return False
    if abs(getattr( row, getVar(tauId, 'JetCSVBtag') ) ) > 0.898:      return False #Veto Btagged taus
    return True


def signalElectronSelection(row, elId):
    '''
    Basic selection for signal electrons (the ones coming from Higgs). No Isolation applied
    '''
    if getattr(row, getVar(elId, 'Pt') ) < 10:                 return False
    if getattr(row, getVar(elId, 'AbsEta') ) > 2.5:            return False
    if abs(getattr(row, getVar(elId, 'DZ') )) > 0.1:           return False
        #if not bool(getattr(row, '%sMVAIDH2TauWP' % elId) ): return False
    return True
    
