import os
import glob
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight

# Determine MC-DATA corrections
is7TeV = bool('7TeV' in os.environ['jobid'])
print "Is 7TeV:", is7TeV


########################################################################
##
##                         PILE UP CORRECTIONS
##
########################################################################

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
pu_distributions  = {
    'doublemu' : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_DoubleMu*pu.root')),
    'doublee'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_DoubleElectron*pu.root')),
    'mueg'     : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_MuEG*pu.root')),}
mc_pu_tag                  = 'S6' if is7TeV else 'S10'
# Hack to use S6 weights for the HWW 7TeV sample we use in 8TeV
if 'HWW3l' in os.environ.get('megatarget', 'NOTSET') and not is7TeV:
    print "Using S6_600bins PU weights for HWW3l"
    mc_pu_tag = 'S6_600bins'

def make_puCorrector(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    if dataset in pu_distributions:
        return PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributions[dataset]))
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')
## def force_pu_distribution(kind):
##     pu_corrector = PileupWeight.PileupWeight( kind, *pu_distributions)



########################################################################
##
##                              TRIGGER
##
########################################################################

#DOUBLE MU
muon_pog_Mu17Mu8_2011      = MuonPOGCorrections.muon_pog_Mu17Mu8_eta_eta_2011
#muon_pog_Mu17Mu8_Mu17_2012 = MuonPOGCorrections.make_muon_pog_Mu17Mu8_Mu17_2012()
#muon_pog_Mu17Mu8_Mu8_2012  = MuonPOGCorrections.make_muon_pog_Mu17Mu8_Mu8_2012()
muon_h2tau_Mu17Mu8_2012    = H2TauCorrections.correct_double_muon_trg_2012

def double_muon_trigger(row,m1,m2):
    'makes scale factor for double mu trigger'
    if is7TeV:
        return muon_pog_Mu17Mu8_2011(getattr(row, '%sEta' % m1), getattr(row, '%sEta' % m2) )
    else:
        return muon_h2tau_Mu17Mu8_2012(getattr(row, '%sPt' % m1), getattr(row, '%sEta' % m1), getattr(row, '%sPt' % m2), getattr(row, '%sEta' % m2))


#MUEG
correct_mueg_mu            = H2TauCorrections.correct_mueg_mu_2011 if is7TeV else H2TauCorrections.correct_mueg_mu_2012
correct_mueg_e             = H2TauCorrections.correct_mueg_e_2011  if is7TeV else H2TauCorrections.correct_mueg_e_2012

#Double electrons scale factors    
correct_double_electron    = H2TauCorrections.correct_double_electron_trg_2011 if is7TeV else H2TauCorrections.correct_double_electron_trg_2012

def double_electron_trigger(row):
    return correct_double_electron( row.e1Pt, row.e1Eta, row.e2Pt, row.e2Eta )



########################################################################
##
##                              ID/ISO
##
########################################################################


#Makes appropriate correction function for electrons or muons according to run period
mu_pog_2011_id             = MuonPOGCorrections.make_muon_pog_PFTight_2011()
mu_pog_2011_iso            = MuonPOGCorrections.make_muon_pog_PFRelIsoDB02_2011()
muon_pog_IsoID             = (lambda pt, eta: mu_pog_2011_id(pt,eta)*mu_pog_2011_iso(pt,eta)) if is7TeV else H2TauCorrections.correct_mu_idiso_2012
electron_corrections       = H2TauCorrections.correct_e_idiso_2011 if is7TeV else H2TauCorrections.correct_e_idiso_2012
electron_tight_corrections = H2TauCorrections.correct_e_TIGHTidiso_2011 if is7TeV else H2TauCorrections.correct_e_TIGHTidiso_2012


def get_muon_corrections(row,*args):
    'makes corrections to iso and id of muons'
    ret = 1.
    for arg in args:
        eta = getattr(row, '%sEta' % arg)
        pt  = getattr(row, '%sPt'  % arg)
        #print muon_pog_IsoID.__name__
        ret *= muon_pog_IsoID( pt, eta)
    return ret

def get_electron_corrections(row,*args):
    'makes corrections to iso and id of electrons'
    ret = 1.
    for arg in args:
        abseta = abs(getattr(row, '%sEta' % arg))
        pt     = getattr(row, '%sPt'  % arg)
        ret   *= electron_corrections(pt,abseta)
    return ret

