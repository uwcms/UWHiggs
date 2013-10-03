import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor, make_corrector_from_th2, build_uncorr_2Droofunctor, FunctorFromMVA
from FinalStateAnalysis.StatTools.VariableScaler import make_scaler
from FinalStateAnalysis.Utilities.smartdict import SmartDict
import itertools
import optimizer
import re
#from optimizer import leading_lepton_iso_tag, subleading_lepton_iso_tag

################################################################################
#### Fitted fake rate functions ################################################
################################################################################

# Get fitted fake rate functions
frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')

def make_simple_mapper(tomap):
    def f_(string):
        ret = string
        for key, val in tomap.iteritems():
            ret = ret.replace(key, val)
        return ret
    return f_

def make_regex_mapper(tomap):
    def f_(string):
        ret = string
        for key, val in tomap.iteritems():
            ret = re.sub(key, val, ret)
        return ret
    return f_


def build_roofunctor_dict(basename, wsname = 'fit_efficiency', functor_name = 'efficiency', mapper = None):
    ret = {}
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret[i] = build_roofunctor(
            basename % lepid,
            wsname, # workspace name
            functor_name
            )
    return ret

def build_2Droofunctor_dict(basename, 
                            xvar, yvar,
                            wsname = 'fit_efficiency', functor_name = 'efficiency', mapper = None):
    ret = {}
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret[i] = build_uncorr_2Droofunctor(
            build_roofunctor(
                basename % (lepid, xvar),
                wsname, # workspace name
                functor_name
            ),
            build_roofunctor(
                basename % (lepid, yvar),
                wsname, # workspace name
                functor_name
            ),
            (basename % (lepid, yvar)).replace('.root','.corrected_inputs.root'),
            )
    return ret


def make_corrector_dict(filename, mapname, mapper=None):
    ret = {}
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret[i] = make_corrector_from_th2(filename % lepid, mapname)
    return ret
    
def make_scaler_dict(filename, mapname):
    ret = {}
    for i in optimizer.lep_id:
        ret[i] = make_scaler(filename % i, 'mass_scale')
    return ret

def make_mva_functor_dict(template, variables, mapper=None):
    ret = SmartDict() #faster than normal dict
    for i in optimizer.lep_id:
        lepid = i
        if mapper:
            lepid = mapper(lepid)
        ret.book(i, FunctorFromMVA,
            template % lepid,
            template % lepid,
            *variables
        )
    return ret


##################
## 1D Muon Func ##
##################

#no changes in muonID in 2013
mapper    = {'eid1[0-9][A-Z][a-z]+_':'', 'idiso02' : 'pfidiso02'}
variables = ['muonJetPt', 'muonPt', 'muonJetCSVBtag'] #'muonPVDXY']#, 'muonJetBtag']
highpt_mu_fr = make_mva_functor_dict(frfit_dir + '/m_wjets_pt20_%s_muonInfo.kNN.weights.xml', variables, mapper=make_regex_mapper(mapper))
lowpt_mu_fr  = make_mva_functor_dict(frfit_dir + '/m_wjets_pt10_%s_muonInfo.kNN.weights.xml', variables, mapper=make_regex_mapper(mapper))

highpt_mu_qcd_fr = make_mva_functor_dict(frfit_dir + '/m_qcd_pt20_%s_muonInfo.kNN.weights.xml', variables, mapper=make_regex_mapper(mapper))
lowpt_mu_qcd_fr  = make_mva_functor_dict(frfit_dir + '/m_qcd_pt10_%s_muonInfo.kNN.weights.xml', variables, mapper=make_regex_mapper(mapper))


#######################
## 1D Electrons Func ##
#######################

variables = ['electronJetPt', 'electronPt']
#EMT
highpt_e_fr = make_mva_functor_dict(frfit_dir + '/e_wjets_pt20_%s_electronInfo.kNN.weights.xml', variables)
lowpt_e_fr  = make_mva_functor_dict(frfit_dir + '/e_wjets_pt10_%s_electronInfo.kNN.weights.xml', variables)

highpt_e_qcd_fr = make_mva_functor_dict(frfit_dir + '/e_qcd_pt20_%s_electronInfo.kNN.weights.xml', variables)
lowpt_e_qcd_fr  = make_mva_functor_dict(frfit_dir + '/e_qcd_pt10_%s_electronInfo.kNN.weights.xml', variables)

#EET
highpt_ee_fr = make_mva_functor_dict(frfit_dir + '/ee_wjetsNoZmass_pt20_%s_electronInfo.kNN.weights.xml', variables) 
lowpt_ee_fr  = make_mva_functor_dict(frfit_dir + '/ee_wjetsNoZmass_pt10_%s_electronInfo.kNN.weights.xml', variables) 

highpt_ee_qcd_fr = make_mva_functor_dict(frfit_dir + '/ee_qcd_pt20_%s_electronInfo.kNN.weights.xml', variables) 
lowpt_ee_qcd_fr  = make_mva_functor_dict(frfit_dir + '/ee_qcd_pt10_%s_electronInfo.kNN.weights.xml', variables) 


##################
## 1D Taus Func ##
##################

tau_fr = build_roofunctor(
    frfit_dir + '/t_ztt_pt20_mvaloose_tauPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

tau_qcd_fr = tau_fr 

e_charge_flip      = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_%s.root", "efficiency_map")         
e_charge_flip_up   = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_%s.root", "efficiency_map_statUp")  
e_charge_flip_down = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_%s.root", "efficiency_map_statDown")
mass_scaler        = make_scaler_dict(frfit_dir+"/charge_flip_prob_map_%s.root", 'mass_scale')
default_scaler     = mass_scaler[mass_scaler.keys()[0]]


e1_charge_flip      = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_e1_%s.root", "efficiency_map")         
e1_charge_flip_up   = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_e1_%s.root", "efficiency_map_statUp")  
e1_charge_flip_down = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_e1_%s.root", "efficiency_map_statDown")

e2_charge_flip      = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_e2_%s.root", "efficiency_map")         
e2_charge_flip_up   = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_e2_%s.root", "efficiency_map_statUp")  
e2_charge_flip_down = make_corrector_dict(frfit_dir+"/charge_flip_prob_map_e2_%s.root", "efficiency_map_statDown")
