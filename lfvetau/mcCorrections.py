import os
import glob
import FinalStateAnalysis.TagAndProbe.HetauCorrection as HetauCorrection
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
from FinalStateAnalysis.PlotTools.decorators import memo, memo_last

@memo
def getVar(name, var):
    return name+var

is7TeV = bool('7TeV' in os.environ['jobid'])
pu_distributions  = {
    'singlee'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleElectron*pu.root'))}
pu_distributionsUp  = {
    'singlee'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleElectron*pu_up.root'))}
pu_distributionsDown  = {
    'singlee'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleElectron*pu_down.root'))}
mc_pu_tag                  = 'S6' if is7TeV else 'S10'


def make_puCorrector(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    weights = []
    if dataset in pu_distributions:# and dataset in pu_distributionsUp and dataset in pu_distributionsDown:
        return PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributions[dataset]))
#        weights = (PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributions[dataset])), PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributionsUp[dataset])), PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributionsDown[dataset])))
#        return weights
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')

def make_puCorrectorUp(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    if dataset in pu_distributions:
        return PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributionsUp[dataset]))
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')

def make_puCorrectorDown(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    if dataset in pu_distributions:
        return PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributionsDown[dataset]))
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')

def make_shifted_weights(default, shifts, functors):
    '''make_shifted_weights(default, shifts, functors) --> functor
    takes as imput the central value functor and two lists
    the name of the shifts and the shifted functors
    the returned functor takes one additional string to
    select the shifted functor. If the shift kwarg is missing
    or does not match any shift tag the central (default)
    fuctor output is returned'''
    #make default fetching faster
    default = default 
    def functor(*args, **kwargs):
        shift = ''
        if 'shift' in kwargs:
            shift = kwargs['shift']
            del kwargs['shift']

            #check if to apply shifts
            for tag, fcn in zip(shifts, functors):
                if tag == shift:
                    return fcn(*args, **kwargs)

        return default(*args, **kwargs)
    return functor

def make_multiple(fcn, indexed=False, shift=0):
    '''make_multiple(fcn, indexed=True, shift=0) --> functor
    takes as imput a weight correction function of pt and eta
    and returns a functor multiple(row,*args) --> weight
    where *args are the base name of the objects upon which 
    compute the correction.

    If indexed is true means that fcn returns a tuple 
    (weight, error) shift +/-1 makes the functor return the 
    weight+/-error in this case'''
    def multiple(row,*args):
        ret = 1.
        for arg in args:
            abseta = getattr( 
                    row, 
                    getVar(arg,'Eta')
                ) 
            pt     = getattr(row, getVar(arg,'Pt'))
            fcn_ret = fcn(pt,abseta)
            if indexed:
                value, err = fcn_ret
                if shift == 1:
                    ret *= (value + err)
                elif shift == -1:
                    ret *= (value - err)
                else:
                    ret *= value
            else:
                ret   *= fcn_ret
        return ret
    return multiple


##put here the trigger correction as in https://github.com/mverzett/UWHiggs/blob/WH_At_Paper/wh/mcCorrectors.py

correct_e             = make_multiple(HetauCorrection.correct_hamburg_e    )
correct_eid13_mva     = make_multiple(HetauCorrection.correct_eid13_mva    )
correct_eid13_p1s_mva = make_multiple(HetauCorrection.correct_eid13_p1s_mva)
correct_eid13_m1s_mva = make_multiple(HetauCorrection.correct_eid13_m1s_mva)

correct_eiso13_mva     = make_multiple(HetauCorrection.correct_eiso13_mva    )
correct_eiso13_p1s_mva = make_multiple(HetauCorrection.correct_eiso13_p1s_mva)
correct_eiso13_m1s_mva = make_multiple(HetauCorrection.correct_eiso13_m1s_mva)

#correct_eid_mva = make_multiple(HetauCorrection.scale_eleId_hww)
#correct_eReco_mva = make_multiple(HetauCorrection.scale_elereco_hww)
#correct_eIso_mva = make_multiple(HetauCorrection.scale_eleIso_hww)
correct_trigger_mva    = make_multiple(HetauCorrection.single_ele_mva, indexed=True)
correct_trigger_mva_up = make_multiple(HetauCorrection.single_ele_mva, indexed=True, shift=1)
correct_trigger_mva_dw = make_multiple(HetauCorrection.single_ele_mva, indexed=True, shift=-1)

efficiency_trigger_mva    = make_multiple(HetauCorrection.single_ele_eff_mva, indexed=True)
efficiency_trigger_mva_up = make_multiple(HetauCorrection.single_ele_eff_mva, indexed=True, shift=1)
efficiency_trigger_mva_dw = make_multiple(HetauCorrection.single_ele_eff_mva, indexed=True, shift=-1)

eiso_correction = make_shifted_weights(
    correct_eiso13_mva, 
    ['eisop1s','eisom1s'], 
    [correct_eiso13_p1s_mva, correct_eiso13_m1s_mva],
)

eid_correction = make_shifted_weights(
    correct_eid13_mva,
    ['eidp1s','eidm1s'],
    [correct_eid13_p1s_mva, correct_eid13_m1s_mva]
)

trig_correction = make_shifted_weights(
    correct_trigger_mva,
    ['trp1s', 'trm1s'],
    [correct_trigger_mva_up, correct_trigger_mva_dw]
)

trig_efficiency = make_shifted_weights(
    efficiency_trigger_mva,
    ['trp1s', 'trm1s'],
    [efficiency_trigger_mva_up, efficiency_trigger_mva_dw]
)
