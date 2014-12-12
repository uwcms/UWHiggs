 #author Mauro Verzetti
'small interface module to deal with optimizization'

import os
import itertools

RUN_OPTIMIZATION = ('RUN_OPTIMIZATION' in os.environ) and eval(os.environ['RUN_OPTIMIZATION'])

_0jets = {
    'tPt'  : range(30,50,10)+[35,45],
    'ePt'  : range(30,70,10)+[35,45],
    'dphi' : [3.14, 3.00, 2.7, 2.4, 2.2],
    'tMtToPfMet' : range(20,50,10)+[35],
}
_0jet_region_template = 'tPt%i_ePt%i_dphi%.2f_tMtToPfMet%i'
def _get_0jet_regions(tPt, ePt, dphi, tMtToPfMet):
    pass_tPt        = [i for i in _0jets['tPt'       ] if tPt        > i] 
    pass_ePt        = [i for i in _0jets['ePt'       ] if ePt        > i] 
    pass_dphi       = [i for i in _0jets['dphi'      ] if dphi       > i] 
    pass_tMtToPfMet = [i for i in _0jets['tMtToPfMet'] if tMtToPfMet < i] 
    return [_0jet_region_template % i for i in itertools.product(pass_tPt, pass_ePt, pass_dphi, pass_tMtToPfMet)]


_1jets = {
    'tPt'  : range(30,60,10)+[35,45],
    'ePt'  : range(30,70,10)+[35,45],
    'tMtToPfMet' : range(20,50,10)+[35],
}
_1jet_region_template = 'tPt%i_ePt%i_tMtToPfMet%i'
def _get_1jet_regions(tPt, ePt, tMtToPfMet):
    pass_tPt        = [i for i in _1jets['tPt'       ] if tPt        > i]
    pass_ePt        = [i for i in _1jets['ePt'       ] if ePt        > i]
    pass_tMtToPfMet = [i for i in _1jets['tMtToPfMet'] if tMtToPfMet < i]
    return [_1jet_region_template % i for i in itertools.product(pass_tPt, pass_ePt, pass_tMtToPfMet)]

_2jets = {
    'tPt'  : range(30,50,10)+[35,45],
    'ePt'  : range(30,70,10)+[35,45],
    'tMtToPfMet' : range(0,50,10)+[35],
    'vbf_mass' : range(400, 600, 100) + [550],
    'vbf_deta' : [2.5, 3.0, 3.5, 4.0],
}
_2jet_region_template = 'tPt%i_ePt%i_tMtToPfMet%i_vbf_mass%i_vbf_deta%.1f'
def _get_2jet_regions(tPt, ePt, tMtToPfMet, vbf_mass, vbf_deta):
    pass_tPt        = [i for i in _2jets['tPt'       ] if tPt        > i]
    pass_ePt        = [i for i in _2jets['ePt'       ] if ePt        > i]
    pass_tMtToPfMet = [i for i in _2jets['tMtToPfMet'] if tMtToPfMet < i]
    pass_vbf_mass = [i for i in _2jets['vbf_mass'] if vbf_mass > i]
    pass_vbf_deta = [i for i in _2jets['vbf_deta'] if vbf_deta > i]
    return [_2jet_region_template % i for i in itertools.product(pass_tPt, pass_ePt, pass_tMtToPfMet, pass_vbf_mass, pass_vbf_deta)]

def empty(*args):
    return []

compute_regions_0jet = _get_0jet_regions if RUN_OPTIMIZATION else empty
compute_regions_1jet = _get_1jet_regions if RUN_OPTIMIZATION else empty
compute_regions_2jet = _get_2jet_regions if RUN_OPTIMIZATION else empty

regions = {'0' : [], '1' : [], '2' : [], '3' : []}
if RUN_OPTIMIZATION:
    regions = {
        '0' : [_0jet_region_template % i for i in itertools.product(_0jets['tPt'], _0jets['ePt'], _0jets['dphi'], _0jets['tMtToPfMet'])],
        '1' : [_1jet_region_template % i for i in itertools.product(_1jets['tPt'], _1jets['ePt'], _1jets['tMtToPfMet'])],
        '2' : [_2jet_region_template % i for i in itertools.product(_2jets['tPt'], _2jets['ePt'], _2jets['tMtToPfMet'], _2jets['vbf_mass'], _2jets['vbf_deta'])],
        '3' : []}

if __name__ == "__main__":
    from pdb import set_trace
    set_trace()
    #print '\n'.join(grid_search.keys())
else:
    print "Running optimization: %s" % RUN_OPTIMIZATION
