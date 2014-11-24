#author Mauro Verzetti
'small interface module to deal with optimizization'

import os
import itertools

def make_0j_selector(tthr, ethr, e_t_DPhi, tMT):
    def selector(row):
        return row.tPt > tthr and \
             row.ePt > ethr and \
             row.e_t_DPhi > e_t_DPhi and \
             row.tMtToPfMet < tMT
    return selector

def make_1j_selector(tthr, ethr, tMT):
    def selector(row):
        return row.tPt > tthr and \
            row.ePt > ethr and \
            row.tMtToPfMet < tMT
    return selector

def make_2j_selector(tthr, ethr, tMT, vbfmass, vbfdeta):
    def selector(row):
        return row.tPt > tthr and \
            row.ePt > ethr and \
            row.tMtToPfMet < tMT and \
            row.vbfMass > vbfmass and \
            row.vbfDeta > vbfdeta
    return selector

RUN_OPTIMIZATION = ('RUN_OPTIMIZATION' in os.environ) and eval(os.environ['RUN_OPTIMIZATION'])

grid_search = {
    0 : {},
    1 : {},
    2 : {},
}

if RUN_OPTIMIZATION:
    t_thrs = range(20,70,10)+[35,45]
    e_thr  = range(20,70,10)+[35,45]

    mt_thr = range(0,50,10)+[35]
    dphi = [3.14, 3.00, 2.7, 2.4, 2.2]
    vbf_mass = range(400, 700, 100) + [550]
    vbf_deta = [3.0, 3.5, 4.0]
    
    gg_template = 't%i_e%i_dp%.1f_mt%i'
    for thresholds in itertools.product(t_thrs, e_thr, dphi, mt_thr):
        grid_search[0][gg_template % thresholds] = make_0j_selector(*thresholds)

    boost_template = 't%i_e%i_mt%i'
    for thresholds in itertools.product(t_thrs, e_thr, mt_thr):
        grid_search[1][boost_template % thresholds] = make_1j_selector(*thresholds)

    vbf_template = 't%i_e%i_mt%i_vbfm%i_vbfeta%i'
    for thresholds in itertools.product(t_thrs, e_thr, mt_thr, vbf_mass, vbf_deta):
        grid_search[2][vbf_template % thresholds] = make_2j_selector(*thresholds)


if __name__ == "__main__":
    from pdb import set_trace
    set_trace()
    #print '\n'.join(grid_search.keys())
else:
    print "Running optimization: %s" % RUN_OPTIMIZATION
