'small interface module to deal with optimizization'

import os
import itertools

RUN_OPTIMIZATION = 'RUN_OPTIMIZATION' in os.environ and bool(os.environ['RUN_OPTIMIZATION'])

print "Running optimization: %s" % RUN_OPTIMIZATION
lep_id = [
    'h2taucuts',
    'h2taucuts020',
    'idiso02'
    ] if RUN_OPTIMIZATION else ['h2taucuts']

LT_cut = [60, 70, 80, 90]

grid_search = {}
if RUN_OPTIMIZATION:
    for cut_combo in itertools.product(lep_id, lep_id, LT_cut):
        cut_name = '%s:%s:%s' % cut_combo
        grid_search[cut_name] = {
            'leading_iso' : cut_combo[0],
            'subleading_iso' : cut_combo[1],
            'LT'   : cut_combo[2],
        }
else:
    grid_search[''] = {
        'leading_iso'    : 'h2taucuts',
        'subleading_iso' : 'h2taucuts',
        'LT'             : 80,
    }
