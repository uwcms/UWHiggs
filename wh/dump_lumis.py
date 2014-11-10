import glob
import os
from FinalStateAnalysis.MetaData.data_views import extract_sample, read_lumi
import re
#from pudb import set_trace; set_trace()

paths = {
    #os.environ['current'],
    'current'  : 'inputs/%s/' % os.environ['jobid'],
    'previous' : 'inputs/2013-Jun-30-8TeV/',
}


for dataset in ['DoubleMu', 'DoubleElectron', 'MuEG']:
    lumi_dict = {}
    for jobid, path in paths.iteritems():
        lumi_dict[jobid] = {}
        for lumifile in glob.glob(path+'/data_'+dataset+'*.lumicalc.sum'):
            lumi_dict[jobid][extract_sample(lumifile)] = read_lumi(lumifile)

    print dataset
    print '%60s%20s%20s' % ('', 'current', 'previous')
    keys = []
    for i in lumi_dict.values():
        keys += i.keys()
    keys = list(set(keys))
    
    total_lumis = dict([(i, {'current' : 0., 'previous' : 0.}) for i in ['2012A', '2012B', '2012C', '2012D', 'TOTAL',]])
    
    for sample in keys:
        curr_l = lumi_dict['current'][sample] if sample in lumi_dict['current'] else 0.
        total_lumis['TOTAL']['current'] += curr_l
        previous_l  = lumi_dict['previous'][sample]  if sample in lumi_dict['previous'] else 0.
        total_lumis['TOTAL']['previous'] += previous_l
        for key in total_lumis:
            if key in sample:
                total_lumis[key]['previous'] += previous_l
                total_lumis[key]['current'] += curr_l
        print '%60s%20s%20s' % (sample, '%.1f' % curr_l, '%.1f' %  previous_l)

    for key, val in total_lumis.iteritems():
        print '%60s%20s%20s' % (key, '%.1f' % val['current'], '%.1f' %  val['previous'])
    print '\n\n'
