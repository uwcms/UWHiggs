#!/usr/bin/env python

from UWHiggs.hzg.datacard.directory_prep import directory_prep
from UWHiggs.hzg.datacard.metadata_association import metadata_association

import ROOT
from ROOT import RooWorkspace

dp = directory_prep('~/HZG_analysis','08JAN2013_vanilla')
dp.build_inputs()

mda = metadata_association(dp)

def process_inputs(meta_data):
    for (sample,chanlist) in meta_data.getAssociation().iteritems():
        for (channel,proclist) in chanlist.iteritems():
            for (process,subproclist) in proclist.iteritems():
                for (subproc,info) in subproclist.iteritems():
                    print sample, channel, process, subproc

                    
                    
                    print '\t',
                    if 'data' not in process:                    
                        if info['isHiggs']:
                            print info['mass'],
                        print info['x_sec'],
                    print info['input_file'].split('/')[-1]

if __name__ == '__main__':
    process_inputs(mda)
                
    

