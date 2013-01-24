#!/usr/bin/env python

'''
File: submit_hzg_analysis.py

Author: L. Gray, FNAL

Description: submit hzg_analyzer jobs by condor

'''

from RecoLuminosity.LumiDB import argparse
import fnmatch
from UWHiggs.hzg.hzg_metainfo import analysis_list, ntupleRoot, ntupleRevision
from FinalStateAnalysis.Utilities.version import fsa_version
import os
import sys


parser = argparse.ArgumentParser(description='Build analysis submission')
parser.add_argument('jobid', help='Job ID identifier')
parser.add_argument('--samples', nargs='+', type=str, required=False,
                    help='Filter samples using list of patterns (shell style)')
parser.add_argument('--subsamples',nargs='+', type=str, required=False,
                    help='Filter subsamples using list of patterns')
parser.add_argument('--vanilla', default=False,action='store_true',
                    help='run with no corrections')
parser.add_argument('--leptonType', type=str, required=True,
                    help='the lepton type')
parser.add_argument('--leptonCor', type=str, required=False,
                    help='the lepton corrections')
parser.add_argument('--photonCor', type=str, required=False,
                    help='the photon corrections')

args = parser.parse_args()

cfg = '%s/src/UWHiggs/hzg/scripts/hzg_analysis.py'\
      %os.environ['CMSSW_BASE']
jobId = args.jobid

def filter_sample(filter,sample):
    if filter:
        passes_wildcard = False
        for pattern in filter:
            if fnmatch.fnmatchcase(sample, pattern):
                passes_wildcard = True
        return passes_wildcard 
    else:
        return True

print " # Job ID: %s Version: %s" % (jobId, fsa_version())
print 'export TERMCAP=screen'
for sample in sorted(analysis_list.keys()):
    subsamples = analysis_list[sample]

    passes_filter = filter_sample(args.samples,sample)
    
    if not passes_filter:
        continue
    for subsample in subsamples.keys():
        if( subsample != 'tupleName' and
            subsample != 'tupleDate' and
            subsample != 'tupleRoot' ):

            passes_filter = filter_sample(args.subsamples,
                                          subsample)    
            if not passes_filter:
                continue

            for tuple in subsamples[subsample]:
                tuple_info = subsamples[subsample][tuple]
                tupleName = tuple
            
                submit_dir_base = "/scratch/{logname}/HZG_analysis/{jobid}"\
                                  "/{sample}.{subsample}.{tupleName}."\
                                  "{leptonType}".format(
                    logname = os.environ['LOGNAME'],
                    jobid = jobId,
                    sample = sample,
                    subsample = subsample,
                    tupleName = tupleName,
                    leptonType = args.leptonType
                )            
                
                dag_directory = os.path.join(submit_dir_base, 'dags')
                # Create the dag directory
                print "mkdir -p %s" % dag_directory
                submit_dir = os.path.join(submit_dir_base, 'submit')

                sys.stderr.write('Building sample submit dir %s\n' % (sample))
                if os.path.exists(submit_dir):
                    sys.stderr.write('=> skipping existing submit'\
                                     ' directory for %s\n' % sample)
                    continue            

                options = ['++leptonType=%s'%(args.leptonType)]
                if '7TeV' in sample:
                    options.append('++runYear=2011')
                    if 'data' in subsample:
                        if '2011A' in subsample:
                            options.append('++runType=A')
                        if '2011B' in subsample:
                            options.append('++runType=B')
                        options.append('++datType=data')
                    else:
                        options.append('++runType=AB')
                        options.append('++datType=mc')
                        options.append('++crossSection=%.5e'%tuple_info['x_sec'])
                else:
                    options.append('++runYear=2012')
                    if 'data' in subsample:
                        options.append('++runType=ABCD')
                        options.append('++datType=data')
                    else:
                        options.append('++runType=ABCD')
                        options.append('++datType=mc')
                        options.append('++crossSection=%.5e'%tuple_info['x_sec'])
                if args.leptonCor:
                    options.append('++leptonCor=%s'%args.leptonCor)

                if args.photonCor:
                    options.append('++photonCor=%s'%args.photonCor)
                    
                if args.vanilla:
                    options.append('++vanilla')
                    
                if tuple_info['vetoIFSR']:
                    options.append('++vetoIFSR')
                options.append("'$inputFileNames'")

                
                farmout_options = []
                farmout_options.append(
                    '--input-dir=root://cmsxrootd.hep.wisc.edu/%s/%s/'\
                    '%s-%s/%s/'%(ntupleRoot,
                                 ntupleRevision,
                                 subsample,
                                 sample,
                                 '.'.join(tuple_info['datasetpath'][1:].split('/')))
                    )
                    
                
                output_dir = '/hdfs/store/user/%s/HZG_analysis/%s/%s/%s/%s'\
                             %(os.environ['LOGNAME'],
                               jobId,
                               args.leptonType,
                               '-'.join([subsample,sample]),
                               '.'.join(tuple_info['datasetpath'][1:].split('/'))
                               )
                
                command = [
                    'farmoutAnalysisJobs',
                    '--job-generates-output-name',
                    '--fwklite',
                    #'--no-shared-fs', # Copy libs to submit dir so we don't kill AFS
                    '--infer-cmssw-path',
                    '--vsize-limit=30000',
                    '--express-queue',
                    '--input-files-per-job=100',
                    '"--output-dir=srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN=%s"' % output_dir,
                    '--submit-dir=%s' % submit_dir,
                    '--output-dag-file=%s/dag.dag' % dag_directory,
                    ]
                
                command.extend(farmout_options)
                
                command.append('-'.join([jobId, tupleName]))
                command.append(cfg)
                command.extend(options)
                print ' '.join(command)
