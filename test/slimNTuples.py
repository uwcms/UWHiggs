#! /bin/env python

import os
import sys
import subprocess
import re
from glob import glob
from fnmatch import fnmatch
import time

if len(sys.argv) < 2 or sys.argv[-1] is '-h' or sys.argv[-1] is '--help':
    print 'Usage ./slimNTuples other_branches_to_keep.list /hdfs/dir/with/files/to/be/skimmed'
    sys.exit(1)

branches_to_keep = sys.argv[-2]
hdfs_path        = sys.argv[-1]
if not os.path.isfile(branches_to_keep):
    print "Error! %s: no such file" % branches_to_keep
    sys.exit(1)

if not os.path.isdir(hdfs_path):
    print "Error! %s: no such directory" % hdfs_path
    sys.exit(1)

print os.getcwd()
print "looking in the code to find used branches..."

#look for row.*
proc = subprocess.Popen("grep 'row\.' %s/src/UWHiggs/?h/*.py" % (os.environ['CMSSW_BASE']), shell=True, stdout=subprocess.PIPE)
stdout, stderr = proc.communicate()
exitc  = proc.wait()
regex = re.compile(r'row\.\w+')
usedBranches = list( set([i.split('.')[-1] for i in regex.findall(stdout)]))

#look for more complex statements (getVar and getattr)
proc = subprocess.Popen("grep getattr %s/src/UWHiggs/?h/*.py | grep row" % (os.environ['CMSSW_BASE']), shell=True, stdout=subprocess.PIPE)
stdout, stderr = proc.communicate()
exitc  = proc.wait()
regex = re.compile(r'getattr\( *row *,[^\'\"]*[\'\"](?P<branch>(?:\w|\%)+)')
otherObjs = [m.group('branch') for m in regex.finditer(stdout) if m]
regex = re.compile(r'getVar\( *\w+ *,[^\'\"]*[\'\"](?P<branch>(?:\w|\%)+)')
otherObjs.extend([m.group('branch') for m in regex.finditer(stdout) if m])
otherObjs = [ i.replace('%s','*') if '%s' in i else '*'+i for i in otherObjs ]

#free some memory
del stdout

#get rid of doubles
otherObjs = list( set(otherObjs) )
#check that our starred branches don't match any non starred, in case of match discard the non starred
usedBranches = filter(lambda x: not any( ( fnmatch(x,y) for y in otherObjs) ), usedBranches) #use generator to be faster
usedBranches.extend(otherObjs)

JOBID   = filter(lambda x: x is not '', hdfs_path.split('/'))[-1]
SAMPLES = [i.split('/')[-1] for i in glob(hdfs_path+'/*')]

sample_tfile = glob('/'.join([hdfs_path,SAMPLES[0],'*.root']))[0]
print 'Using %s as sample file to find trees and matching branches...' % sample_tfile

import ROOT
def GetContent(dir):
    tempList = dir.GetListOfKeys()
    retList = []
    for it in range(0,tempList.GetSize()):
       retList.append(tempList.At(it).ReadObj())
    return retList

def FindTrees( directory, dirName, objectList ):
    dirContent = GetContent(directory)
    for entry in dirContent:
        if type(entry) is ROOT.TDirectory or type(entry) is ROOT.TDirectoryFile:
            subdirName = os.path.join(dirName,entry.GetName())
            FindTrees(entry, subdirName,objectList)
        elif entry.InheritsFrom('TTree'):
            pathname = os.path.join(dirName,entry.GetName())
            objectList.append(pathname)

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

#Open sample file
tfile = ROOT.TFile.Open(sample_tfile)

#find trees inside
forest = []
FindTrees( tfile, '', forest )

timestamp = str(int(time.mktime(time.gmtime())))
listsdir  = '/scratch/mverzett/%s' % timestamp
## make_dir(listsdir)

#define matching function, use generator to be faster (should stop at the first match)
match  = lambda x:  any( ( fnmatch(x,y) for y in usedBranches) )
## with open(listsdir+'/trees_location.list','w') as f:
##     f.write('\n'.join(forest))

tot_branches  = 0
kept_branches = 0
for tree_name in forest:
    if 'Ntuple' in tree_name:
        tree = tfile.Get(tree_name)
        #Get All the branches
        matching_branches = [branch.GetName() for branch in tree.GetListOfBranches()]
        tot_branches     += len(matching_branches)
        matching_branches = filter(match, matching_branches)
        kept_branches    += len(matching_branches)
        ## with open(listsdir+'/%s.list' % tree_name.replace('/','_'),'w') as f:
        ##     f.write('\n'.join(matching_branches))

compression_ratio = float(kept_branches) / float(tot_branches)
#os.path.getsize(
print 'Compression ratio: %s' % compression_ratio
print "Beginning multi threaded copy/complession/merge..."
approx_chunk_size = 100*10**6 #~100MB
for sample in SAMPLES:
    print 'Merging %s...' % sample
    chunk = []
    size  = 0
    for fname in glob('/'.join([hdfs_path,sample,'*.root'])):
        isize = size + os.path.getsize(fname)
        if isize*compression_ratio > approx_chunk_size:
            print "Executing chenneso %s" % chunk.__repr__()
            size  = os.path.getsize(fname)
            chunk = [fname]
        else:
            size = isize
            chunk.append(fname)

