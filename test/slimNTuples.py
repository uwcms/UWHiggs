#! /bin/env python

import os
import sys
import subprocess
import re
from glob import glob
from fnmatch import fnmatch
from hashlib import md5
import time
from progressbar import ETA, ProgressBar, FormatLabel, Bar

if len(sys.argv) < 2 or '-h' in sys.argv or '--help' in sys.argv:
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
make_dir(listsdir)

#define matching function, use generator to be faster (should stop at the first match)
match  = lambda x:  any( ( fnmatch(x,y) for y in usedBranches) )
with open(listsdir+'/trees_location.list','w') as f:
    f.write('\n'.join(forest))

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
        with open(listsdir+'/%s.list' % tree_name.replace('/','_'),'w') as f:
            f.write('\n'.join(matching_branches))

compression_ratio = float(kept_branches) / float(tot_branches)
#os.path.getsize(
print 'Compression ratio: %s' % compression_ratio
print "Beginning multi threaded copy/complession/merge..."
approx_chunk_size = 70*10**6 #~50MB
procs   = {}
threads = os.environ['megaworkers'] if 'megaworkers' in os.environ else 2
with open("copyMergeSlim_"+timestamp+".status",'w') as status_file:
    logdir = os.getcwd()+'/'+timestamp
    make_dir(logdir)
    make_dir('/'.join(['','scratch',os.environ['USER'],'data',JOBID]) )
    files_dict   = dict( [(sample, glob('/'.join([hdfs_path,sample,'*.root'])))  for sample in SAMPLES] )
    tot_numfiles = reduce(lambda x, y: x+y, map(len,files_dict.itervalues()) )
    pbar  = ProgressBar(
        widgets = [
            FormatLabel(
                'Copied %(value)i/' + str(tot_numfiles) + ' files. '),
                ETA(),
                Bar('>')],
        maxval = tot_numfiles ).start()
    files_already_merged = 0
    for sample in SAMPLES:
        print 'Merging %s...' % sample
        sample_dir = '/'.join(['','scratch',os.environ['USER'],'data',JOBID,sample])
        make_dir( sample_dir )
        chunk = []
        size  = 0
        for fname in files_dict[sample]:
            isize = size + os.path.getsize(fname)
            if isize*compression_ratio > approx_chunk_size:
                hasher  = md5()
                hasher.update(' '.join(chunk))
                outName = hasher.hexdigest()
                errlog = open('/'.join([logdir,outName+'.err']), 'w')
                outlog = open('/'.join([logdir,outName+'.out']), 'w')
                command = ['slimAndMergeNtuple',listsdir,'/'.join([sample_dir,outName+'.root']) ]
                command.extend(chunk)
                outlog.write("Executed: \n")
                outlog.write(' '.join(command))
                outlog.write('\n\n')
                proc   = subprocess.Popen(command, stdout=outlog, stderr=errlog)
                procs[outName] = {
                    'proc' : proc,
                    'out'  : outlog,
                    'err'  : errlog,
                    'nfile': len(chunk),
                    }
                size  = os.path.getsize(fname)
                chunk = [fname]
            else:
                size = isize
                chunk.append(fname)
            while len(procs) >= threads: #wait the jobs to finish
                time.sleep(3) #wait
                todel = []
                for name, info in procs.iteritems():
                    info['err'].flush()
                    info['out'].flush()
                    ret_code = info['proc'].poll()
                    if ret_code is not None: #One is done!
                        val = 'OK' if ret_code == 0 else 'ERROR'
                        status_file.write('%s MERGE_%s\n' % (name, val) ) #write on the status if it's done
                        info['err'].close()
                        info['out'].close()
                        files_already_merged += info['nfile']
                        pbar.update(files_already_merged)
                        todel.append(name)
                for name in todel:
                    del procs[name]

    for name, info in procs.iteritems():
        ret_code = info['proc'].communicate() #now wait everythin is done
        if ret_code is not None: #One is done!
            print "%s DONE!" % name
            val = 'OK' if ret_code == 0 else 'ERROR'
            status_file.write('%s MERGE_%s\n' % (name, val) ) #write on the status if it's done
            info['err'].close()
            info['out'].close()
            files_already_merged += info['nfile']
            pbar.update(files_already_merged)
            todel.append(name)
    pbar.finish()
