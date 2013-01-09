import os
from os import path
import subprocess

#this class expects to be given the directory which contains
#the root of a directory structure containing all signal samples,
#background samples, and real data samples
#
# like: /path/to/dirs/<sample dir name>/<channels>/<sample groups>/<samples>
#
def get_leaf_dirs(arg,dirname,names):
    for i,name in enumerate(names):        
        if not path.isdir(path.join(dirname,name)):
            del names[i]
    if len(names) == 0:
        arg.append(dirname)           

class directory_prep:
    def __init__(self, path_to_root, sample_name):
        self._path_to_root = path_to_root
        self._sample_name = sample_name    

    def build_inputs(self):
        root_path = path.join(self._path_to_root,self._sample_name)
        root_path = path.expanduser(root_path)
        
        sample_dirs = []
        path.walk(root_path,get_leaf_dirs,sample_dirs)

        proc_groups = {}

        for dir in sample_dirs:
            dirroot,name = path.split(dir)
            dontcare,procname = path.split(dirroot)

            if procname not in proc_groups:
                proc_groups[procname] = []
            
            print 'Creating input file %s.root'%name
            command = 'hadd -v 0 -f %s.root %s/*.root'%(dir,dir)            
            subprocess.call(command,shell=True)

            proc_groups[procname].append('%s.root'%dir)

        return proc_groups

            
