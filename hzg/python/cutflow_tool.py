#!/usr/bin/env python

import cPickle as pickle
import gzip
import sys

if len(sys.argv) != 3:
    sys.exit('this program accepts up to two arguments:\n'\
             '\targ 1: the dataset type: mm, ee, mmg, eeg and\n'\
             '\targ 2: a file to compare to (with run lumi event in each line of the file)!')

def load(filename):    
    file = gzip.GzipFile(filename, 'rb')
    buffer = ""
    while 1:
        data = file.read()
        if data == "":
            break
        buffer += data
    object = pickle.loads(buffer)
    file.close()
    return object

def load_compare(filename):
    file = open(filename, 'rb')
    event_list = []
    for line in file.readline():
        temp = line.split()
        run = int(temp[0])
        lumi = int(temp[1])
        event = long(temp[2])
        event_list.append((run,lumi,event))

mmg_cuts = {}

def do_cutflow_compare(events,compare):
    for event in compare:
        if event not in events:
            print event,'is not in raw sample, somehow you have selected an event with min(DR all objects) < 0.3'
        else:
            pass

options = {'mmg':['mmg_event_dump.pkl.gz',mmg_cuts]}

if __name__ == "__main__":
    opt = sys.argv[1]
    
    events = load(options[opt][0])
    options[opt][1](events)
    
    
