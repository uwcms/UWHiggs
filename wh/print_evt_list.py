#! /bin/env python

import ROOT
import sys

tfile_path = sys.argv[-2]
tree_path  = sys.argv[-1]

tfile = ROOT.TFile.Open(tfile_path)
tree  = tfile.Get(tree_path)

for row in tree:
    print '%i:%i:%i' % (row.run, row.lumi, row.evt)


