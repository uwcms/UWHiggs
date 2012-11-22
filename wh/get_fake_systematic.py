#!/usr/bin/env python

'''

Get fake systematic.

'''

from rootpy import io
import sys

the_filename = sys.argv[1]
base_dir = sys.argv[2]
syst_name = sys.argv[3]

the_file = io.open(the_filename, 'READ')
nom = the_file.Get(base_dir + '/fakes').Integral()
q = the_file.Get(base_dir + '_q/fakes').Integral()

print "%s fakes %s %0.2f" % (base_dir, syst_name, 1 + (nom - q) / nom)
