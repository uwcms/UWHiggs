#! /bin/env python

import sys
import os
import re

"Usage pull files get_unc_val_from_pull.py pulls.txt nuisance initial_value header"

header   = sys.argv[-1]
prefit   = float(sys.argv[-2])
nuisance = sys.argv[-3]
filename = sys.argv[-4]

if not os.path.isfile(filename):
    print "%s not found!" % filename
    sys.exit(1)

lines = open(filename).readlines()
parser    = re.compile(r'(?P<nuisance>\w+)\s+(?P<pull>(?:[\*\!] )?[\-+]\d+\.\d+)\,\s(?P<constraint>\d+\.\d+)')

for line in lines:
    match = parser.match(line)
    if match:
        if match.group('nuisance') != nuisance:
            continue
        pull       = float(match.group('pull').replace('* ','').replace('! ',''))
        constraint = float(match.group('constraint'))
        max_dev    = max(abs(pull), abs(constraint))
        new_dev    = abs(1-prefit)*max_dev + 1
        print '%s %s %.2f' % (header, nuisance, new_dev)
        
