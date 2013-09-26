#!/usr/bin/env python

'''

Get fake systematic.

'''

from rootpy import io
import rootpy.plotting.views as views
import sys

the_filename = sys.argv[1]
base_dirs = sys.argv[2].split(',')
syst_name = sys.argv[3]

the_file = io.open(the_filename, 'READ')
nom_view = views.SumView(
    *[ views.SubdirectoryView( the_file, dirname ) for dirname in base_dirs]
)

qcd_view = views.SumView(
    *[ views.SubdirectoryView( the_file, dirname+'_q' ) for dirname in base_dirs]
)

nom = nom_view.Get('fakes').Integral()
qcd = qcd_view.Get('fakes').Integral()

print "%s fakes %s %0.2f" % (','.join(base_dirs), syst_name, 1 + abs(nom - qcd) / nom)
