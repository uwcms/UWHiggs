#from mauro plotters

#Set logging before anything to override rootpy very verbose defaults
import sys
import logging
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

import os
import ROOT
from pdb import set_trace
from FinalStateAnalysis.PlotTools.MegaBase import make_dirs
from FinalStateAnalysis.MetaData.data_styles import data_styles
from FinalStateAnalysis.PlotTools.BlindView import BlindView,  blind_in_range
from FinalStateAnalysis.PlotTools.SubtractionView      import SubtractionView, PositiveView
import itertools
import glob
import sys
from BasePlotter import BasePlotter
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument('--no-plots', dest='no_plots', action='store_true',
                    default=False, help='Does not print plots')
parser.add_argument('--no-shapes', dest='no_shapes', action='store_true',
                    default=False, help='Does not create shapes for limit computation')
args = parser.parse_args()

jobid = os.environ['jobid']
#jobid = 'MCntuples_3March' 
channel = 'et'
import rootpy.plotting.views as views
        
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

print "\nPlotting %s for %s\n" % (channel, jobid)

#check if blind
blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
blind_region=[100, 150] if blind else None

embedded = True

plotter = BasePlotter(blind_region,use_embedded=embedded)

signal_region = 'os/gg/ept30/%s/selected'
jets_names = [
   ('0', 'gg0etau'  , 1),
   ('1', 'boostetau', 1),
   ('2', 'vbfetau'  , 5),
]
pjoin = os.path.join
for njets, cat_name, rebin in jets_names:
   output_path = plotter.base_out_dir
   tfile = ROOT.TFile(pjoin(output_path, 'shapes.%s.root' % njets), 'recreate')
   output_dir = tfile.mkdir(cat_name)
   unc_conf_lines, unc_vals_lines = plotter.write_shapes_for_yields( 
      signal_region % njets, 'h_collmass_pfmet', output_dir, rebin=rebin)
   logging.warning('shape file %s created' % tfile.GetName()) 
   tfile.Close()
   with open(pjoin(output_path, 'unc.%s.conf' % njets), 'w') as conf:
      conf.write('\n'.join(unc_conf_lines))
   with open(pjoin(output_path, 'unc.%s.vals' % njets), 'w') as vals:
      vals.write('\n'.join(unc_vals_lines))
         
with open(pjoin(output_path,'.shapes_timestamp'),'w') as stamp:
   stamp.write('no use')



