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

plotter = BasePlotter(blind_region)

signs = ['os', 'ss']
jets = ['0', '1', '2']
processtype = ['gg']
threshold = ['ept30']

histo_info = [('tPt', 'p_{T}(#tau) (GeV)', 5), ('tEta', '#eta(#tau)', 2),  ('tPhi', '#phi(#tau)', 5), 
             ('ePt', 'p_{T}(e) (GeV)', 5), ('eEta', '#eta(e)', 2),  ('ePhi', '#phi(e)', 5), 
             ('et_DeltaPhi', 'e#tau #Delta#phi', 1), ('et_DeltaR', 'e#tau #Delta R', 1),
             ('h_collmass_pfmet', 'M_{coll}(e#tau) (GeV)', 1), ('h_vismass', 'M_{vis} (GeV)', 1),
             ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1) , ('ePFMET_Mt', 'M_{T} e-PFMET', 5), 
             ('tPFMET_Mt', 'M_{T} #tau-PFMET', 5) ]

logging.debug("Starting plotting")
for sign, proc, thr, njet in itertools.product(signs, processtype, threshold, jets):
        path = os.path.join(sign, proc, thr, njet)
        plotter.set_subdir(path)
        
        for var, xlabel, rebin in histo_info:
                logging.debug("Plotting %s/%s" % (path, var) )
                plotter.pad.SetLogy(False)
                #plotter.plot_without_uncert(foldername,h[0], rebin=int(h[2]), xaxis=h[1], leftside=False, show_ratio=True, ratio_range=0.5, sort=True, obj=['e'])
                plotter.plot_with_bkg_uncert(path, var, rebin, xlabel,
                                             leftside=False, show_ratio=True, ratio_range=0.5, 
                                             sort=True, obj=['e'])
                
                plotter.save(var)


