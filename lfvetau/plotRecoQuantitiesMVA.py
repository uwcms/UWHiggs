
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
print 'blind?', blind
blind_region=[100, 150] if blind else None
#blind_region=[100, 200] if blind else None

embedded = True

plotter = BasePlotter(blind_region,use_embedded=embedded)
if not args.no_plots:
   signs = ['os']
   jets = ['0',
      '1',
      '2'
   ]
   processtype = ['gg']
   threshold = ['ept30']
   
   histo_info = [
      ('h_collmass_pfmet', 'M_{coll}(e#tau) (GeV)', 1)##,
      ##('tPt', 'p_{T}(#tau) (GeV)', 1), 
      ##('tEta', '#eta(#tau)', 1),  
      ##('tPhi', '#phi(#tau)', 1), 
      ##('ePt', 'p_{T}(e) (GeV)', 1), 
      ##('eEta', '#eta(e)', 1),  
      ##('ePhi', '#phi(e)', 1), 
      ##('e_t_DPhi', 'e#tau #Delta#phi', 1), 
      ##('e_t_DR', 'e#tau #Delta R', 1),
      ##('e_t_Mass', 'M_{vis} (GeV)', 1),
      ##('jetVeto30', 'number of jets (p_{T} > 30 GeV)', 1) , 
      ##('eMtToPfMet', 'M_{T} e-PFMET', 1), 
      ##('tMtToPfMet', 'M_{T} #tau-PFMET', 1) , 
      ##('type1_pfMet_Et', 'pfMet', 1),
      ##('vbfMass', 'M(j_{1},j_{2}) (GeV)', 1),
      ##('vbfDeta', '#Delta#eta (j_{1}, j_{2})', 1)
      
   ]
   
   logging.debug("Starting plotting")
   for sign, proc, thr, njet in itertools.product(signs, processtype, threshold, jets):
      path = os.path.join(sign, proc, thr, njet)
      
      plotter.set_subdir(os.path.join('embedded',path)) if embedded else plotter.set_subdir(path)
      
      for var, xlabel, rebin in histo_info:
         logging.debug("Plotting %s/%s" % (path, var) )
         plotter.pad.SetLogy(False)
         ## if int(njet)==2: 
         ##   if 'collmass' in var or 'Mass' in var: 
         ##      rebin=rebin
         ##   elif not 'Eta' in var and not 'jet' in var: 
         ##       rebin = rebin*2
         plotter.plot_with_bkg_uncert(path, var, rebin, xlabel,
                                      leftside=False, show_ratio=True, ratio_range=1., 
                                      sort=True, obj=['e'],  plot_data=True)
         
            
         plotter.save(var,dotroot=True)
         
      plotter.set_subdir(os.path.join('embedded', path+'/selected'))if embedded else plotter.set_subdir(path+'/selected')

      for var, xlabel, rebin in histo_info:
         ##if int(njet)==1: 
         ##   if not 'Eta' in var and not 'jet' in var: rebin = rebin
         ##if int(njet) ==2: 
         ##   if 'collmass' in var or 'Mass' in var: rebin=rebin
         ##   if 'Pt' in var or 'Mt' in var or 'pfMet' in var : rebin=rebin*4         
            
         logging.debug("Plotting %s/%s" % (path, var) )
         plotter.pad.SetLogy(False)
         plotter.plot_with_bkg_uncert(path+'/selected', var, rebin, xlabel,
                                      leftside=False, show_ratio=True, ratio_range=1., 
                                      sort=True, obj=['e'], plot_data=True)
      

         plotter.save(var,dotroot=True)

#make shapes for limit setting
if not args.no_shapes:
   signal_region = 'os/gg/ept30/%s/selected'
   ##signal_region = 'os/gg/ept30/%s'
   jets_names = [
           ('0', 'gg0etau'  , 1),
           ('1', 'boostetau', 1),#was 2
           ('2', 'vbfetau'  , 1),#was 5
   ]
   pjoin = os.path.join
   for njets, cat_name, rebin in jets_names:
      output_path = plotter.base_out_dir
      tfile = ROOT.TFile(pjoin(output_path, 'shapes.%s.root' % njets), 'recreate')
      output_dir = tfile.mkdir(cat_name)
      unc_conf_lines, unc_vals_lines = plotter.write_shapes( 
         signal_region % njets, 'h_collmass_pfmet', output_dir, rebin=rebin,
         br_strenght=1, last=300)
      logging.warning('shape file %s created' % tfile.GetName()) 
      tfile.Close()
      with open(pjoin(output_path, 'unc.%s.conf' % njets), 'w') as conf:
         conf.write('\n'.join(unc_conf_lines))
      with open(pjoin(output_path, 'unc.%s.vals' % njets), 'w') as vals:
         vals.write('\n'.join(unc_vals_lines))

   with open(pjoin(output_path,'.shapes_timestamp'),'w') as stamp:
      stamp.write('no use')



