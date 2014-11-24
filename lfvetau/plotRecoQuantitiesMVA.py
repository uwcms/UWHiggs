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

embedded = True

plotter = BasePlotter(blind_region,use_embedded=embedded)

signs = ['os','ss']
jets = ['0','1','2']
processtype = ['gg']
threshold = ['ept30']

histo_info = [
   ('tPt', 'p_{T}(#tau) (GeV)', 1), 
   ('tEta', '#eta(#tau)', 1),  
   ('tPhi', '#phi(#tau)', 1), 
   ('ePt', 'p_{T}(e) (GeV)', 1), 
   ('eEta', '#eta(e)', 1),  
   ('ePhi', '#phi(e)', 1), 
   ('e_t_DPhi', 'e#tau #Delta#phi', 1), 
   ('e_t_DR', 'e#tau #Delta R', 1),
   ('h_collmass_pfmet', 'M_{coll}(e#tau) (GeV)', 1), 
   ('e_t_Mass', 'M_{vis} (GeV)', 1),
   ('jetVeto30', 'number of jets (p_{T} > 30 GeV)', 1) , 
   ('eMtToPfMet', 'M_{T} e-PFMET', 1), 
   ('tMtToPfMet', 'M_{T} #tau-PFMET', 1) , 
   ('pfMet_Et', 'pfMet', 1)
]

logging.debug("Starting plotting")
for sign, proc, thr, njet in itertools.product(signs, processtype, threshold, jets):
   path = os.path.join(sign, proc, thr, njet)
   
   plotter.set_subdir(os.path.join('embedded',path)) if embedded else plotter.set_subdir(path)
           
   for var, xlabel, rebin in histo_info:
      logging.debug("Plotting %s/%s" % (path, var) )
      plotter.pad.SetLogy(False)
      #plotter.plot_without_uncert(foldername,h[0], rebin=int(h[2]), xaxis=h[1], leftside=False, show_ratio=True, ratio_range=0.5, sort=True, obj=['e'])
      if int(njet)==2: 
         if 'collmass' in var or 'Mass' in var: 
            rebin=rebin*5
         elif not 'Eta' in var and not 'jet' in var: 
            rebin = rebin*2

      plotter.plot_with_bkg_uncert(path, var, rebin, xlabel,
                                    leftside=False, show_ratio=True, ratio_range=1., 
                                    sort=True, obj=['e'])


      plotter.save(var)

   plotter.set_subdir(os.path.join('embedded', path+'/selected'))if embedded else plotter.set_subdir(path+'/selected')

   for var, xlabel, rebin in histo_info:
      if int(njet)==1: 
         if not 'Eta' in var and not 'jet' in var: rebin = rebin*2
      if int(njet) ==2: 
         if 'collmass' in var or 'Mass' in var: rebin=rebin*5
         if 'Pt' in var or 'Mt' in var or 'pfMet' in var : rebin=rebin*4
      logging.debug("Plotting %s/%s" % (path, var) )
      plotter.pad.SetLogy(False)
      plotter.plot_with_bkg_uncert(path+'/selected', var, rebin, xlabel,
                                   leftside=False, show_ratio=True, ratio_range=1., 
                                   sort=True, obj=['e'])
      
      
      plotter.save(var,dotroot=False)

#make shapes for limit setting
signal_region = 'os/gg/ept30/%s/selected'
jets_names = [
        ('0', 'gg0etau'  ),
        ('1', 'boostetau'),
        ('2', 'vbfetau'  ),
]
pjoin = os.path.join
for njets, cat_name in jets_names:
   output_path = plotter.base_out_dir
   tfile = ROOT.TFile(pjoin(output_path, 'shapes.%s.root' % njets), 'recreate')
   output_dir = tfile.mkdir(cat_name)
   unc_conf_lines, unc_vals_lines = plotter.write_shapes( signal_region % njets, 'h_collmass_pfmet', output_dir)
   logging.warning('shape file %s created' % tfile.GetName()) 
   tfile.Close()
   with open(pjoin(output_path, 'unc.%s.conf' % njets), 'w') as conf:
      conf.write('\n'.join(unc_conf_lines))
   with open(pjoin(output_path, 'unc.%s.vals' % njets), 'w') as vals:
      vals.write('\n'.join(unc_vals_lines))

with open(pjoin(output_path,'.shapes_timestamp'),'w') as stamp:
   stamp.write('no use')



