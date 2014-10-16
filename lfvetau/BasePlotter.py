'''

Base class which makes nice plots.

Author: Evan K. Friis, UW

'''

import fnmatch
import re
import os
import math
import rootpy.plotting.views as views
from rootpy.plotting.hist import HistStack
import rootpy.plotting as plotting
from FinalStateAnalysis.MetaData.data_views import data_views
from FinalStateAnalysis.PlotTools.RebinView import RebinView
from FinalStateAnalysis.PlotTools.BlindView import BlindView,  blind_in_range
from FinalStateAnalysis.Utilities.struct import struct
import FinalStateAnalysis.Utilities.prettyjson as prettyjson
from FinalStateAnalysis.MetaData.data_styles import data_styles
from FinalStateAnalysis.PlotTools.Plotter  import Plotter
from FinalStateAnalysis.PlotTools.SubtractionView      import SubtractionView, PositiveView

import ROOT
import glob
from pdb import set_trace

def create_mapper(mapping):
    def _f(path):
        for key, out in mapping.iteritems():
            if key == path:
                path = path.replace(key,out)
                print 'path', path
        return path
    return _f


class BasePlotter(Plotter):
    def __init__ (self, channel, files, lumifiles, outputdir,forceLumi=-1): 
        cwd = os.getcwd()
        self.channel = channel
        jobid = os.environ['jobid']
        period = '7TeV' if '7TeV' in jobid else '8TeV'
        self.period = period
        self.sqrts = 7 if '7TeV' in jobid else 8
        self.mc_samples = [
            'ggHiggsToETau',
            'vbfHiggsToETau',
            'GluGluToHToTauTau_M-125_8TeV-powheg-pythia6',
            'VBF_HToTauTau_M-125_8TeV-powheg-pythia6',
            #'Zjets_M50', 
            'Zjets_M50_skimmedLL',
            'Z1jets_M50_skimmedLL',
            'Z2jets_M50_skimmedLL',
            'Z3jets_M50_skimmedLL',
            'Z4jets_M50_skimmedLL',
            'Zjets_M50_skimmedTT',
            'Z1jets_M50_skimmedTT',
            'Z2jets_M50_skimmedTT',
            'Z3jets_M50_skimmedTT',
            'Z4jets_M50_skimmedTT',
            'TTJets*',
            'T_t*',
            'Tbar_t*', 
            'WplusJets_madgraph_skimmed',
            'Wplus1Jets_madgraph',
            # 'Wplus1Jets_madgraph_tapas',
            'Wplus2Jets_madgraph',
            # 'Wplus2Jets_madgraph_tapas',
            'Wplus3Jets_madgraph',
            'Wplus4Jets_madgraph',
            'WWJets*',
            'WZJets*',
            'ZZJets*',
            'Fake*',
            'data_Zetau*',
            'data_Single*'
            
        ]
        
##        files = []
##        lumifiles = []
##  
##        for x in mc_samples:
##            files.extend(glob.glob('results/%s/LFVHETauAnalyzerMVA/%s.root' % (jobid, x)))
##            lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))
##        self.outputdir = 'plots/%s/LFVHETauAnalyzerMVA/%s' % (jobid, channel)
##        self.base_out_dir = self.outputdir
##        if not os.path.exists(self.outputdir):
##            os.makedirs(self.outputdir)
##
        self.blindstart=100
        self.blindend=150

        blinder = None
        blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
        print '\n\nRunning Blind: %s\n\n' % blind
        self.blind = blind

        if blind:
            # Don't look at the SS all pass region
            blinder = lambda x: BlindView(x, "os/.*",blind_in_range(100, 150))
        super(BasePlotter, self).__init__(files, lumifiles, outputdir,   blinder)
    def simpleplot_mc(self, folder, signal, variable, rebin=1, xaxis='',
                      leftside=True, xrange=None, preprocess=None, sort=True,forceLumi=-1):
        ''' Compare Monte Carlo signal to bkg '''
        #path = os.path.join(folder, variable)
        #is_not_signal = lambda x: x is not signame
        #set_trace()
        
        mc_stack_view = self.make_stack(rebin, preprocess, folder)
        mc_stack=mc_stack_view.Get(variable)
        
        #self.canvas.SetLogy(False)
        mc_stack.SetTitle('')
        mc_stack.Draw()
        
        mc_stack.GetHistogram().GetXaxis().SetTitle(xaxis)

        if xrange:
            mc_stack.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            mc_stack.Draw()

        signalview=[]
        mymax=0
        for sig in signal:
                signal_view=self.get_view(sig)
                if preprocess:
                    signal_view=preprocess(signal_view)
                signal_view=self.get_wild_dir(
                    self.rebin_view(signal_view,rebin),folder)
                signal=signal_view.Get(variable)
                signalview.append(signal)
                signal.Draw("SAME")
                if signal.GetBinContent(signal.GetMaximumBin()) > mymax:
                    mymax = signal.GetBinContent(signal.GetMaximumBin())
                    
                self.keep.append(signal)
        if mymax > mc_stack.GetMaximum():
            mc_stack.SetMaximum(mymax*1.2)

        self.keep.append(mc_stack)
        
        all_var=[]
        all_var.extend([mc_stack]) 
        all_var.extend(signalview) 
            
        self.add_legend(all_var, leftside, entries=len(mc_stack)+len(signalview))
 
    def plot_data(self, folder, variable, rebin=1, xaxis='',
                  leftside=True, xrange=None, preprocess=None):
        ''' Compare Monte Carlo signal to bkg '''
        #path = os.path.join(folder, variable)
        #is_not_signal = lambda x: x is not signame
        #set_trace()
        data_view = self.get_view('data')
        if preprocess:
            data_view = preprocess( data_view )
        data_view = self.get_wild_dir(
            self.rebin_view(data_view, rebin),
            folder
            )
        data = data_view.Get(variable)
        data.Draw()

         
        data.GetXaxis().SetTitle(xaxis)

        if xrange:
            data.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            data.Draw()
        self.keep.append(data)
        
        
        
    def compare_data(self, folder, variable, folder2, rebin=1, xaxis='',
                     leftside=True, xrange=None, preprocess=None, show_ratio=False, ratio_range=0.2,rescale=1. ):
        ''' Compare Monte Carlo signal to bkg '''
        #path = os.path.join(folder, variable)
        #is_not_signal = lambda x: x is not signame
        #set_trace()
        data_view = self.get_view('data')
        data_view2 = self.get_view('data')
        if preprocess:
            data_view = preprocess( data_view )
        data_view = self.get_wild_dir(
            self.rebin_view(data_view, rebin),
            folder
            )
        data_view2 = self.get_wild_dir(
            self.rebin_view(data_view2, rebin),
            folder2
            )
        data = data_view.Get(variable)
        data2 = data_view2.Get(variable)
         
        data.GetXaxis().SetTitle(xaxis)
        #if (data.Integral()!=0 and data2.Integral()!=0) :
        #    data.Scale(1./data.Integral())
        #    data2.Scale(1./data2.Integral())
        data.Draw()
        data.GetYaxis().SetRangeUser(0, data.GetBinContent(data.GetMaximumBin())*1.2)
        if xrange:
            data.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            data.Draw()

        data2.Draw("SAME")
        data2.SetMarkerColor(2)
        data2.Draw("SAME")
       

        if show_ratio:
            self.add_sn_ratio_plot(data2, data, xrange, ratio_range)

        self.keep.append(data)
        self.keep.append(data2)
        
        
    def add_ratio_bandplot(self, data_hist, mc_stack, err_hist,  x_range=None, ratio_range=0.2):
        #resize the canvas and the pad to fit the second pad
        self.canvas.SetCanvasSize( self.canvas.GetWw(), int(self.canvas.GetWh()*1.3) )
        self.canvas.cd()
        self.pad.SetPad(0, 0.33, 1., 1.)
        self.pad.Draw()
        self.canvas.cd()
        #create lower pad
        self.lower_pad = plotting.Pad('low', 'low', 0, 0., 1., 0.33)
        self.lower_pad.Draw()
        self.lower_pad.cd()

        
        mc_hist    = None
        if isinstance(mc_stack, plotting.HistStack):
            mc_hist = sum(mc_stack.GetHists())
        else:
            mc_hist = mc_stack
        data_clone = data_hist.Clone()
        data_clone.Divide(mc_hist)
        
        band = err_hist.Clone()
        
        err = []
        ibin =1 
        while ibin < band.GetXaxis().GetNbins()+1:
            if mc_hist.GetBinContent(ibin) <> 0 : 
                err.append((ibin, band.GetBinError(ibin)/band.GetBinContent(ibin)))
            ibin+=1

        band.Divide(mc_hist.Clone())
        #print err
        for ibin in err:
            band.SetBinError(ibin[0], ibin[1])

        if not x_range:
            nbins = data_clone.GetNbinsX()
            x_range = (data_clone.GetBinLowEdge(1), 
                       data_clone.GetBinLowEdge(nbins)+data_clone.GetBinWidth(nbins))
        else:
            data_clone.GetXaxis().SetRangeUser(*x_range)


        ref_function = ROOT.TF1('f', "1.", *x_range)
        ref_function.SetLineWidth(2)
        ref_function.SetLineStyle(2)
        
        data_clone.Draw()
 
        if ratio_range:
            data_clone.GetYaxis().SetRangeUser(1-ratio_range, 1+ratio_range)
        ref_function.Draw('same')
        band.SetMarkerStyle(0)
        band.SetLineColor(1)
        band.SetFillStyle(3001)
        band.SetFillColor(1)

        band.Draw('samee2')
       
        self.keep.append(data_clone)
        self.keep.append(band)
        self.keep.append(ref_function)
        self.pad.cd()
        return data_clone

 

    def get_isoid_unc (self, folder, variable, rebin, preprocess, obj=['e1', 'e2']):
       # 'e1idp1s/','e1idm1s/',  'e1isop1s/','e1isom1s/','e2idp1s/','e2idm1s/',  'e2isop1s/','e2isom1s/'
        unc_list = []
        for n, lep in enumerate(obj):
            
            bkg_lep_iso_p1s_stack_view = self.make_stack(rebin, preprocess, lep+'isop1s/'+folder)
            bkg_lep_iso_p1s_stack= bkg_lep_iso_p1s_stack_view.Get(variable)
            histo_lep_iso_p1s =  bkg_lep_iso_p1s_stack.GetStack().Last().Clone()
            
            bkg_lep_iso_m1s_stack_view = self.make_stack(rebin, preprocess, lep+'isom1s/'+folder)
            bkg_lep_iso_m1s_stack=bkg_lep_iso_m1s_stack_view.Get(variable)
            histo_lep_iso_m1s = bkg_lep_iso_m1s_stack.GetStack().Last().Clone()

            bkg_lep_id_p1s_stack_view = self.make_stack(rebin, preprocess, lep+'idp1s/'+folder)
            bkg_lep_id_p1s_stack= bkg_lep_id_p1s_stack_view.Get(variable)
            histo_lep_id_p1s =  bkg_lep_id_p1s_stack.GetStack().Last().Clone()
        
            bkg_lep_id_m1s_stack_view = self.make_stack(rebin, preprocess, lep+'idm1s/'+folder)
            bkg_lep_id_m1s_stack=bkg_lep_id_m1s_stack_view.Get(variable)
            histo_lep_id_m1s = bkg_lep_id_m1s_stack.GetStack().Last().Clone()
        
            #bkg_obj2_iso_p1s_stack_view = self.make_stack(rebin, preprocess, obj2+'isop1s/'+folder)
            #bkg_obj2_iso_p1s_stack= bkg_obj2_iso_p1s_stack_view.Get(variable)
            #histo_obj2_iso_p1s =  bkg_obj2_iso_p1s_stack.GetStack().Last().Clone()
        
            #bkg_obj2_iso_m1s_stack_view = self.make_stack(rebin, preprocess, obj2+'isom1s/'+folder)
            #bkg_obj2_iso_m1s_stack=bkg_obj2_iso_m1s_stack_view.Get(variable)
            #histo_obj2_iso_m1s = bkg_obj2_iso_m1s_stack.GetStack().Last().Clone()

            #bkg_obj2_id_p1s_stack_view = self.make_stack(rebin, preprocess, obj2+'idp1s/'+folder)
            #bkg_obj2_id_p1s_stack= bkg_obj2_id_p1s_stack_view.Get(variable)
            #histo_obj2_id_p1s =  bkg_obj2_id_p1s_stack.GetStack().Last().Clone()
        
            #bkg_obj2_id_m1s_stack_view = self.make_stack(rebin, preprocess, obj2+'idm1s/'+folder)
            #bkg_obj2_id_m1s_stack=bkg_obj2_id_m1s_stack_view.Get(variable)
            #histo_obj2_id_m1s = bkg_obj2_id_m1s_stack.GetStack().Last().Clone()
        
            self.keep.append(histo_lep_iso_p1s)
            self.keep.append(histo_lep_iso_m1s)
            self.keep.append(histo_lep_id_p1s)
            self.keep.append(histo_lep_id_m1s)
            #self.keep.append(histo_obj2_iso_p1s)
            #self.keep.append(histo_obj2_iso_m1s)
            #self.keep.append(histo_obj2_id_p1s)
            #self.keep.append(histo_obj2_id_m1s)

            unc_list.extend([(histo_lep_iso_p1s, histo_lep_iso_m1s), (histo_lep_id_p1s, histo_lep_id_m1s)])
        return unc_list
#[(histo_obj1_iso_p1s, histo_obj1_iso_m1s), (histo_obj1_id_p1s, histo_obj1_id_m1s), (histo_obj2_iso_p1s, histo_obj2_iso_m1s), (histo_obj2_id_p1s, histo_obj2_id_m1s)]
    def add_ratio_diff(self, data_hist, mc_stack, err_hist,  x_range=None, ratio_range=0.2):
        #resize the canvas and the pad to fit the second pad
        self.canvas.SetCanvasSize( self.canvas.GetWw(), int(self.canvas.GetWh()*1.3) )
        self.canvas.cd()
        self.pad.SetPad(0, 0.33, 1., 1.)
        self.pad.Draw()
        self.canvas.cd()
        #create lower pad
        self.lower_pad = plotting.Pad('low', 'low', 0, 0., 1., 0.33)
        self.lower_pad.Draw()
        self.lower_pad.cd()

        
        mc_hist    = None
        if isinstance(mc_stack, plotting.HistStack):
            mc_hist = sum(mc_stack.GetHists())
        else:
            mc_hist = mc_stack
        data_clone = data_hist.Clone()
        data_clone.Sumw2()
        data_clone.Add(mc_hist, -1)
        data_clone.Divide(mc_hist)
        
        band = err_hist.Clone()
        
        err = []
        ibin =1 
        while ibin < band.GetXaxis().GetNbins()+1:
            if mc_hist.GetBinContent(ibin) <> 0 : 
                err.append((ibin, band.GetBinError(ibin)/band.GetBinContent(ibin)))
                
            ibin+=1

        band.Divide(mc_hist.Clone())
        #print err
        for ibin in err:
            band.SetBinError(ibin[0], ibin[1])
        blind=None
        blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
        if blind<>None:
            for ibin in err:
                if ibin>=band.FindBin(self.blindstart) and ibin <= band.FindBin(self.blindend):
                    band.SetBinError(ibin, 10)

        if not x_range:
            nbins = data_clone.GetNbinsX()
            x_range = (data_clone.GetBinLowEdge(1), 
                       data_clone.GetBinLowEdge(nbins)+data_clone.GetBinWidth(nbins))
        else:
            data_clone.GetXaxis().SetRangeUser(*x_range)


        ref_function = ROOT.TF1('f', "0.", *x_range)
        ref_function.SetLineWidth(2)
        ref_function.SetLineStyle(2)
        
        data_clone.Draw()
 
        if ratio_range:
            data_clone.GetYaxis().SetRangeUser(-ratio_range, +ratio_range)
        ref_function.Draw('same')
        band.SetMarkerStyle(0)
        band.SetLineColor(1)
        band.SetFillStyle(3001)
        band.SetFillColor(1)

        band.Draw('samee2')
       
        self.keep.append(data_clone)
        self.keep.append(band)
        self.keep.append(ref_function)
        self.pad.cd()
        return data_clone

    def plot_with_bkg_uncert (self, folder, variable, rebin=1, xaxis='',
                        leftside=True, xrange=None, preprocess=None,
                              show_ratio=False, ratio_range=0.2, sort=True, obj=['e1', 'e2']):
        
        

        mc_stack_view = self.make_stack(rebin, preprocess, folder, sort)

        mc_stack = mc_stack_view.Get(variable)
        mc_stack.Draw()
        
        self.canvas.SetLogy(True)
        mc_stack.GetHistogram().GetXaxis().SetTitle(xaxis)
        if xrange:
            mc_stack.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            mc_stack.Draw()
        self.keep.append(mc_stack)
        

        finalhisto= mc_stack.GetStack().Last().Clone()
        finalhisto.Sumw2()

        histlist = mc_stack.GetHists();
        bkg_stack = mc_stack_view.Get(variable)
        bkg_stack.GetStack().RemoveLast()## mettere il check se c'e` il fake altirmenti histo=mc_stack.GetStack().Last().Clone()
        histo=bkg_stack.GetStack().Last().Clone()
        histo.Sumw2()
        
        fake_p1s_histo=None
      
        if not folder.startswith('tLoose'):
            ##if folder.startswith('os'):
            ## newfolder=folder#.replace('os','ss')
                ##      print 'newfolder', newfolder
                ##      #tau fake error
                ##  
                ##fake_stack_view = self.make_stack(rebin, preprocess, 'tLoose/'+newfolder, sort)
                ##      fake_p1s_stack_view = self.make_stack(rebin, preprocess, 'tLooseUp/'+newfolder, sort)
                ##      fake_m1s_stack_view = self.make_stack(rebin, preprocess, 'tLooseDown/'+newfolder, sort)
                ##                
                ##  #if (type(fake_stack_view)==ROOT.THStack):
                ##  #print 'adding tau fakerate uncertainty'
                ##  
                ##  try:
                ##fake_stack = fake_stack_view.Get(variable)
                ##fakehisto = fake_stack.GetStack().Last().Clone()
          ##      
          ##      fake_p1s_stack = fake_p1s_stack_view.Get(variable) 
          ##      fake_p1s_histo = fake_p1s_stack.GetStack().Last().Clone()
          ##      fake_m1s_stack = fake_m1s_stack_view.Get(variable) 
          ##      fake_m1s_histo = fake_m1s_stack.GetStack().Last().Clone()
          ##      
          ##  except:
          ##      print 'no fake taus'
          ##      fake_diffup = 0
          ##      fake_diffdown = 0
                
            
            print folder
            print self.mc_samples
            isFakesIn= False
            if 'Fakes' in self.mc_samples:
                self.mc_samples.remove('Fakes')  
                self.mc_samples.remove('finalDYLL')
                isFakesIn=True
            
            bkg_p1s_stack_view = self.make_stack(rebin, preprocess, 'p1s/'+folder, sort)
            bkg_p1s_stack=bkg_p1s_stack_view.Get(variable)
            histo_p1s = bkg_p1s_stack.GetStack().Last().Clone()
        
            bkg_m1s_stack_view = self.make_stack(rebin, preprocess, 'm1s/'+folder, sort)
            bkg_m1s_stack=bkg_m1s_stack_view.Get(variable)
            histo_m1s = bkg_m1s_stack.GetStack().Last().Clone()
        
            tr_p1s_stack_view = self.make_stack(rebin, preprocess, 'trp1s/'+folder, sort)
            tr_p1s_stack=tr_p1s_stack_view.Get(variable)
            histotr_p1s = tr_p1s_stack.GetStack().Last().Clone()
        
            tr_m1s_stack_view = self.make_stack(rebin, preprocess, 'trm1s/'+folder, sort)
            tr_m1s_stack=tr_m1s_stack_view.Get(variable)
            histotr_m1s = tr_m1s_stack.GetStack().Last().Clone()
            #if os.path.exists(folder+'tLoose/') :
            #    fr_tau_stack_view = self.make_stack(rebin, preprocess, 'tLoose/'+folder, sort)
            #    fr_tau_stack=fr_tau_stack_view.Get(variable)
            #    histofr_tau = fr_tau_stack.GetStack().Last().Clone()
            #    histofr_tau.SetTitle('Fakes')
            
            xsec_unc_mapper = {
                'EWK Dibosons': 0.056,
                'SM Higgs': 0.0, #still to add
                't#bar{t}' : 0.026,
                'Single Top' : 0.041,
                'DY (#rightarrow #tau#tau)  + jets' : 0., 
                'DY (#rightarrow ll)  + jets' : 0.032,
                'W + jets': 0.035,
                'QCD' : 0.0,
                'Fakes' :0.0,
            }
            
        ibin =1
            
        iso_id_unc = self.get_isoid_unc( folder, variable, rebin, preprocess, obj)
        if isFakesIn:
            self.mc_samples.append('Fakes')
            self.mc_samples.append('finalDYLL')
          
        met_histos=[]
        if 'Met' in variable or 'MET' in variable:
            if not 'mva' in variable or not 'MVA' in variable:
                if 'type1_' in variable or '_Ty1' in variable:
                    ##print variable
                    newMetName=variable.replace('type1_', '') if 'type1_' in variable  else variable.replace('_Ty', '')
                    met_histos.extend([mc_stack_view.Get(newMetName+'_jes'),mc_stack_view.Get(newMetName+'_mes'),mc_stack_view.Get(newMetName+'_tes'),mc_stack_view.Get(newMetName+'_ues'), mc_stack_view.Get(newMetName+'_ees') ]) ##add the ees
                else:
                    met_histos.extend([mc_stack_view.Get(variable+'_jes'),mc_stack_view.Get(variable+'_mes'),mc_stack_view.Get(variable+'_tes'),mc_stack_view.Get(variable+'_ues'),mc_stack_view.Get(variable+'_ees')]) ##add the ees
            
        while ibin < histo.GetXaxis().GetNbins()+1: 
            ##print ibin,histo.GetXaxis().GetNbins()
            err2=0

            #if isFakesIn:
                #diff_fake_up= abs(fakehisto.GetBinContent(ibin)-fake_p1s_histo.GetBinContent(ibin))
                #diff_fake_down= abs(fakehisto.GetBinContent(ibin)-fake_m1s_histo.GetBinContent(ibin))
                #print ibin, fakehisto.GetBinContent(ibin), fake_p1s_histo.GetBinContent(ibin), fake_m1s_histo.GetBinContent(ibin)
                #diff_fake = diff_fake_up if  diff_fake_up >  diff_fake_down else  diff_fake_down

                #err2+= pow(diff_fake,2)
                ##err2+=pow(fakehisto.GetBinContent(ibin) * 0.3)
                #print 'TAU SYSTEMATICS:', fakehisto.GetBinContent(ibin), diff_fake 

            for h in histlist:
                err2 += h.GetBinError(ibin)**2
                
                err2 += h.GetBinContent(ibin)*xsec_unc_mapper[h.GetTitle()]
              
                     
                    
                #tau fake rake, I approximate with 5%. Redo it in the correct way
                #diff_fake = abs(histo.GetBinContent(ibin) - finalhisto.GetBinContent(ibin))
                #err2 += pow(diff_fake*0.05,2)
                
                #why is not working?
                #diff_fake = abs(fake_diffup.GetBinContent(ibin)) if abs(fake_diffup.GetBinContent(ibin)) > abs(fake_diffdown.GetBinContent(ibin)) else abs(fake_diffdown.GetBinContent(ibin))
 
            if not folder.startswith('tLoose'):
                diff_p = abs(histo.GetBinContent(ibin) - histo_p1s.GetBinContent(ibin))
                diff_m = abs(histo.GetBinContent(ibin) - histo_m1s.GetBinContent(ibin))
                err2 += diff_p**2  if diff_p > diff_m  else diff_m**2  
                diff_tr_p = abs(histo.GetBinContent(ibin) - histotr_p1s.GetBinContent(ibin))
                diff_tr_m = abs(histo.GetBinContent(ibin) - histotr_m1s.GetBinContent(ibin))
                err2 += diff_tr_p**2  if diff_tr_p > diff_tr_m  else diff_tr_m**2  
                err2 += (0.026*h.GetBinContent(ibin))**2 #lumi error

 
                for h in iso_id_unc :
                    diff_isoid_p = abs(histo.GetBinContent(ibin) - h[0].GetBinContent(ibin))
                    diff_isoid_m = abs(histo.GetBinContent(ibin) - h[1].GetBinContent(ibin))
                    err2 += diff_isoid_p**2  if diff_isoid_p > diff_isoid_m  else diff_isoid_m**2  
 
                #add met unc https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Systematics_Tools only to met plots (no cut on the met)
            if 'Met' in variable  or 'MET' in variable:
                if not 'mva' in variable or not 'MVA' in variable:
                    if not '_jes' in variable and not '_mes' in variable and not '_tes' in variable \
                       and not '_ees' in variable  and not '_ues' in variable:
                        if not 'type1_' in variable and not '_Ty1' in variable:
                            for h in met_histos:
                                diff = abs(histo.GetBinContent(ibin) - h.GetStack().Last().GetBinContent(ibin))
                                ##print variable, diff, histo.GetBinContent(ibin), h.GetStack().Last().GetBinContent(ibin), h.GetName()
                                err2+=diff**2
                        else: 
                            myh = mc_stack_view.Get(variable.replace('type1_', '')) if 'type1_' in variable else variable.replace('_Ty', '')
                            for h in met_histos:
                                diff = abs(h.GetStack().Last().GetBinContent(ibin) - myh.GetStack().Last().GetBinContent(ibin))
                                err2+=diff**2
            
            #histo.SetBinError(ibin, math.sqrt(err2))
            #finalhisto.SetBinError(ibin, math.sqrt(err2))
            


            ibin+=1


        finalhisto.Draw('samee2')
        finalhisto.SetMarkerStyle(0)
        finalhisto.SetLineColor(1)
        finalhisto.SetFillStyle(3001)
        finalhisto.SetFillColor(1)

        
        self.keep.append(finalhisto)
        # Draw data
        data_view = self.get_view('data')
        if preprocess:
            data_view = preprocess( data_view )
        data_view = self.get_wild_dir(
            self.rebin_view(data_view, rebin),
            folder
            )
        data = data_view.Get(variable)
        data.Draw('same')
        print 'data', data.Integral()
        self.keep.append(data)
        ## Make sure we can see everything
        if data.GetMaximum() > mc_stack.GetMaximum():
            mc_stack.SetMaximum(1.2*data.GetMaximum()) 
            #mc_stack.SetMinimum(0.00001*data.GetMaximum())
            
            # # Ad legend
        #allentries = [data]
        #allentries.append( x for x in mc_stack.GetHists() )
        self.add_legend([data, mc_stack], leftside, entries=len(mc_stack.GetHists())+1)
        if show_ratio:
           self.add_ratio_diff(data, mc_stack, finalhisto, xrange, ratio_range)
                       
##-----

    def plot_without_uncert (self, folder, variable, rebin=1, xaxis='',
                        leftside=True, xrange=None, preprocess=None,
                              show_ratio=False, ratio_range=0.2, sort=True, obj=['e1', 'e2']):
        
        

        mc_stack_view = self.make_stack(rebin, preprocess, folder, sort)

        mc_stack = mc_stack_view.Get(variable)
        mc_stack.Draw()
        
        self.canvas.SetLogy(True)
        mc_stack.GetHistogram().GetXaxis().SetTitle(xaxis)
        if xrange:
            mc_stack.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            mc_stack.Draw()
        self.keep.append(mc_stack)
        

        finalhisto= mc_stack.GetStack().Last().Clone()
        finalhisto.Sumw2()

        histlist = mc_stack.GetHists();
        bkg_stack = mc_stack_view.Get(variable)
        bkg_stack.GetStack().RemoveLast()## mettere il check se c'e` il fake altirmenti histo=mc_stack.GetStack().Last().Clone()
        histo=bkg_stack.GetStack().Last().Clone()
        histo.Sumw2()
        
        fake_p1s_histo=None
      
        if not folder.startswith('tLoose'):
            isFakesIn= False
            if 'Fakes' in self.mc_samples:
                self.mc_samples.remove('Fakes')  
                self.mc_samples.remove('finalDYLL')
                isFakesIn=True
            
            
        ibin =1
            
        if isFakesIn:
            self.mc_samples.append('Fakes')
            self.mc_samples.append('finalDYLL')
          


        finalhisto.Draw('samee2')
        finalhisto.SetMarkerStyle(0)
        finalhisto.SetLineColor(1)
        finalhisto.SetFillStyle(3001)
        finalhisto.SetFillColor(1)

        
        self.keep.append(finalhisto)
        # Draw data
        data_view = self.get_view('data')
        if preprocess:
            data_view = preprocess( data_view )
        data_view = self.get_wild_dir(
            self.rebin_view(data_view, rebin),
            folder
            )
        data = data_view.Get(variable)
        data.Draw('same')
        print 'data', data.Integral()
        self.keep.append(data)
        ## Make sure we can see everything
        if data.GetMaximum() > mc_stack.GetMaximum():
            mc_stack.SetMaximum(1.2*data.GetMaximum()) 
            
            # # Add legend
        self.add_legend([data, mc_stack], leftside, entries=len(mc_stack.GetHists())+1)
        if show_ratio:
            self.add_ratio_diff(data, mc_stack, finalhisto, xrange, ratio_range)
            


    
 
  
