'''

Base class which makes nice plots.

Author: Evan K. Friis, UW

'''

import fnmatch
import re
import os
import math
import rootpy.plotting.views as views
from pdb import set_trace
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
from FinalStateAnalysis.PlotTools.MedianView     import MedianView
from FinalStateAnalysis.PlotTools.SystematicsView     import SystematicsView
from FinalStateAnalysis.StatTools.quad     import quad
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

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )

def histo_diff_quad(mc_err, *systematics):
    nbins = mc_err.GetNbinsX()
    clone = mc_err.Clone()
    sys_up = [i for i, _ in systematics]
    sys_dw = [i for _, i in systematics]

    #bin loop
    for ibin in range(nbins+2): #from uflow to oflow
        content = clone.GetBinContent(ibin)
        error   = clone.GetBinError(ibin)

        shifts_up = [abs(i.GetBinContent(ibin) - content) for i in sys_up]
        shifts_dw = [abs(i.GetBinContent(ibin) - content) for i in sys_dw]
        max_shift = [max(i, j) for i, j in zip(shifts_up, shifts_dw)]
        
        new_err = quad(error, *max_shift)
        clone.SetBinError(ibin, new_err)

    return clone

class BasePlotter(Plotter):
    def __init__ (self, blind_region=None, forceLumi=-1, use_embedded=False): 
        cwd = os.getcwd()
        self.period = '8TeV'
        self.sqrts  = 8
        jobid = os.environ['jobid']
        self.use_embedded = use_embedded
        print "\nPlotting e tau for %s\n" % jobid

        files     = glob.glob('results/%s/LFVHETauAnalyzerMVA/*.root' % jobid)
        lumifiles = glob.glob('inputs/%s/*.lumicalc.sum' % jobid)

        outputdir = 'plots/%s/lfvet' % jobid
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)

        samples = [os.path.basename(i).split('.')[0] for i in files]
        
        self.blind_region=blind_region
        if self.blind_region:
            # Don't look at the SS all pass region
            blinder = lambda x: BlindView(x, "os/.*",blind_in_range(*self.blind_region))

        super(BasePlotter, self).__init__(files, lumifiles, outputdir, blinder, forceLumi=forceLumi)

        self.views['fakes'] = {'view' : self.make_fakes()}
        self.mc_samples = [
            'GluGluToHToTauTau_M-125*', 
            'VBF_HToTauTau_M-125*',
            'TTJets*',
            'T*_t*',
            '[WZ][WZ]Jets',
            #'Wplus*Jets_madgraph*', #superseded by fakes
            'Z*jets_M50_skimmedLL',
            'Z*jets_M50_skimmedTT',
        ]

        if use_embedded:
            self.mc_samples.pop()
            self.views['Ztt_embedded'] = {'view' : self.make_embedded('os/gg/ept30/h_collmass_pfmet')}
        
    def make_fakes(self):
        '''Sets up the fakes view'''
        data_view = self.get_view('data')
        central_fakes = views.SubdirectoryView(data_view, 'tLoose')
        up_fakes = views.SubdirectoryView(data_view, 'tLooseUp')
        style = data_styles['Fakes*']
        return views.TitleView(
            views.StyleView(
                MedianView(highv=up_fakes, centv=central_fakes),
                **remove_name_entry(style)
            ),
            style['name']
        )
        
    def make_embedded(self, normalization_path):
        '''Configures the embedded view'''
        embedded_view = self.get_view('ZetauEmbedded_Run2012*', 'unweighted_view')
        zjets_view = self.get_view('Z*jets_M50_skimmedTT')

        embedded_histo = embedded_view.Get(normalization_path)
        zjets_histo = zjets_view.Get(normalization_path)

        embed_int = embedded_histo.Integral()
        zjets_int = zjets_histo.Integral()

        scale_factor = zjets_int / embed_int
        scaled_view = views.ScaleView(embedded_view, scale_factor)
        return scaled_view

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

        if self.blind_region:
            for ibin in err:
                if ibin >= band.FindBin(self.blind_region[0]) and ibin <= band.FindBin(self.blind_region[1]):
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
        band.SetFillStyle('x')
        band.SetFillColor(1)

        band.Draw('samee2')
       
        self.keep.append(data_clone)
        self.keep.append(band)
        self.keep.append(ref_function)
        self.pad.cd()
        return data_clone

    def add_shape_systematics(self, histo, path, view, folder_systematics = [], name_systematics = []):
        '''Adds shape systematics
        add_shape_systematics(self, histo, path, view, folder_systematics = [], name_systematics = []) --> histo
        histo is the central value 
        path is the path if the central value histo
        view contains all the systematics. 
        folder_systematics is a list of tuples with the folders containing shifts (up, down)
        name_systematics is a list of tuples containing the postfix to obtain the shifts (up, down)
        '''
        systematics = []
        for sys_up, sys_dw in folder_systematics:
            h_up = view.Get(os.path.join(sys_up, path))
            h_dw = view.Get(os.path.join(sys_dw, path))
            systematics.append(
                (h_up, h_dw)
                )

        #check if we have to apply met uncertainties
        for sys_up, sys_dw in name_systematics:
            h_up = view.Get(path + sys_up)
            h_dw = view.Get(path + sys_dw)
            systematics.append(
                (h_up, h_dw)
            )
        
        #ADD systematics
        return histo_diff_quad(histo, *systematics)
        

    def plot_with_bkg_uncert (self, folder, variable, rebin=1, xaxis='',
                              leftside=True, xrange=None, preprocess=None,
                              show_ratio=False, ratio_range=0.2, sort=True, obj=['e1', 'e2']):
        #xsection uncertainties
        #names must match with what defined in self.mc_samples
        xsec_unc_mapper = {
            'TTJets*' : 0.026,
            'T*_t*' : 0.041,
            '[WZ][WZ]Jets' : 0.056, #diboson
            'Wplus*Jets_madgraph*' : 0.035, #WJets
            'Z*jets_M50_skimmedLL' : 0.032,
            'Z*jets_M50_skimmedTT' : 0.032,
        }

        path = os.path.join(folder,variable)

        #make MC views with xsec error
        mc_views_nosys = self.mc_views(rebin, preprocess)
        mc_views = []
        for view, name in zip(mc_views_nosys, self.mc_samples):
            new = SystematicsView(
                view,
                xsec_unc_mapper.get(name, 0.) #default to 0
            )
            mc_views.append(new)

        #make MC stack
        mc_stack_view = views.StackView(*mc_views, sorted=sort) 
        mc_stack = mc_stack_view.Get( path )

        #make histo clone for err computation
        mc_sum_view = views.SumView(*mc_views)
        mc_err = mc_sum_view.Get( path )
        
        #Add MC-only systematics
        folder_systematics = [
            ('p1s', 'm1s'), #PU correction
        ]

        met_systematics = [
            ('_jes', '_jes_minus'), 
            ('_mes', '_mes_minus'), 
            ('_ees', '_ees_minus'), 
            ('_tes', '_tes_minus'), 
            ('_ues', '_ues_minus'), 
        ]

        name_systematics = [] #which are not MET
        #add MET sys if necessary
        if 'collmass' in variable.lower() or \
           'met' in variable.lower():
            c = 1
            #name_systematics.extend(met_systematics)

        #add them
        mc_err = self.add_shape_systematics(
            mc_err, 
            path, 
            mc_sum_view, 
            folder_systematics,
            name_systematics)

        #check if we are using the embedded sample
        if self.use_embedded:
            embed_view = self.get_view('Ztt_embedded')
            if preprocess:
                embed_view = preprocess(embed_view)
            embed_view = RebinView( embed_view, rebin)

            embed = embed_view.Get(path)

            #add xsec error
            embed = SystematicsView.add_error( embed, xsec_unc_mapper['Z*jets_M50_skimmedTT'])
            
            #add them to backgrounds
            mc_stack.Add(embed)
            mc_err += embed
            mc_sum_view = views.SumView(mc_sum_view, embed_view)

        #Add MC and embed systematics
        folder_systematics = [
            ('trp1s', 'trm1s'), #trig scale factor
        ]

        #Add as many eid sys as requested
        for name in obj:
            folder_systematics.extend([
                ('%sidp1s'  % name, '%sidm1s'  % name), #eID scale factor
                ('%sisop1s' % name, '%sisop1s' % name), #e Iso scale factor
            ])

        mc_err = self.add_shape_systematics(
            mc_err, 
            path, 
            mc_sum_view, 
            folder_systematics)

        #add lumi uncertainty
        mc_err = SystematicsView.add_error( mc_err, 0.026 )
        
        #get fakes
        fakes_view = self.get_view('fakes')
        if preprocess:
            fakes_view = preprocess(fakes_view)
        fakes_view = RebinView(fakes_view, rebin)

        fakes = fakes_view.Get(path)

        #add them to backgrounds
        mc_stack.Add(fakes)
        mc_err += fakes
        #set_trace()

        #draw stack
        mc_stack.Draw()
        self.keep.append(mc_stack)
        
        #set cosmetics
        self.canvas.SetLogy(True)
        mc_stack.GetHistogram().GetXaxis().SetTitle(xaxis)
        if xrange:
            mc_stack.GetXaxis().SetRangeUser(xrange[0], xrange[1])
            mc_stack.Draw()
              
        #set cosmetics and draw error
        mc_err.SetMarkerStyle(0)
        mc_err.SetLineColor(1)
        mc_err.SetFillStyle('x')
        mc_err.SetFillColor(1)
        mc_err.Draw('pe2 same')
        self.keep.append(mc_err)

        #Get signal
        signals = [
            'ggHiggsToETau',
            'vbfHiggsToETau',
        ]
        for name in signals:
            sig_view = self.get_view(name)
            if preprocess:
                sig_view = preprocess(sig_view)
            sig_view = RebinView(sig_view, rebin)
            
            histogram = sig_view.Get(path)
            histogram.Draw('same')
            self.keep.append(histogram)

        # Draw data
        data_view = self.get_view('data')
        if preprocess:
            data_view = preprocess( data_view )
        data_view = self.rebin_view(data_view, rebin)
        data = data_view.Get(path)

        data.Draw('same')
        print 'data', data.Integral()
        self.keep.append(data)

        ## Make sure we can see everything
        if data.GetMaximum() > mc_stack.GetMaximum():
            mc_stack.SetMaximum(1.2*data.GetMaximum()) 

        self.add_legend([data, mc_stack], leftside, entries=len(mc_stack.GetHists())+1)
        if show_ratio:
           self.add_ratio_plot(data, mc_err, xrange, ratio_range, True) # add_ratio_diff(data, mc_stack, mc_err, xrange, ratio_range)
                       
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
            


    
 
  
