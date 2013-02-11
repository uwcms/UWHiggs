import rootpy.plotting.views as views
import rootpy.plotting as plotting
#from rootpy.tree import TreeChain
import rootpy.io
from FinalStateAnalysis.PlotTools.Plotter import Plotter
from FinalStateAnalysis.PlotTools.BlindView import BlindView
from FinalStateAnalysis.PlotTools.PoissonView import PoissonView
from FinalStateAnalysis.PlotTools.HistToTGRaphErrors import HistToTGRaphErrors, HistStackToTGRaphErrors
from FinalStateAnalysis.PlotTools.InflateErrorView import InflateErrorView
from FinalStateAnalysis.MetaData.data_styles import data_styles
from FinalStateAnalysis.PlotTools.THBin import zipBins
#from RecoLuminosity.LumiDB import argparse''' 
import sys
import os
import glob
import pprint
import ROOT

import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)

def quad(*xs):
    return math.sqrt(sum(x*x for x in xs))


class WHChargeFlipControlPlotsPlotter(Plotter):
    def __init__(self, channel):
        self.samples = [ 'Zjets*', 'WZ*', 'ZZ*']#'WH*',
        self.samples += ['data*']
        self.jobid = os.environ['jobid']
        self.channel = channel
        self.period = '7TeV' if '7TeV' in self.jobid else '8TeV'
        self.sqrts = 7 if '7TeV' in self.jobid else 8
        files = []
        lumifiles = []
        for x in self.samples:
            files += glob.glob('results/%s/FakeRates%s/%s.root' % (self.jobid, channel, x))
            lumifiles += glob.glob('inputs/%s/%s.lumicalc.sum' % (self.jobid, x))
        self.outputdir = os.path.join('results', self.jobid, 'fakerate_fits')
        pprint.pprint(files)
        super(WHChargeFlipControlPlotsPlotter, self).__init__(files, lumifiles, self.outputdir, None)

    def plot_mc_vs_data(self, folder, variable, rebin=1, xaxis='', leftside=True, xrange=None):
        super(ZHPlotterBase, self).plot_mc_vs_data(folder, variable, rebin, xaxis, leftside, xrange)
        self.add_cms_blurb(self.sqrts)

    def make_views(self, rebin):
        ''' Make signal views with FR background estimation '''

        wz_view = InflateErrorView(
            views.SubdirectoryView(
                self.rebin_view(self.get_view('WZ*'), rebin),
                'charge'
                ),
            0.16
            )
                
        zz_view = InflateErrorView(
            views.SubdirectoryView(
                self.rebin_view(self.get_view('ZZ*'), rebin),
                'charge'
                ),
            0.16
            )

        zj_view = views.SubdirectoryView(
            self.rebin_view(self.get_view('Zjets*'), rebin),
            'charge'
            )

        diboson_view = views.StyleView( views.SumView(wz_view,zz_view), fillcolor=(ROOT.EColor.kRed-7), drawstyle='hist' )
        all_data_view =views.SubdirectoryView( self.rebin_view(self.get_view('data'), rebin), 'charge')
            
        output = {
            'wz'      : wz_view,
            'zz'      : zz_view,
            'zjets'   : zj_view,
            'data'    : all_data_view,
            'diboson' : diboson_view,
        }

        return output

    def apply_error_from_sysHist(self, hlow, hcent, hhigh):
        if hlow.GetNbinsX() == hcent.GetNbinsX() == hhigh.GetNbinsX():
            for binl, binm, binh in zipBins(hlow, hcent, hhigh):
                binm.error = max( abs(binm.content-binl.content), abs(binm.content-binh.content) )
        else:
            raise BinError("The tree input histograms don't have the same binning")

    @staticmethod
    def hists_to_graph_asymm_err(hlow, hcent, hhigh):
        'convert three histograms into one TGraphAsymmError'
        if hlow.GetNbinsX() == hcent.GetNbinsX() == hhigh.GetNbinsX():
            binsx = []
            binsy = []
            errx  = []
            erryl = []
            erryh = []
            for binl, binm, binh in zipBins(hlow, hcent, hhigh):
                binsx.append( binm.center )
                binsy.append( binm.content)
                errx.append(  binm.width/2. )
                erryl.append( binm.content - binl.content)
                erryh.append( binh.content - binm.content)
            
            return ROOT.TGraphAsymmErrors(len(binsx),
                                          array('f',binsx),
                                          array('f',binsy),
                                          array('f',errx),
                                          array('f',errx),
                                          array('f',erryl),
                                          array('f',erryh)
                )
        else:
            raise BinError("The tree input histograms don't have the same binning")

    def make_xcheck_plot(self, var, sample, **kwargs):
        self.canvas.cd()
        rebin = kwargs['rebin'] if 'rebin' in kwargs else 1
        sample_views = self.make_views(rebin)
        if sample not in sample_views:
            print "Sample not found! Available: \nSkipping..." % sample_views.keys().__repr__()
            return

        testH       = sample_views[sample].Get('SS_'+var)
        testH.drawstyle = 'ep'

        weightedH   = sample_views[sample].Get('OS_weight_'+var)
        weightedHup = sample_views[sample].Get('OS_weightSysUp_'+var)
        weightedHdw = sample_views[sample].Get('OS_weightSysDwn_'+var)
        print "%s   %s: %i Events, %i Entries. Weighted from sideband: %.1f Events." % (sample, var+' observed', testH.Integral(), testH.GetEntries(),  weightedH.Integral())
        self.apply_error_from_sysHist(weightedHdw, weightedH, weightedHup)
        diboson = None

        graph = HistToTGRaphErrors(weightedH)
        graph.SetFillColor(ROOT.kAzure -4)
        self.keep.append(graph)
        if 'xrange' in kwargs:
            graph.GetXaxis().SetRangeUser(kwargs['xrange'][0],kwargs['xrange'][1])
        if 'yrange' in kwargs:
            graph.GetYaxis().SetRangeUser(kwargs['yrange'][0],kwargs['yrange'][1])

        graph.Draw("2A")
        graph.GetYaxis().SetRangeUser(0, max(testH.GetMaximum(),weightedH.GetMaximum())*1.2)
        testH.Draw(' same')



##########################################
##    EE Control Plots
##########################################
plotter = WHChargeFlipControlPlotsPlotter('EE')

plotter.make_xcheck_plot('TrkMass','data',rebin=2)
plotter.save('EE_Charge_Flip_xcheck_trk_invMass')

plotter.make_xcheck_plot('SCMass','data',rebin=2)
plotter.save('EE_Charge_Flip_xcheck_SC_invMass')

plotter.make_xcheck_plot('ePt','data',rebin=2)
plotter.save('EE_Charge_Flip_xcheck_ePt')

plotter.make_xcheck_plot('eAbsEta','data',rebin=2)
plotter.save('EE_Charge_Flip_xcheck_eAbsEta')

plotter.make_xcheck_plot('SCDPhi','data', rebin=6)
plotter.save('EE_Charge_Flip_xcheck_eSCDPhi')

plotter.make_xcheck_plot('SCEnergy','data',rebin=5, xrange=[0,600])
plotter.save('EE_Charge_Flip_xcheck_eSCEnergy')


##########################################
##    Zjets Closure test
##########################################

plotter.make_xcheck_plot('TrkMass','zjets',rebin=2)
plotter.save('EE_Charge_Flip_closure_test_zjets_trk_invM')

plotter.make_xcheck_plot('SCMass','data',rebin=2)
plotter.save('EE_Charge_Flip_closure_test_zjets_SC_invMass')

plotter.make_xcheck_plot('ePt','zjets',rebin=2)
plotter.save('EE_Charge_Flip_closure_test_zjets_ePt')

plotter.make_xcheck_plot('eAbsEta','zjets',rebin=2)
plotter.save('EE_Charge_Flip_closure_test_zjets_eAbsEta')

plotter.make_xcheck_plot('SCDPhi','zjets', rebin=6)
plotter.save('EE_Charge_Flip_closure_test_zjets_eSCDPhi')

plotter.make_xcheck_plot('SCEnergy','zjets', rebin=5, xrange=[0,600])
plotter.save('EE_Charge_Flip_closure_test_zjets_SCEnergy')
