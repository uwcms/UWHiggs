'''

Make inclusive Z->mumu control plots

'''

import os
import glob
from FinalStateAnalysis.PlotTools.Plotter import Plotter
import rootpy.plotting.views as views
import math
from FinalStateAnalysis.PlotTools.HistToTGRaphErrors import HistStackToTGRaphErrors
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors
import ROOT

ROOT.gROOT.SetBatch(True)

def quad(*xs):
    return math.sqrt(sum(x * x for x in xs))

class MedianView(object):
    ''' Takes high and low, returns median assigning half the diff as error. '''
    def __init__(self, highv, lowv):
        self.highv = highv
        self.centv = views.ScaleView( views.SumView(lowv, self.highv) , 0.5)

    def Get(self, path):
        central_hist = self.centv.Get(path)
        high_hist    = self.highv.Get(path)

        ret_hist = central_hist.Clone()
        for bin in range(1, high_hist.GetNbinsX() + 1):
            error = quad(
                central_hist.GetBinError(bin),
                (high_hist.GetBinContent(bin) - central_hist.GetBinContent(bin))
            )
            ret_hist.SetBinError(bin, error)
        return ret_hist

class ControlZEEPlotter(Plotter):
    def __init__(self):
        self.jobid = os.environ['jobid']
        self.output_dir = os.path.join('results', self.jobid, 'plots', 'zee')
        samples = [
            'Zjets_M50',
            #'WZ*',
            #'ZZ*',
            #'WW*',
            'TT*',
            #'WplusJets*',
            "data_DoubleE*",
            ]
        files = []
        lumifiles = []
        self.sqrts = 7 if '7TeV' in self.jobid else 8
        for x in samples:
            files.extend(glob.glob('results/%s/ControlZEE/%s.root' % (self.jobid, x)))
            lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (self.jobid, x)))
        super(ControlZEEPlotter, self).__init__(files, lumifiles, self.output_dir)
        self.mc_samples = ['Zjets_M50']

    def get_flip_data(self, rebin=1, xaxis=''):
        data_view = self.get_view('data')
        data_view = self.rebin_view(data_view, rebin) if rebin != 1 else data_view
        #Get ss/p1p2 views
        ss_p1p2_view = views.SubdirectoryView( data_view, 'ss/p1p2')
        ss_p1p2_view = views.TitleView( views.StyleView(ss_p1p2_view, **data_styles['data*']), 'observed;%s' % xaxis)

        def make_fakes_view(sign, weight_type):
            # View of weighted obj1-fails data
            obj1_view = views.SubdirectoryView(data_view, '%s/f1p2/%s' % (sign, weight_type))
            # View of weighted obj2-fails data
            obj2_view = views.SubdirectoryView(data_view, '%s/p1f2/%s' % (sign, weight_type))
            # View of weighted obj1&2-fails data
            obj12_view = views.SubdirectoryView(data_view, '%s/f1f2/%s' % (sign, weight_type))
            # Give the individual object views nice colors
            subtract_obj12_view = views.ScaleView(obj12_view, -1)
            return obj1_view, obj2_view, subtract_obj12_view

        #Get fakes according to WJets or QCD
        ss_f1p2_qcd_view, ss_p1f2_qcd_view, ss_f1f2_qcd_view = make_fakes_view('ss','qcd_w')
        ss_f1p2_wje_view, ss_p1f2_wje_view, ss_f1f2_wje_view = make_fakes_view('ss','wjet_w')

        ss_fakes_1   = MedianView(ss_f1p2_qcd_view, ss_f1p2_wje_view)
        ss_fakes_2   = MedianView(ss_p1f2_qcd_view, ss_p1f2_wje_view)
        ss_fakes_12  = MedianView(ss_f1f2_qcd_view, ss_f1f2_wje_view)
        ss_fakes_est = views.SumView(ss_fakes_1, ss_fakes_2, ss_fakes_12)
        ss_fakes_est = views.TitleView( views.StyleView(ss_fakes_est, **data_styles['Zjets*']), 'Fakes;%s' % xaxis)

        os_f1p2_qcd_views, os_p1f2_qcd_views, os_f1f2_qcd_views = zip(make_fakes_view('os', 'qcd_w/charge_weightSysUp'), make_fakes_view('os', 'qcd_w/charge_weightSysDwn'))
        os_f1p2_wje_views, os_p1f2_wje_views, os_f1f2_wje_views = zip(make_fakes_view('os', 'wjet_w/charge_weightSysUp'),make_fakes_view('os', 'wjet_w/charge_weightSysDwn'))

        os_f1p2_qcd_view = MedianView( *os_f1p2_qcd_views )
        os_p1f2_qcd_view = MedianView( *os_p1f2_qcd_views )
        os_f1f2_qcd_view = MedianView( *os_f1f2_qcd_views )

        os_f1p2_wje_view = MedianView( *os_f1p2_wje_views )
        os_p1f2_wje_view = MedianView( *os_p1f2_wje_views )
        os_f1f2_wje_view = MedianView( *os_f1f2_wje_views )

        os_fakes_1   = MedianView(os_f1p2_qcd_view, os_f1p2_wje_view)
        os_fakes_2   = MedianView(os_p1f2_qcd_view, os_p1f2_wje_view)
        os_fakes_12  = MedianView(os_f1f2_qcd_view, os_f1f2_wje_view)
        os_fakes_est = views.SumView(os_fakes_1, os_fakes_2, os_fakes_12)
        neg_os_fakes = views.ScaleView(os_fakes_est, -1)

        os_flip_est_up  = views.SubdirectoryView( data_view, 'os/p1p2/charge_weightSysUp')
        os_flip_est_dwn = views.SubdirectoryView( data_view, 'os/p1p2/charge_weightSysDwn')
        os_flip_est     = MedianView(os_flip_est_up, os_flip_est_dwn)
        os_flip_est_nofake = os_flip_est#views.SumView(os_flip_est, neg_os_fakes)
        os_flip_est_nofake = views.TitleView( views.StyleView(os_flip_est_nofake, **data_styles['WZ*']), 'charge-fakes;%s' % xaxis)
        return ss_p1p2_view, ss_fakes_est, os_flip_est_nofake
            
    def make_charge_flip_control_plot(self, variable, xaxis='', rebin=1):
        ss_p1p2_view, ss_fakes_est, os_flip_est_nofake = self.get_flip_data(rebin,xaxis)
        events_estimate = views.StackView( ss_fakes_est,os_flip_est_nofake)
        
        obs_hist       = ss_p1p2_view.Get(variable)
        estimate_hist  = events_estimate.Get(variable)
        estimate_error = HistStackToTGRaphErrors( estimate_hist )
        estimate_error.SetFillStyle(3013)
        estimate_error.SetFillColor(ROOT.EColor.kBlack)
        estimate_error.SetTitle('Error on estimate')

        print variable, estimate_hist.GetMaximum(), max(list(obs_hist)), max( [ estimate_hist.GetMaximum(), max(list(obs_hist)) ] )
        hmax = max( [ estimate_hist.GetMaximum(), max(list(obs_hist)) ] )
        obs_hist.GetYaxis().SetRangeUser(0,hmax*1.3)

        obs_hist.Draw()
        estimate_hist.Draw('same')
        self.canvas.Update()
        estimate_error.Draw('2 same')
        obs_hist.Draw('same')
        self.keep.extend([
            estimate_hist,
            estimate_error,
            obs_hist
            ])

        legend = self.add_legend([obs_hist], leftside=False, entries=4)
        legend.AddEntry(estimate_hist,'f')
        #legend.AddEntry(estimate_error,'f')        
        legend.Draw()
        self.add_cms_blurb(self.sqrts)



plotter = ControlZEEPlotter()

plotter.make_charge_flip_control_plot('TrkMass','Tracker Inv Mass (GeV)',2)
plotter.save('EE_Charge_Flip_xcheck_trk_invMass')

plotter.make_charge_flip_control_plot('TrkMass_NOSCALE','Tracker Inv Mass (GeV)',2)
plotter.save('EE_Charge_Flip_xcheck_trk_invMass_NoScale')

plotter.make_charge_flip_control_plot('SCMass','SuperCluster Inv Mass (GeV)',2)
plotter.save('EE_Charge_Flip_xcheck_SC_invMass')

plotter.make_charge_flip_control_plot('ePt','electron p_{T}',2)
plotter.save('EE_Charge_Flip_xcheck_ePt')

plotter.make_charge_flip_control_plot('eAbsEta','electron |#eta|',4)
plotter.save('EE_Charge_Flip_xcheck_eAbsEta')

plotter.make_charge_flip_control_plot('SCDPhi','Super Cluster #Delta#phi', 6)
plotter.save('EE_Charge_Flip_xcheck_eSCDPhi')

plotter.make_charge_flip_control_plot('SCEnergy','Super cluster energy (GeV)',5)
plotter.save('EE_Charge_Flip_xcheck_eSCEnergy')


## plotter.plot_mc_vs_data('zmm', 'm1m2Mass', rebin=2, xaxis='m_{#mu#mu} (GeV)')
## plotter.add_cms_blurb(sqrts)
## plotter.save('mass')
## 
## plotter.plot_mc_vs_data('zmm', 'm1m2Mass', rebin=6, xaxis='m_{#mu#mu} (GeV)')
## plotter.add_cms_blurb(sqrts)
## plotter.save('mass_rebin')
## 
## plotter.plot_mc_vs_data('zmm', 'm1Pt')
## plotter.save('m1Pt')
## plotter.plot_mc_vs_data('zmm', 'm1Pt', 5)
## plotter.save('m1Pt_rebin')
## plotter.plot_mc_vs_data('zmm', 'm2Pt')
## plotter.save('m2Pt')
## plotter.plot_mc_vs_data('zmm', 'm2Pt', 5)
## plotter.save('m2Pt_rebin')
## 
## plotter.plot_mc_vs_data('zmm', 'm1AbsEta')
## plotter.save('m1AbsEta')
## plotter.plot_mc_vs_data('zmm', 'm2AbsEta')
## plotter.save('m2AbsEta')
## 
## plotter.plot_mc_vs_data('zmm', 'm1AbsEta', 5)
## plotter.save('m1AbsEta_rebin')
## plotter.plot_mc_vs_data('zmm', 'm2AbsEta', 5)
## plotter.save('m2AbsEta_rebin')
## 
## plotter.plot_mc_vs_data('zmm', 'nvtx')
## plotter.save('nvtx')
