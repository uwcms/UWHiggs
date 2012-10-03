'''

Base class to do WH plotting.

Author: Evan K. Friis, UW

Takes as input a set of ROOT [files] with analysis histgrams, and the corresponding
lumicalc.sum [lumifiles] that hve the effective lumi for each sample.

If [blind] is true, data in the p1p2p3 region will not be plotted.

'''

import rootpy.plotting.views as views
from FinalStateAnalysis.PlotTools.Plotter import Plotter
from FinalStateAnalysis.PlotTools.BlindView import BlindView
from FinalStateAnalysis.PlotTools.PoissonView import PoissonView
from FinalStateAnalysis.MetaData.data_styles import data_styles
#from RecoLuminosity.LumiDB import argparse''' 
import sys
import os
import glob
from THBin import zipBins
import ROOT

import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)

def quad(*xs):
    return math.sqrt(sum(x*x for x in xs))

class BackgroundErrorView(object):
    ''' Compute the total background error in each bin. '''
    def __init__(self, fakes, wz, zz, wz_error=0.1, zz_error=0.04, fake_error=0.3):
        self.fakes = fakes
        self.wz = wz
        self.zz = zz
        self.fake_error = fake_error
        self.wz_error = wz_error
        self.zz_error = zz_error

    def Get(self, path):
        fakes = self.fakes.Get(path)
        wz = self.wz.Get(path)
        zz = self.zz.Get(path)

        bkg_error = wz.Clone()
        bkg_error.SetTitle("Bkg. Unc.")
        bkg_error.Reset()
        for bin_err, bin_wz, bin_zz, bin_fake in zipBins(bkg_error, wz, zz, fakes): #range(1, bkg_error.GetNbinsX() + 1):
            bin_err.error = quad(
                bin_fake.error,
                bin_fake.content*self.fake_error,
                bin_wz.content*self.wz_error,
                bin_zz.content*self.zz_error
            )
            bin_err.content = (
                bin_fake.content(bin) +
                bin_wz.content(bin) +
                bin_zz.content(bin)
            )
        bkg_error.SetMarkerSize(0)
        bkg_error.SetFillColor(1)
        bkg_error.SetFillStyle(3013)
        bkg_error.legendstyle = 'f'
        return bkg_error

def guess_var_content(channel):
    chSet  = set(channel)
    chSort = sorted(chSet)
    ret    = []
    for char in chSort:
        nc = chSort.count(char)
        if nc == 1:
            ret.append( char.lower() )
        else:
            for i in range(1,nc+1):
                ret.append( '%s%s' % (char.lower(), i) )
    return ret

def get_mass_histos(channel):
    varcont = guess_var_content(channel)
    for pos, obj1 in enumerate( varcont ):
        for obj2 in varcont[pos+1:]:
            yield "%s_%s_Mass" % (obj1, obj2)

products_map = { #a finction might be cooler but is totally a bummer
                 'MMMT' : (('m1', 'm2'),('m3','t')),
                 'MMET' : (('m1', 'm2'),('e','t')),
                 'MMEM' : (('m1', 'm2'),('e','m3')),
                 'MMTT' : (('m1', 'm2'),('t1','t2')),
                 'EEMT' : (('e1', 'e2'),('m','t')),
                 'EEET' : (('e1', 'e2'),('e3','t')),
                 'EEEM' : (('e1', 'e2'),('e3','m')),
                 'EETT' : (('e1', 'e2'),('t1','t2')),
                 }

naming_map = {
    'm' : '#mu',
    'e' : 'e',
    't' : '#tau',
    }

var_map = {
    'Mass'  : 'M_{%s%s}',
    'Pt'    : 'p_{T%s}',
    'Abseta': '|#eta_{%s}|',
    }


class ZHPlotterBase(Plotter):
    def __init__(self, channel, blind=False):
        self.samples = [ 'Zjets_M50', 'WplusJets_madgraph', 'WZJetsTo3LNu*', 'ZZ*', 'WW*', 'VH*', 'WH*', 'TTplusJets_madgraph']
        self.samples += ['data_DoubleMu*'] if channel[:2] == 'MM' else ['data_DoubleElectron*']
        self.jobid = os.environ['jobid']
        self.channel = channel
        self.period = '7TeV' if '7TeV' in jobid else '8TeV'
        self.sqrts = 7 if '7TeV' in jobid else 8
        files = []
        self.blind = blind
        lumifiles = []
        for x in self.samples:
            files += glob.glob('results/%s/ZHAnalyze%s/%s.root' % (jobid, channel, x))
            lumifiles += glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x))
        self.outputdir = 'results/%s/plots/%s' % (jobid, channel.lower() )
        blinder = None
        if blind:
            # Don't look at the SS all pass region
            blinder = lambda x: BlindView(x, "ss/p1p2p3/.*")
        super(ZHPlotterBase, self).__init__(files, lumifiles, self.outputdir, blinder)
        self.general_histos = ["nTruePU", "weight", "rho", "nvtx"]
        self.general_histos += ['doubleMuPrescale'] if channel[:2] == 'MM' else []
        self.kin_histos     = ["%sPt", "%sJetPt", "%sAbsEta"]
        self.mass_histos    = [h for h in get_mass_histos(channel)]

    def plot_mc_vs_data(self, folder, variable, rebin=1, xaxis='', leftside=True, xrange=None):
        super(ZHPlotterBase, self).plot_mc_vs_data(folder, variable, rebin, xaxis, leftside, xrange)
        self.add_cms_blurb(self.sqrts)


    def make_signal_views(self, rebin, unblinded=True):
        ''' Make signal views with FR background estimation '''

        wz_view = views.SubdirectoryView(
            self.rebin_view(self.get_view('WZJetsTo3LNu*'), rebin),
            'os/All_Passed/'
        )
        zz_view = views.SubdirectoryView(
            self.rebin_view(self.get_view('ZZJetsTo4L*'), rebin),
            'os/All_Passed/'
        )
        all_data_view =self.rebin_view(self.get_view('data'), rebin)
        if unblinded and self.blind:
            all_data_view =self.rebin_view(
                self.get_view('data', 'unblinded_view'), rebin)
            
        data_view = views.SubdirectoryView(all_data_view, 'os/All_Passed/')
        #Categories (to match Abdollah's naming convention)
        probes = [p+'IsoFailed' for p in products_map[self.channel][1]]
        cat0   = 'os/'+ '_'.join(probes) + '/all_weights_applied/'
        cat1   = 'os/'+ probes[0] + '/obj1_weight/'
        cat2   = 'os/'+ probes[1] + '/obj2_weight/'

        # View of weighted obj1-fails data
        cat1_view = views.SubdirectoryView(all_data_view, cat1)
        # View of weighted obj2-fails data
        cat2_view = views.SubdirectoryView(all_data_view, cat2)
        # View of weighted obj1&2-fails data
        cat0_view = views.SubdirectoryView(all_data_view, cat0)

        subtract_cat0_view = views.ScaleView(cat0_view, -1)
        # Corrected fake view
        fakes_view = views.SumView(cat1_view, cat2_view, subtract_cat0_view)
        fakes_view = views.TitleView(
            views.StyleView(fakes_view, **data_styles['Zjets*']), 'Non-prompt')

        ## charge_fakes = views.TitleView(
        ##     views.StyleView(
        ##         views.SubdirectoryView(all_data_view, 'os/p1p2p3/c1'),
        ##         **data_styles['TT*']), 'Charge mis-id')

        output = {
            'wz' : wz_view,
            'zz' : zz_view,
            'data' : data_view,
            'cat1' : cat1_view,
            'cat2' : cat2_view,
            'fakes' : fakes_view,
            ## 'charge_fakes' : charge_fakes,
        }

        # Add signal
        for mass in [110, 120, 130, 140]:
            vh_view = views.SubdirectoryView(
               self.rebin_view(self.get_view('VH_*%i' % mass), rebin),
                'os/All_Passed/'
            )
            output['vh%i' % mass] = vh_view
            ww_view = views.SubdirectoryView(
               self.rebin_view(self.get_view('WH_%i*' % mass), rebin),
                'os/All_Passed/'
            )
            output['vh%i_hww' % mass] = ww_view
            output['signal%i' % mass] = views.SumView(ww_view, vh_view)

        return output

    def write_shapes(self, variable, rebin, outdir, unblinded=False):
        ''' Write final shape histos for [variable] into a TDirectory [outputdir] '''
        sig_view = self.make_signal_views(rebin, unblinded)
        outdir.cd()
        wz = sig_view['wz'].Get(variable)
        zz = sig_view['zz'].Get(variable)
        obs = sig_view['data'].Get(variable)
        fakes = sig_view['fakes'].Get(variable)

        wz.SetName('wz')
        zz.SetName('zz')
        obs.SetName('data_obs')
        fakes.SetName('fakes')

        for mass in [110, 120, 130, 140]:
            vh = sig_view['vh%i' % mass].Get(variable)
            vh.SetName('VH%i' % mass)
            vh.Write()
            ww = sig_view['vh%i_hww' % mass].Get(variable)
            ww.SetName('VH_hww%i' % mass)
            ww.Write()

        wz.Write()
        zz.Write()
        obs.Write()
        fakes.Write()

    def write_cut_and_count(self, variable, outdir, unblinded=False):
        ''' Version of write_shapes(...) with only one bin.

        Equivalent to a cut & count analysis.
        '''
        sig_view = self.make_signal_views(1, unblinded)
        nbins = sig_view['wz'].Get(variable).GetNbinsX()
        return self.write_shapes(variable, nbins, outdir, unblinded)


    def plot_final(self, variable, rebin=1, xaxis='', maxy=10, show_error=False, magnifyHiggs=5 ):
        ''' Plot the final output - with bkg. estimation '''
        
        sig_view = self.make_signal_views(rebin)
        vh_nx = views.TitleView(
            views.StyleView(
                views.ScaleView(sig_view['signal120'], magnifyHiggs),
                **data_styles['VH*']
            ),
            "(%s#times) m_{H} = 120" % magnifyHiggs
        )

        stack = views.StackView(
            sig_view['wz'],
            sig_view['zz'],
            sig_view['fakes'],
            vh_10x,
        )
        histo = stack.Get(variable)
        histo.Draw()
        histo.GetHistogram().GetXaxis().SetTitle(xaxis)
        histo.SetMaximum(maxy)
        self.keep.append(histo)

        # Add legend
        legend = self.add_legend(histo, leftside=False, entries=4)

        if show_error:
            bkg_error_view = BackgroundErrorView(
                sig_view['fakes'],
                sig_view['wz'],
                sig_view['zz'],
            )
            bkg_error = bkg_error_view.Get(variable)
            self.keep.append(bkg_error)
            bkg_error.Draw('pe2,same')
            legend.AddEntry(bkg_error)

        # Use poisson error bars on the data
        sig_view['data'] = PoissonView(sig_view['data'], x_err=False)

        data = sig_view['data'].Get(variable)
        data.Draw('pe,same')
        self.keep.append(data)

        #legend.AddEntry(data)
        legend.Draw()

        
if __name__ <> "__main__":
    sys.exit(0)
jobid = os.environ['jobid']

channels = filter( lambda x: len(x) == 4 and x.upper() == x, sys.argv) #['MMMT']#[ Z+H for Z in ['EE','MM'] for H in ['MT', 'ET', 'TT', 'EM']] #
print channels
for channel in channels:
    print "Plotting %s for %s" % (channel, jobid)
    Zprod, Hprod = products_map[channel]
    texZprod = tuple( naming_map[p[0]] for p in Zprod )
    texHprod = tuple( naming_map[p[0]] for p in Hprod )
    plotter = ZHPlotterBase(channel)
    
    ###########################################################################
    ##  Z control plots #####################################################
    ###########################################################################
    plotter.plot_mc_vs_data('os/All_Passed', '%s_%s_Pt' % Zprod, rebin=10, xaxis='p_{T%s%s} (GeV)' % texZprod, leftside=False)
    plotter.save('mcdata-os-all_passed_ZPt')

    plotter.plot_mc_vs_data('os/All_Passed', '%s_%s_Mass' % Zprod, rebin=10, xaxis='M_{%s%s} (GeV)' % texZprod, leftside=False)
    plotter.save('mcdata-os-all_passed_ZMass')

    ###########################################################################
    ##  H control plots #####################################################
    ###########################################################################
    plotter.plot_mc_vs_data('os/All_Passed', '%s_%s_Pt' % Hprod, rebin=10, xaxis='p_{T%s%s} (GeV)' % texHprod, leftside=False)
    plotter.save('mcdata-os-all_passed_HPt')

    plotter.plot_mc_vs_data('os/All_Passed', '%s_%s_Mass' % Hprod, rebin=10, xaxis='M_{%s%s} (GeV)' % texHprod, leftside=False)
    plotter.save('mcdata-os-all_passed_HMass')

    ###########################################################################
    ##  Making shape file     #################################################
    ###########################################################################

    shape_file = ROOT.TFile( os.path.join(plotter.outputdir, '%s_shapes_%s.root' % (channel.lower(), plotter.period)), 'RECREATE')
    shape_dir  = shape_file.mkdir( channel.lower() )
    plotter.write_shapes('%s_%s_Mass' % Hprod, 20, shape_dir, unblinded=True)
    #plotter.write_cut_and_count('subMass', shape_dir, unblinded=True)
    shape_file.Close()

    ## plotter.plot_mc_vs_data('os/p1p2f3', 'nTruePU', rebin=1, xaxis='True PU')
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('mcdata-os-p1p2f3-nTruePU')

    ## plotter.plot('Zjets_M50', 'os/p1p2f3/nTruePU', 'nTruePU', rebin=1, xaxis='True PU')
    ## plotter.save('zjets-os-p1p2f3-nTruePU')


    ## plotter.plot_mc_vs_data('os/p1p2f3', 'bCSVVeto', rebin=1, xaxis='bveto')
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('mcdata-os-p1p2f3-bveto')

    ## plotter.plot_mc_vs_data('os/p1p2f3/w3', 'emMass', 10)
    ## plotter.save('mcdata-os-p1p2f3-w3-emMass')

    ## plotter.plot_mc_vs_data('os/p1f2p3', 'emMass', 10)
    ## plotter.save('mcdata-os-p1f2p3-emMass')

    ## plotter.plot_mc_vs_data('os/f1p2p3', 'emMass', 10)
    ## plotter.save('mcdata-os-p1f2p3-emMass')

    ## # Check PU variables
    ## #plotter.plot_mc_vs_data('os/p1p2f3', 'rho')
    ## #plotter.save('mcdata-os-p1p2f3-rho')

    ## #plotter.plot_mc_vs_data('os/p1p2f3', 'nvtx')
    ## #plotter.save('mcdata-os-p1p2f3-nvtx')

    ## # Lower stat but closer to signal region
    ## #plotter.plot_mc_vs_data('os/p1p2p3', 'rho')
    ## #plotter.save('mcdata-os-p1p2p3-rho')

    ## #plotter.plot_mc_vs_data('os/p1p2p3', 'nvtx')
    ## #plotter.save('mcdata-os-p1p2p3-nvtx')

    ## # Make Z->mumu + tau jet control
    ## def make_styler(color, format=None):
    ##     def unsuck(x):
    ##         x.SetFillStyle(0)
    ##         x.SetLineColor(color)
    ##         x.SetLineWidth(2)
    ##         x.SetMaximum(1.5*x.GetMaximum())
    ##         if format:
    ##             x.format = format
    ##     return unsuck

    ## weighted = plotter.plot('data', 'os/p1p2f3/w3/emMass',  'hist', rebin=20, styler=make_styler(2, 'hist'), xaxis='m_{e#mu} (GeV)')
    ## unweighted = plotter.plot('data', 'os/p1p2p3/emMass', 'same', rebin=20, styler=make_styler(1), xaxis='m_{e#mu} (GeV)')
    ## weighted.SetTitle('e^{+}#mu^{-} + fake #tau_{h} est.')
    ## weighted.legendstyle = 'l'
    ## unweighted.SetTitle('e^{+}#mu^{-} + fake #tau_{h} obs.')
    ## unweighted.legendstyle = 'pe'
    ## plotter.add_legend([weighted, unweighted])
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('ztt-os-fr-control')

    ## #plotter.plot('data', 'os/p1p2p3/prescale', styler=make_styler(1))
    ## #plotter.save('ztt-os-prescale-check')

    ## #plotter.plot('Zjets_M50', 'os/p1p2f3/weight')
    ## #plotter.save('ztt-mc-event-weights')
    ## ## Check MC weights
    ## #plotter.plot('Zjets_M50', 'os/p1p2f3/weight_nopu')
    ## #plotter.save('ztt-mc-event-weight_nopu')


    ## ###########################################################################
    ## ##  FR sideband MC-vs-Data ################################################
    ## ###########################################################################

    ## plotter.plot_mc_vs_data('ss/p1f2p3', 'mPt', 5, '#mu_{1} p_{T}', leftside=False)
    ## plotter.save('mcdata-ss-p1f2p3-mPt')

    ## plotter.plot_mc_vs_data('ss/p1f2p3', 'subMass', 20, 'Subleading mass (GeV)', leftside=False)
    ## plotter.save('mcdata-ss-p1f2p3-subMass')

    ## plotter.plot_mc_vs_data('ss/p1p2f3', 'subMass', 20, 'Subleading mass (GeV)', leftside=False)
    ## plotter.save('mcdata-ss-p1p2f3-subMass')

    ## plotter.plot_mc_vs_data('ss/f1p2p3', 'subMass', 20, 'Subleading mass (GeV)', leftside=False)
    ## plotter.save('mcdata-ss-f1p2p3-subMass')

    ## plotter.plot_mc_vs_data('ss/p1f2p3/w2', 'mPt', 5, '#mu_{1} p_{T}', leftside=False)
    ## plotter.save('mcdata-ss-p1f2p3-w2-mPt')

    ## plotter.plot_mc_vs_data('ss/p1f2p3', 'ePt', 5, 'Electron p_{T}', leftside=False)
    ## plotter.save('mcdata-ss-p1f2p3-ePt')

    ## plotter.plot_mc_vs_data('ss/p1f2p3/w2', 'ePt', 5, 'Electron p_{T}', leftside=False)
    ## plotter.save('mcdata-ss-p1f2p3-w2-ePt')

    ## plotter.plot_mc_vs_data('ss/f1p2p3', 'ePt', 5, 'Electron p_{T}', leftside=False)
    ## plotter.save('mcdata-ss-f1p2p3-ePt')

    ## plotter.plot_mc_vs_data('ss/f1p2p3/w1', 'ePt', 5, 'Electron p_{T}', leftside=False)
    ## plotter.save('mcdata-ss-f1p2p3-w2-ePt')


