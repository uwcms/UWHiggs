'''

Base class to do WH plotting.

Author: Evan K. Friis, UW

Takes as input a set of ROOT [files] with analysis histgrams, and the
corresponding lumicalc.sum [lumifiles] that hve the effective lumi for each
sample.

If [blind] is true, data in the p1p2p3 region will not be plotted.

'''

import rootpy.plotting.views as views
from FinalStateAnalysis.PlotTools.Plotter        import Plotter
from FinalStateAnalysis.PlotTools.BlindView      import BlindView
from FinalStateAnalysis.PlotTools.PoissonView    import PoissonView
from FinalStateAnalysis.PlotTools.MedianView     import MedianView
from FinalStateAnalysis.PlotTools.ProjectionView import ProjectionView
from FinalStateAnalysis.PlotTools.FixedIntegralView import FixedIntegralView
from FinalStateAnalysis.PlotTools.DifferentialView import DifferentialView
from FinalStateAnalysis.PlotTools.SubtractionView import SubtractionView, PositiveView
from FinalStateAnalysis.PlotTools.RebinView  import RebinView
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors
from FinalStateAnalysis.PlotTools.decorators import memo
from FinalStateAnalysis.MetaData.datacommon  import br_w_leptons, br_z_leptons
from FinalStateAnalysis.Utilities.floatformatting import smart_float_format
from pdb import set_trace
from optparse import OptionParser
import os
import ROOT
import glob
import math
import logging
from fnmatch import fnmatch

ROOT.gStyle.SetOptTitle(0)
#ROOT.gROOT.SetStyle('Plain')

parser = OptionParser(description=__doc__)
parser.add_option('--dry-run', action='store_true', dest='dry_run', default = False,
                  help='produces only shape file and minimal histograms')

parser.add_option('--no-mc-data', action='store_true', dest='no_mc_data', default = False)
parser.add_option('--no-wz', action='store_true', dest='no_wz', default = False)
parser.add_option('--no-signal', action='store_true', dest='no_signal', default = False)
parser.add_option('--no-f3', action='store_true', dest='no_f3', default = False)
parser.add_option('--no-shapes', action='store_true', dest='no_shapes', default = False)

parser.add_option('--prefix', metavar='label', type=str, dest='prefix', default = '',
                  help='prefix to eppend before histogram name o be used to make the shapes' )
parser.add_option('--prefixes', metavar='label', type=str, dest='prefixes', default = '',
                  help='prefix to eppend before histogram name o be used to make the shapes' )

project_x = lambda x: ProjectionView(x, 'X', [0, 650])

def get_chi_square(hdata, hexp):
    chi2  = 0.
    nbins = 0.
    for i in xrange(1, hdata.GetNbinsX()+1):
        data  = hdata.GetBinContent(i)
        edata = hdata.GetBinError(i)
        exp   = hexp.GetBinContent(i)
        eexp  = hexp.GetBinError(i)
        if data > 0:
            chi2  += (data - exp)**2/(edata**2 + eexp**2)
            nbins += 1
    return chi2, nbins
        

def quad(*xs):
    return math.sqrt(sum(x * x for x in xs))

def create_mapper(mapping):
    def _f(path):
        for key, out in mapping.iteritems():
            if key == path:
                path = path.replace(key,out)
        return path
    return _f

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )

def fake_rate_estimate(histograms): #always, ALWAYS give as 1,2,0
    cat1 = histograms.next()
    cat2 = histograms.next()
    cat0 = histograms.next()
    ret  = cat0.Clone()
    ret.Reset()
    integral = sum( i.Integral()   for i in [cat0, cat1, cat2])
    entries  = sum( i.GetEntries() for i in [cat0, cat1, cat2])
    weight   = integral / float(entries)

    for i in range(ret.GetNbinsX() + 2):
        c1 = cat1.GetBinContent(i)
        c2 = cat2.GetBinContent(i)
        c0 = cat0.GetBinContent(i)
        e1 = cat1.GetBinError(i)
        e2 = cat2.GetBinError(i)
        e0 = cat0.GetBinError(i)
        content = c1+c2-c0
        #if content <= 0:
        #    print "bin [%.1f, %.1f]" % (cat0.GetBinLowEdge(i), cat0.GetBinLowEdge(i)+cat0.GetBinWidth(i))
        #    print "content %s c0: %s c1: %s c2: %s" % (content, c0, c1, c2)
        if content <= 0 and c0 > 0:
            #print "setting content to %s +/- %s" % (c0, e0)
            ret.SetBinContent(i,c0)
            ret.SetBinError(i,e0)
        elif content <= 0:
            #print "setting content to 10**-5 +/- %s" % (1.8*weight)
            ret.SetBinContent(i,10**-5)
            ret.SetBinError(i,1.8*weight)
        else:
            ret.SetBinContent(i,c1+c2-c0)
            ret.SetBinError(i,quad(e1,e2,e0))
    return ret

def make_empty_bin_remover(weight):
    def mc_empty_bin_remover(histogram):
        ret  = histogram.Clone()
        for i in range(ret.GetNbinsX() + 2):
            if ret.GetBinContent(i) <= 0:
                ret.SetBinContent(i, 0.9200*weight) #MEAN WEIGHT
                ret.SetBinError(i, 1.8*weight)
        return ret
    return mc_empty_bin_remover

class MultyErrorView(object): #FIXME
    '''takes a StackView and adds systematics'''
    def __init__(self, stackview, systematics_map):
        self.stack_view = stackview
        self.systematics_map = systematics_map

    def Get(self, path):
        stack_components = self.stack_view.Get(path).GetHists()
        stack_sum = sum(stack_components)
        
        for bin in range(1, stack_sum.GetNbinsX() + 1):
            #sys_errors = [ sum(hist.GetBinContent(bin)*value for hist in stack_components if fnmatch(hist.GetTitle(), key)) for key, value in self.systematics_map ]
            sys_errors = []
            for key, value in self.systematics_map.iteritems():
                it = 0
                for hist in stack_components:
                    if fnmatch(hist.GetTitle(), key):
                        print '%s matches %s' % (hist.GetTitle(), key)
                        it += hist.GetBinContent(bin)*value
                sys_errors.append(it)

            error = quad( stack_sum.GetBinError(bin), *sys_errors )
            stack_sum.SetBinError(bin, error)

        stack_sum.SetMarkerSize(0)
        stack_sum.SetFillColor(1)
        stack_sum.SetFillStyle(3013)
        stack_sum.legendstyle = 'f'
        return stack_sum


class BackgroundErrorView(object):
    ''' Compute the total background error in each bin. '''
    def __init__(self, fakes, wz, zz, charge_fake, wz_error=0.1, zz_error=0.04,
                 fake_error=0.):
        self.fakes = fakes
        self.wz = wz
        self.zz = zz
        self.fake_error = fake_error
        self.wz_error = wz_error
        self.zz_error = zz_error
        self.charge_fake = charge_fake

    def Get(self, path):
        fakes = self.fakes.Get(path)
        wz = self.wz.Get(path)
        zz = self.zz.Get(path)
        charge_fake = self.charge_fake.Get(path)

        bkg_error = wz.Clone()
        bkg_error.SetTitle("Bkg. Unc.")
        bkg_error.Reset()
        for bin in range(1, bkg_error.GetNbinsX() + 1):
            error = quad(
                fakes.GetBinError(bin),
                fakes.GetBinContent(bin) * self.fake_error,
                wz.GetBinContent(bin) * self.wz_error,
                zz.GetBinContent(bin) * self.zz_error,
                charge_fake.GetBinError(bin)
            )
            total = (
                fakes.GetBinContent(bin) +
                wz.GetBinContent(bin) +
                zz.GetBinContent(bin) +
                charge_fake.GetBinContent(bin)
            )
            bkg_error.SetBinContent(bin, total)
            bkg_error.SetBinError(bin, error)
        bkg_error.SetMarkerSize(0)
        bkg_error.SetFillColor(1)
        bkg_error.SetFillStyle(3013)
        bkg_error.legendstyle = 'f'
        return bkg_error


class QCDCorrectionView(object):
    ''' Get the estimate of the fake rate, correcting for QCD contamination

    The amount of QCD contamination "q" in each control region (cr) bin
    is cr_qcd_est/cr.

    The final estimate is then:

        q*cr_qcd_weight + (1-q)*cr_ewk_weight

    '''
    def __init__(self, data_view, cr, cr_qcd_est, cr_ewk_weight,
                 cr_qcd_weight):
        self.cr = views.SubdirectoryView(data_view, cr)
        self.cr_qcd_est = views.SubdirectoryView(data_view, cr_qcd_est)
        self.cr_ewk_weight = views.SubdirectoryView(data_view, cr_ewk_weight)
        self.cr_qcd_weight = views.SubdirectoryView(data_view, cr_qcd_weight)

    def Get(self, path):
        ''' Get the fake est histo corrected for the QCD effect '''

        # Compute the fraction of QCD in each bin
        qcd_fraction = self.cr_qcd_est.Get(path).Clone()
        qcd_fraction.Divide(self.cr.Get(path))
        # A histogram where each bin is initially 1.0
        ewk_fraction = self.cr_qcd_est.Get(path).Clone()
        ewk_fraction.Reset()
        # Make sure each bin is <= 1
        for bin in range(0, qcd_fraction.GetNbinsX() + 1):
            ewk_fraction.SetBinContent(bin, 1)
            ewk_fraction.SetBinError(bin, 0)
            #print bin, qcd_fraction.GetBinContent(bin)
            if qcd_fraction.GetBinContent(bin) > 1:
                qcd_fraction.SetBinContent(bin, 1)

        ewk_fraction.Sumw2()
        ewk_fraction.Add(qcd_fraction, -1)
        for bin in range(0, qcd_fraction.GetNbinsX() + 1):
            pass
            #print ewk_fraction.GetBinContent(bin)
        print "QCD correction in %s" % path

        cr_ewk_weight_corrected = self.cr_ewk_weight.Get(path).Clone()
        print "ewk", cr_ewk_weight_corrected.Integral()
        cr_ewk_weight_corrected.Multiply(ewk_fraction)

        cr_qcd_weight_corrected = self.cr_qcd_weight.Get(path).Clone()
        print "qcd", cr_qcd_weight_corrected.Integral()
        cr_qcd_weight_corrected.Multiply(qcd_fraction)

        total = cr_ewk_weight_corrected.Clone()
        total.Add(cr_qcd_weight_corrected)
        print "total", total.Integral()

        return total


def make_styler(color, format=None):
    def unsuck(x):
        x.SetFillStyle(0)
        x.SetLineColor(color)
        x.SetLineWidth(2)
        x.SetMaximum(1.5*x.GetMaximum())
        if format:
            x.drawstyle = format
    return unsuck


class WHPlotterBase(Plotter):
    def __init__(self, channel, obj1_charge_mapper={}, obj2_charge_mapper={}):
        cwd = os.getcwd()
        #os.chdir( os.path.dirname(__file__) )
        self.channel = channel
        jobid = os.environ['jobid']
        print "\nPlotting %s for %s\n" % (channel, jobid)

        # Figure out if we are 7 or 8 TeV
        period = '7TeV' if '7TeV' in jobid else '8TeV'
        self.period = period
        self.sqrts = 7 if '7TeV' in jobid else 8
        samples = [
            'Zjets_M50',
            'WplusJets_madgraph',
            'WZJetsTo3LNu*',
            'ZZ*',
            'VH*',
            'WW*',
            #'TTplusJets_madgraph',
            "data_*",
        ]

        files = []
        lumifiles = []

        for x in samples:
            files.extend(glob.glob('results/%s/WHAnalyze%s/%s.root' % (jobid, channel, x)))
            lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))

        self.outputdir = 'results/%s/plots/%s' % (jobid, channel.lower())
        self.base_out_dir = self.outputdir
        if not os.path.exists(self.outputdir):
            os.makedirs(self.outputdir)

        blinder = None
        blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
        print '\n\nRunning Blind: %s\n\n' % blind
        self.blind = blind
        self.obj1_charge_mapper=obj1_charge_mapper #maps scaled charge flip histograms to the ones of the normal categories
        self.obj2_charge_mapper=obj2_charge_mapper #maps scaled charge flip histograms to the ones of the normal categories
        if blind:
            # Don't look at the SS all pass region
            blinder = lambda x: BlindView(x, "ss/tau_os/p1p2p3/.*")
        super(WHPlotterBase, self).__init__(files, lumifiles, self.outputdir,
                                            blinder)
        self.defaults = {} #allows to set some options and avoid repeating them each function call
        self.mc_samples = filter(lambda x: not x.startswith('data'), samples)
        #os.chdir( cwd )
        #create a fake wiew summing up all HWW
        #self.views['VH_hww_sum'] = {
        #    'unweighted_view' : views.SumView(
        #        *[item['unweighted_view'] for name, item in self.views.iteritems() if fnmatch(name, 'VH_*_HWW*')]
        #    )
        #}

    def set_subdir(self, folder):
        self.outputdir = '/'.join([self.base_out_dir, folder])

    def apply_to_dict(self, dictionary, viewtype, *args, **kwargs): #project, project_axis, rebin):
        ret = {}
        for key, val in dictionary.iteritems():
            if isinstance(val, dict):
                ret[key] = self.apply_to_dict(val, viewtype, *args, **kwargs) #project, project_axis, rebin)
            elif val is None:
                ret[key] = None
            else:
                ret[key] = viewtype(val, *args, **kwargs) # self.rebin_view( ProjectionView(val, project_axis, project), rebin )
        return ret

    def make_additional_fakes_view(self, unblinded=False, rebin=1, 
                                   project=None, project_axis=None, tau_charge='tau_os' ):
        other_tau_sign = 'tau_os' if tau_charge == 'tau_ss' else 'tau_ss'
        def preprocess(view):
            ret = view
            if project and project_axis:
                ret = ProjectionView(ret, project_axis, project)
            return RebinView( ret, rebin )

        zjets_view    = preprocess( self.get_view('Zjets*') )
        all_data_view = preprocess( self.get_view('data') )
        zjets_fakes   = views.SumView(
            views.SubdirectoryView( zjets_view, 'ss/%s/f1p2p3/w1/' % tau_charge ),
            views.SubdirectoryView( zjets_view, 'ss/%s/p1f2p3/w2/' % tau_charge )
            #,views.SubdirectoryView( zjets_view, 'ss/%s/f1f2p3/w12/' % tau_charge ) #FIXME: check that effect is small
        )
        tau_fakes  = views.SubdirectoryView(all_data_view, 'ss/%s/p1p2f3/w3/' % tau_charge)
        full_fakes = views.SumView(tau_fakes, zjets_fakes)

        return {
            'fakes'      : full_fakes,
            'zjet_fakes' : zjets_fakes,
            'tau_fakes'  : tau_fakes,
        }

    def make_signal_views(self, unblinded=False, qcd_weight_fraction=0, #MARK
                          rebin=1, project=None, project_axis=None, tau_charge='tau_os' ):
        ''' Make signal views with FR background estimation '''
        other_tau_sign = 'tau_os' if tau_charge == 'tau_ss' else 'tau_ss'

        def preprocess(view):
            ret = view
            if project and project_axis:
                ret = ProjectionView(ret, project_axis, project)
            return RebinView( ret, rebin )

        all_wz_view_tautau = preprocess( self.get_view('WZJetsTo3LNu*ZToTauTau*') )
        wz_view_tautau = views.SubdirectoryView(
                all_wz_view_tautau,
                'ss/%s/p1p2p3/' % tau_charge
        )

        tomatch = 'WZJetsTo3LNu' if self.sqrts == 7 else 'WZJetsTo3LNu_pythia'
        all_wz_view_3l = preprocess( self.get_view(tomatch) )
        wz_view_3l = views.SubdirectoryView(
                all_wz_view_3l,
                'ss/%s/p1p2p3/' % tau_charge
        )
        all_wz_view = views.SumView(all_wz_view_tautau, all_wz_view_3l)

        zz_view = preprocess(
            views.SubdirectoryView(
                self.get_view('ZZJetsTo4L*'),
                'ss/%s/p1p2p3/' % tau_charge
            )
        )

        all_data_view = self.get_view('data')
        #if unblinded:
        #    all_data_view = self.get_view('data', 'unblinded_view')
        all_data_view = preprocess(all_data_view)
            
        data_view = views.SubdirectoryView(all_data_view, 'ss/%s/p1p2p3/' % tau_charge)

        def make_fakes(qcd_fraction):
            def make_fakes_view(weight_type, scale):
                scaled_bare_data = views.ScaleView(all_data_view, scale)
                scaled_wz_data   = views.ScaleView(all_wz_view, scale)
                scaled_data      = SubtractionView(scaled_bare_data, scaled_wz_data, restrict_positive=True)

                # View of weighted obj1-fails data
                obj1_view = views.SubdirectoryView(
                    scaled_data, 'ss/%s/f1p2p3/%s1' % (tau_charge, weight_type))
                # View of weighted obj2-fails data
                obj2_view = views.SubdirectoryView(
                    scaled_data, 'ss/%s/p1f2p3/%s2' % (tau_charge, weight_type))
                # View of weighted obj1&2-fails data
                obj12_view = views.SubdirectoryView(
                    scaled_data, 'ss/%s/f1f2p3/%s12' % (tau_charge, weight_type))

                # Give the individual object views nice colors
                obj1_view = views.TitleView(
                    views.StyleView(obj1_view, **remove_name_entry(data_styles['TT*'])),
                    'Reducible bkg. 1')
                obj2_view = views.TitleView(
                    views.StyleView(obj2_view, **remove_name_entry(data_styles['QCD*'])),
                    'Reducible bkg. 2')
                obj12_view = views.TitleView(
                    views.StyleView(obj12_view, **remove_name_entry(data_styles['WW*'])),
                    'Reducible bkg. 12')

                return obj1_view, obj2_view, obj12_view

            qcd1, qcd2, qcd12 = make_fakes_view('q', qcd_fraction)
            wjet1, wjet2, wjet12 = make_fakes_view(
                'w', 1 - qcd_fraction)

            obj1_view = views.SumView(qcd1, wjet1)
            obj2_view = views.SumView(qcd2, wjet2)
            obj12_view = views.SumView(qcd12, wjet12)

            # Corrected fake view
            fakes_view = views.MultiFunctorView(fake_rate_estimate, obj1_view, obj2_view, obj12_view)
            fakes_view = views.TitleView(
                views.StyleView(fakes_view, **remove_name_entry(data_styles['Zjets*'])), 'Reducible bkg.')
            return obj1_view, obj2_view, obj12_view, fakes_view

        obj1_view, obj2_view, obj12_view, fakes_view = make_fakes(qcd_weight_fraction)
        fakes_view_05 = make_fakes(0.5)[-1]
        fakes_view_0  = make_fakes(0)[-1]
        fakes_view_1  = make_fakes(1)[-1]
 
        charge_fakes = views.TitleView( 
            views.StyleView(
                views.SumView(
                    views.PathModifierView(
                        views.SubdirectoryView(all_data_view, 'os/%s/p1p2p3/c1' % other_tau_sign), #FIXME: needs to be fixed for charge 3 region
                        create_mapper(self.obj1_charge_mapper)
                        ),
                    views.PathModifierView(
                        views.SubdirectoryView(all_data_view, 'os/%s/p1p2p3/c2' % tau_charge),
                        create_mapper(self.obj2_charge_mapper)
                        ),
                    ),
                **remove_name_entry(data_styles['TT*'])),
            'Charge mis-id')

        charge_fakes_sysup = views.TitleView( 
            views.StyleView(
                views.SumView(
                    views.PathModifierView(
                        views.SubdirectoryView(all_data_view, 'os/%s/p1p2p3/c1_sysup' % other_tau_sign),
                        create_mapper(self.obj1_charge_mapper)
                        ),
                    views.PathModifierView(
                        views.SubdirectoryView(all_data_view, 'os/%s/p1p2p3/c2_sysup' % tau_charge),
                        create_mapper(self.obj2_charge_mapper)
                        ),
                    ),
                **remove_name_entry(data_styles['TT*'])),
            'Charge mis-id')

        #charge_fakes = MedianView(highv=charge_fakes_sysup, centv=charge_fakes)

        output = {
            'wz': wz_view_tautau,
            'wz_3l': wz_view_3l,
            'zz': zz_view,
            'data': data_view,
            'obj1': obj1_view,
            'obj2': obj2_view,
            'obj12': obj12_view,
            'fakes': fakes_view,
            'weighted_fakes' : {
                0.  : fakes_view_0,
                0.5 : fakes_view_05,
                1.  : fakes_view_1,
            },
            'charge_fakes': {
                'central' : charge_fakes,
                'sys_up'  : charge_fakes_sysup,
                }
        }

        # Add signal
        data_total_lumi = self.views['data']['intlumi']
        for mass in range(90, 165, 5):
            try:
                vh_base_name =  'VH_%s' % mass if self.sqrts == 7 else 'VH_H2Tau_M-%s' % mass
                vh_base_name =  'VHtautau_lepdecay_%s' % mass \
                                if 'VHtautau_lepdecay_%s' % mass in self.views else \
                                vh_base_name
                #print 'using %s' %  vh_base_name
                vh_base_view =  self.views[vh_base_name]['view'] 

                vh_view = views.SubdirectoryView(
                    vh_base_view, #self.get_view('VH_*%i' % mass),
                    'ss/%s/p1p2p3/' % tau_charge 
                )
                output['vh%i' % mass] = preprocess(vh_view)
            except KeyError:
                #logging.warning('No sample found matching VH_*%i' % mass)
                continue

            if mass % 10 == 0 and mass < 150:
                # Only have 10 GeV steps for WW
                try:
                    ww_view = views.SubdirectoryView(
                        self.get_view('VH_%i_HWW*' % mass),
                        'ss/%s/p1p2p3/' % tau_charge
                    )
                    output['vh%i_hww' % mass] = preprocess(ww_view)
                except KeyError:
                    #logging.warning('No sample found matching VH_%i_HWW*' % mass)
                    continue
                #output['signal%i' % mass] = views.SumView(ww_view, vh_view) if ww_view else vh_view

        return output

    def make_qcd_proj_views(self, control_region, rebin):
        ''' Make views when obj1 or obj2 fails, projecting QCD in

        QCD comes from the triple fake region
        '''
        all_data_view = self.rebin_view(self.get_view('data'), rebin)

        mapping = {
            1: {
                'obs': 'ss/tau_os/f1p2p3',
                'qcd': 'ss/tau_os/f1f2f3/w23',
            },
            2: {
                'obs': 'ss/tau_os/p1f2p3',
                'qcd': 'ss/tau_os/f1f2f3/w13',
            },
        }

        data_view = views.TitleView(views.SubdirectoryView(
            all_data_view, mapping[control_region]['obs']),
            "Anti-iso obj %i" % control_region)

        qcd_view = views.TitleView(views.SubdirectoryView(
            all_data_view, mapping[control_region]['qcd']), "QCD")

        qcd_view = views.StyleView(
            qcd_view, drawstyle='hist', linecolor=colors['red'],
            fillstyle=0, legendstyle='l'
        )
        return {'obs': data_view, 'qcd': qcd_view}

    @memo
    def make_obj3_fail_cr_views(self, qcd_correction=False,
                                qcd_weight_fraction=0, tau_charge='tau_os'):
        ''' Make views when obj3 fails, estimating the bkg in obj1 pass using
            f1p2f3 '''
        other_tau_sign = 'tau_os' if tau_charge == 'tau_ss' else 'tau_ss'

        all_wz_ztt_view = self.get_view('WZJetsTo3LNu*ZToTauTau*')
        wz_view = views.SubdirectoryView(
            all_wz_ztt_view,
            'ss/%s/p1p2f3/' % tau_charge
        )

        tomatch = 'WZJetsTo3LNu' if self.sqrts == 7 else 'WZJetsTo3LNu_pythia'
        all_wz_3l_view = self.get_view(tomatch)
        wz_view_3l = views.SubdirectoryView(
            all_wz_3l_view,
            'ss/%s/p1p2f3/' % tau_charge
        )

        all_wz_view = views.SumView(all_wz_ztt_view, all_wz_3l_view)

        zz_view = views.SubdirectoryView(
            self.get_view('ZZJetsTo4L*'),
            'ss/%s/p1p2f3/' % tau_charge
        )
        all_data_view = self.get_view('data')
        data_view = views.SubdirectoryView(all_data_view, 'ss/%s/p1p2f3/' % tau_charge)

        def make_fakes(qcd_fraction):
            def make_fakes_view(weight_type, scale):
                scaled_bare_data = views.ScaleView(all_data_view, scale)
                scaled_wz_data   = views.ScaleView(all_wz_view,   scale)
                scaled_data      = SubtractionView(scaled_bare_data, scaled_wz_data, restrict_positive=True) 
                # View of weighted obj1-fails data
                obj1_view = views.SubdirectoryView(
                    scaled_data, 'ss/%s/f1p2f3/%s1' % (tau_charge, weight_type))
                # View of weighted obj2-fails data
                obj2_view = views.SubdirectoryView(
                    scaled_data, 'ss/%s/p1f2f3/%s2' % (tau_charge, weight_type))
                # View of weighted obj1&2-fails data
                obj12_view = views.SubdirectoryView(
                    scaled_data, 'ss/%s/f1f2f3/%s12' % (tau_charge, weight_type))

                # Give the individual object views nice colors
                obj1_view = views.TitleView(
                    views.StyleView(obj1_view, **remove_name_entry(data_styles['TT*'])),
                    'Reducible bkg. 1')
                obj2_view = views.TitleView(
                    views.StyleView(obj2_view, **remove_name_entry(data_styles['QCD*'])),
                    'Reducible bkg. 2')
                obj12_view = views.TitleView(
                    views.StyleView(obj12_view, **remove_name_entry(data_styles['WW*'])),
                    'Reducible bkg. 12')

                subtract_obj12_view = views.ScaleView(obj12_view, -1)
                return obj1_view, obj2_view, obj12_view, subtract_obj12_view

            qcd1, qcd2, qcd12, negqcd12 = make_fakes_view('q', qcd_fraction)
            wjet1, wjet2, wjet12, negwjet12 = make_fakes_view(
                'w', 1 - qcd_fraction)

            obj1_view  = views.SumView(qcd1, wjet1)
            obj2_view  = views.SumView(qcd2, wjet2)
            obj12_view = views.SumView(qcd12, wjet12)
            subtract_obj12_view = views.SumView(negqcd12, negwjet12)

            # Corrected fake view
            fakes_view = views.SumView(obj1_view, obj2_view, subtract_obj12_view)
            fakes_view = views.TitleView(
                views.StyleView(fakes_view, **remove_name_entry(data_styles['Zjets*'])), 'Reducible bkg.')
            return obj1_view, obj2_view, obj12_view, fakes_view

        obj1_view, obj2_view, obj12_view, fakes_view = make_fakes(qcd_weight_fraction)
        fakes_view_05 = make_fakes(0.5)[-1]
        fakes_view_0  = make_fakes(0)[-1]
        fakes_view_1  = make_fakes(1)[-1]

        style_dict_no_name = remove_name_entry(data_styles['TT*'])
        charge_fakes = views.TitleView(  #FIXME
            views.StyleView(
                views.SumView(
                    views.PathModifierView(
                        views.SubdirectoryView(all_data_view, 'os/%s/p1p2f3/c1' % other_tau_sign),
                        create_mapper(self.obj1_charge_mapper)
                        ),
                    views.PathModifierView(
                        views.SubdirectoryView(all_data_view, 'os/%s/p1p2f3/c2' % tau_charge),
                        create_mapper(self.obj2_charge_mapper)
                        ),
                    ),
                **style_dict_no_name),
            'Charge mis-id')
        
        charge_fakes_sysup = views.TitleView(
                 views.StyleView(
                     views.SumView(
                         views.PathModifierView(
                             views.SubdirectoryView(all_data_view, 'os/%s/p1p2f3/c1_sysup'  % other_tau_sign),
                             create_mapper(self.obj1_charge_mapper)
                         ),
                         views.PathModifierView(
                             views.SubdirectoryView(all_data_view, 'os/%s/p1p2f3/c2_sysup' % tau_charge),
                             create_mapper(self.obj2_charge_mapper)
                         ),
                     ),
                 **style_dict_no_name),
            'Charge mis-id')

        #charge_fakes = MedianView(highv=charge_fakes_sysup, centv=charge_fakes)
                            
        output = {
            'wz': wz_view,
            'wz_3l': wz_view_3l,
            'zz': zz_view,
            'data': data_view,
            'obj1': obj1_view,
            'obj2': obj2_view,
            'obj12': obj12_view,
            'fakes': fakes_view,
            'weighted_fakes' : {
                0.  : fakes_view_0,
                0.5 : fakes_view_05,
                1.  : fakes_view_1,
            },
            'charge_fakes': {
                'central' : charge_fakes,
                'sys_up'  : charge_fakes_sysup,
                }
        }

        #Add signal @ 120, just mo make bkg fitting easier
        mass = 120
        vh_view = views.SubdirectoryView(
            self.get_view('VH_*%i' % mass),
            'ss/tau_os/p1p2f3/'
        )
        output['vh%i' % mass] = vh_view
        try:
            ww_view = views.SubdirectoryView(
                self.get_view('VH_%i_HWW*' % mass),
                'ss/tau_os/p1p2f3/'
            )
        except KeyError:
            #logging.warning('No sample found matching VH_%i_HWW*' % mass)
            ww_view = None
        output['vh%i_hww' % mass] = ww_view
        output['signal%i' % mass] = views.SumView(ww_view, vh_view) if ww_view else vh_view

        return output

    def make_wz_cr_views(self, rebin=1, project=None, project_axis=None):
        ''' Make WZ control region views with FR background estimation '''

        def preprocess(view):
            ret = view
            if project and project_axis:
                ret = ProjectionView(ret, project_axis, project)
            return RebinView( ret, rebin )

        wz_view_tautau_all = preprocess( self.get_view('WZJetsTo3LNu*ZToTauTau*') )
        wz_view_tautau     = views.SubdirectoryView(wz_view_tautau_all, 'ss/tau_os/p1p2p3_enhance_wz/')

        tomatch = 'WZJetsTo3LNu' if self.sqrts == 7 else 'WZJetsTo3LNu_pythia'
        wz_view_3l_all = preprocess( self.get_view(tomatch) )
        wz_view_3l     = views.SubdirectoryView(wz_view_3l_all, 'ss/tau_os/p1p2p3_enhance_wz/')
        wz_view_all    = views.SumView(wz_view_tautau_all, wz_view_3l_all)

        zz_view_all = preprocess( self.get_view('ZZJetsTo4L*') )
        zz_view     = views.SubdirectoryView(zz_view_all, 'ss/tau_os/p1p2p3_enhance_wz/')
            
        all_data_view = preprocess( self.get_view('data') )
        data_view = views.SubdirectoryView(
            all_data_view, 'ss/tau_os/p1p2p3_enhance_wz/')

        # View of weighted obj2-fails data
        fakes_view = views.SubdirectoryView(
            all_data_view, 'ss/tau_os/p1f2p3_enhance_wz/w2')
        fakes_view = views.StyleView(fakes_view, **remove_name_entry(data_styles['Zjets*']))

        # Correct
        wz_in_fakes_view = views.SubdirectoryView(wz_view_all, 'ss/tau_os/p1f2p3_enhance_wz/w2')
        zz_in_fakes_view = views.SubdirectoryView(zz_view_all, 'ss/tau_os/p1f2p3_enhance_wz/w2')

        diboson_view = views.SumView(wz_in_fakes_view, zz_in_fakes_view)
        inverted_diboson_view = views.ScaleView(diboson_view, -1)

        fakes_view = views.SumView(fakes_view, inverted_diboson_view)
        fakes_view = views.TitleView(fakes_view, 'Reducible bkg.')

        output = {
            'wz_ztt': wz_view_tautau,
            'wz_3l' : wz_view_3l,
            'zz'    : zz_view,
            'data'  : data_view,
            'fakes' : fakes_view
        }
        return output
        # Add signal
        #for mass in [110, 120, 130, 140]:
        #    vh_view = views.SubdirectoryView(
        #        self.rebin_view(self.get_view('VH_*%i' % mass, 'unweighted_view'), rebin),
        #        'ss/tau_os/p1p2p3/'
        #    )
        #    output['vh%i' % mass] = vh_view
        #
        #return output

    def write_shapes(self, variable, rebin, outdir,
                     qcd_fraction=0., #[1., 0., -1.], 
                     show_charge_fakes=False,
                     project=None, project_axis=None, different_fakes=False):
        ''' Write final shapes for [variable] into a TDirectory [outputdir] '''
        show_charge_fakes = show_charge_fakes if 'show_charge_fakes' not in self.defaults else self.defaults['show_charge_fakes']
        sig_view = self.make_signal_views(unblinded=(not self.blind),
                                          qcd_weight_fraction=qcd_fraction,
                                          rebin=rebin, project=project, project_axis=project_axis)

        different_fakes_views = self.make_additional_fakes_view( unblinded=(not self.blind), rebin=rebin, 
                                   project=project, project_axis=project_axis)
        
        
        outdir.cd()
        wz_weight = self.get_view('WZJetsTo3LNu*ZToTauTau*', 'weight')
        zz_weight = self.get_view('ZZJetsTo4L*', 'weight')
        print "wz_weight: %s" % wz_weight
        print "zz_weight: %s" % zz_weight

        wz = views.FunctorView( views.SumView(sig_view['wz'], sig_view['wz_3l']), make_empty_bin_remover(wz_weight)).Get(variable)
        zz = views.FunctorView( sig_view['zz'], make_empty_bin_remover(zz_weight)).Get(variable)
        obs = sig_view['data'].Get(variable)
        fakes = sig_view['fakes'].Get(variable) if not different_fakes else different_fakes_views['fakes'].Get(variable)

        fakes_down = different_fakes_views['fakes'].Get(variable)
        fakes_up   = PositiveView(
            views.SumView( views.ScaleView(sig_view['fakes'], 2.), 
                           views.ScaleView(different_fakes_views['fakes'], -1.) 
                       )
        ).Get(variable)

        wz.SetName('wz')
        zz.SetName('zz')
        obs.SetName('data_obs')
        fakes.SetName('fakes')
        fakes_down.SetName('fakes_CMS_vhtt_%s_fakeshape_%sTeVDown' % (outdir.GetName(), self.sqrts))
        fakes_up.SetName('fakes_CMS_vhtt_%s_fakeshape_%sTeVUp' % (outdir.GetName(), self.sqrts))
        #for mass in [110, 115, 120, 125, 130, 135, 140]:
        #set_trace()
        for mass in range(90, 165, 5):
            try:
                vh = None
                if mass == 90 and self.sqrts == 8:
                    vh = views.ScaleView(sig_view['vh100'], 1.3719).Get(variable)
                elif mass == 95 and self.sqrts == 8:
                    vh = views.ScaleView(sig_view['vh100'], 1.1717).Get(variable)
                else:
                    vh = sig_view['vh%i' % mass].Get(variable)
                vh.SetName('WH%i' % mass)
                vh.SetLineColor(0)
                vh.Write()
            except KeyError:
                #logging.warning('No sample found matching VH_*%i' % mass)
                continue

            if mass % 10 == 0 and mass < 150:
                # Only have 10 GeV steps for WW
                if 'vh%i_hww' % mass in sig_view:
                    ww = sig_view['vh%i_hww' % mass].Get(variable)
                    ww.SetName('WH_hww%i' % mass)
                    ww.Write()

        wz.Write()
        zz.Write()
        obs.Write()
        fakes.Write()
        fakes_down.Write()
        fakes_up.Write()
        #charge_fakes_CMS_vhtt_emt_chargeFlip_8TeVUpx
        if show_charge_fakes:
            logging.info('adding charge fakes shape errors')
            charge_fakes          = sig_view['charge_fakes']['central'].Get(variable)
            charge_fakes_sys_up   = sig_view['charge_fakes']['sys_up' ].Get(variable) #shift='up') 
            charge_fakes_sys_down = charge_fakes+charge_fakes - charge_fakes_sys_up
            charge_fakes.SetName('charge_fakes')
            charge_fakes_sys_up.SetName('charge_fakes_CMS_vhtt_%s_chargeFlip_%sTeVUp' % (self.channel.lower(), self.sqrts))
            charge_fakes_sys_down.SetName('charge_fakes_CMS_vhtt_%s_chargeFlip_%sTeVDown' % (self.channel.lower(), self.sqrts))
            charge_fakes.Write()
            charge_fakes_sys_up.Write()
            charge_fakes_sys_down.Write()
        
    def write_cut_and_count(self, variable, outdir):
        ''' Version of write_shapes(...) with only one bin.

        Equivalent to a cut & count analysis.
        '''
        sig_view = self.make_signal_views( unblinded=(not self.blind))
        nbins = sig_view['wz'].Get(variable).GetNbinsX()
        return self.write_shapes(variable, nbins, outdir, unblinded)

    def write_f3_shapes(self, variable, rebin, outdir,
                     qcd_fraction=0, show_charge_fakes=False,
                     project=None, project_axis=None):
        ''' Write final shapes for [variable] into a TDirectory [outputdir] '''
        show_charge_fakes = show_charge_fakes if 'show_charge_fakes' not in self.defaults else self.defaults['show_charge_fakes']
        sig_view = self.make_obj3_fail_cr_views(False, qcd_weight_fraction=qcd_fraction)
        
        if project and project_axis:
            sig_view = self.apply_to_dict( sig_view, ProjectionView, project_axis, project )
        sig_view = self.apply_to_dict( sig_view, RebinView, rebin )

        outdir.cd()
        wz = views.SumView(sig_view['wz'], sig_view['wz_3l']).Get(variable)
        zz = sig_view['zz'].Get(variable)
        obs = sig_view['data'].Get(variable)
        fakes = sig_view['fakes'].Get(variable)

        wz.SetName('wz')
        zz.SetName('zz')
        obs.SetName('data_obs')
        fakes.SetName('fakes')
        
        wz.Write()
        zz.Write()
        obs.Write()
        fakes.Write()

        mass = 120
        vh = sig_view['vh%i' % mass].Get(variable)
        vh.SetName('WH%i' % mass)
        vh.Write()
        # Only have 10 GeV steps for WW
        if sig_view['vh%i_hww' % mass]:
            ww = sig_view['vh%i_hww' % mass].Get(variable)
            ww.SetName('WH_hww125')
            ww.Write()

        #charge_fakes_CMS_vhtt_emt_chargeFlip_8TeVUpx
        if show_charge_fakes:
            logging.info('adding charge fakes shape errors')
            charge_fakes          = sig_view['charge_fakes']['central'].Get(variable)
            charge_fakes_sys_up   = sig_view['charge_fakes']['sys_up' ].Get(variable) #shift='up') 
            charge_fakes_sys_down = charge_fakes+charge_fakes - charge_fakes_sys_up
            charge_fakes.SetName('charge_fakes')
            charge_fakes_sys_up.SetName('charge_fakes_CMS_vhtt_%s_chargeFlip_%sTeVUp' % (self.channel.lower(), self.sqrts))
            charge_fakes_sys_down.SetName('charge_fakes_CMS_vhtt_%s_chargeFlip_%sTeVDown' % (self.channel.lower(), self.sqrts))
            charge_fakes.Write()
            charge_fakes_sys_up.Write()
            charge_fakes_sys_down.Write()

    def plot_final(self, variable, rebin=1, xaxis='', maxy=24,
                   show_error=False, qcd_correction=False, stack_higgs=True, 
                   qcd_weight_fraction=0., x_range=None, show_charge_fakes=False,
                   leftside_legend=False, higgs_xsec_multiplier=1, project=None, 
                   project_axis=None, differential=False, yaxis='Events', tau_charge='tau_os', **kwargs):
        ''' Plot the final output - with bkg. estimation '''        
        show_charge_fakes = show_charge_fakes if 'show_charge_fakes' not in self.defaults else self.defaults['show_charge_fakes']
        sig_view = self.make_signal_views(unblinded=(not self.blind), 
                                          qcd_weight_fraction=qcd_weight_fraction,
                                          rebin=rebin, project=project, 
                                          project_axis=project_axis, tau_charge=tau_charge)

        if differential:
            sig_view = self.apply_to_dict(sig_view, DifferentialView)

        vh_10x = views.TitleView(
            views.StyleView(
                views.ScaleView(sig_view['vh125'], higgs_xsec_multiplier),
                **remove_name_entry(data_styles['VH*'])
            ),
            "(%i#times) m_{H} = 125" % higgs_xsec_multiplier
        )
        charge_fakes_view = MedianView(highv=sig_view['charge_fakes']['sys_up'], centv=sig_view['charge_fakes']['central'])

        # Fudge factor to go from 120->125 - change in xsec*BR
        #vh_10x = views.ScaleView(vh_10x), .783)
        tostack = [sig_view['wz_3l'], sig_view['zz'], sig_view['wz'], sig_view['fakes'], vh_10x] if stack_higgs else \
            [sig_view['wz_3l'], sig_view['zz'], sig_view['wz'], sig_view['fakes']]
        if show_charge_fakes:
            tostack = tostack[:2]+[charge_fakes_view]+tostack[2:]

        vh_hww = views.ScaleView(sig_view['vh120_hww'], .783) if 'vh120_hww' in sig_view else None
        if vh_hww:
            tostack = tostack[:-1] + [vh_hww] + tostack[-1:]

        stack = views.StackView( *tostack )
        histo = stack.Get(variable)
        
        histo.Draw()
        histo.GetHistogram().GetXaxis().SetTitle(xaxis)
        histo.GetHistogram().GetYaxis().SetTitle(yaxis)

        if x_range:
            histo.GetHistogram().GetXaxis().SetRangeUser(x_range[0], x_range[1])
        self.keep.append(histo)

        # Add legend
        entries = len(tostack)+1
        if show_error:
            entries += 1
        legend  = self.add_legend(histo, leftside=leftside_legend, entries=entries)

        if show_error:
            #correct_qcd_view = None
            #if qcd_weight_fraction == 0:
            #    fakes05 = sig_view['weighted_fakes'][1.]
            #    correct_qcd_view = MedianView(lowv=fakes05, centv=sig_view['fakes'])
            #
            #elif qcd_weight_fraction == 0.5:
            #    fakes1 = sig_view['weighted_fakes'][1.]
            #    correct_qcd_view = MedianView(highv=fakes1, centv=sig_view['fakes'])
            #
            #elif qcd_weight_fraction == 1:
            #    fakes05 = sig_view['weighted_fakes'][0.5]
            #    correct_qcd_view = MedianView(lowv=fakes05, centv=sig_view['fakes'])

            bkg_error_view = BackgroundErrorView(
                sig_view['fakes'], #correct_qcd_view, #sig_view['fakes'],
                views.SumView( sig_view['wz'], sig_view['wz_3l']),
                sig_view['zz'],
                charge_fakes_view,
                fake_error=0.3,
                **kwargs
            )
            bkg_error = bkg_error_view.Get(variable)
            self.keep.append(bkg_error)
            bkg_error.Draw('pe2,same')
            legend.AddEntry(bkg_error)

        # Use poisson error bars on the data
        sig_view['data'] = PoissonView(sig_view['data'], x_err=False, is_scaled=differential) #PoissonView(, x_err=False)

        data = sig_view['data'].Get(variable)
        ymax = histo.GetMaximum()
        if not self.blind or tau_charge != 'tau_os':
            #print "drawing", data.Integral()
            data.Draw('pe,same')
            legend.AddEntry(data)
            ymax = max(ymax, data.GetMaximum())
        self.keep.append(data)

        if isinstance(maxy, (int, long, float)):
            #print "setting maxy to %s" % maxy
            histo.SetMaximum(maxy)
            self.canvas.Update()
        else:
            histo.SetMaximum(ymax*1.2)

        if not stack_higgs:
            higgs_plot = vh_10x.Get(variable)
            higgs_plot.Draw('same')
            self.keep.append(higgs_plot)
            
        legend.Draw()

    def plot_final_wz(self, variable, rebin=1, xaxis='', maxy=None, project=None, project_axis=None):
        ''' Plot the final WZ control region - with bkg. estimation '''
        sig_view = self.make_wz_cr_views(rebin, project, project_axis)
 
        stack = views.StackView(
            sig_view['wz_ztt'],
            sig_view['zz'],
            sig_view['fakes'],
            sig_view['wz_3l']
        )
        histo = stack.Get(variable)
        histo.Draw()
        histo.GetHistogram().GetXaxis().SetTitle(xaxis)
        data = sig_view['data'].Get(variable)
        data.Draw('same')
        if maxy is not None:
            histo.SetMaximum(float(maxy))
        else:
            histo.SetMaximum(1.2 * max(histo.GetMaximum(), data.GetMaximum()))
        self.keep.append(data)
        self.keep.append(histo)

        # Add legend
        self.add_legend(histo, leftside=False, entries=len(histo))

    def plot_final_f3(self, variable, rebin=1, xaxis='', maxy=None,
                      show_error=True, qcd_correction=False,
                      qcd_weight_fraction=0., x_range=None, #):
                      show_chi2=False,project=None, 
                      project_axis=None, differential=False, 
                      yaxis='Events', tau_charge='tau_os', show_ratio=False,
                      ratio_range=None, fit=None, **kwargs):
        ''' Plot the final F3 control region - with bkg. estimation '''
        show_chi2 = False #broken
        sig_view = self.make_obj3_fail_cr_views(
            qcd_correction, qcd_weight_fraction, tau_charge)
        if project and project_axis:
            sig_view = self.apply_to_dict( sig_view, ProjectionView, project_axis, project ) 
        sig_view = self.apply_to_dict( sig_view, RebinView, rebin ) #Rebin
        if differential:
            sig_view = self.apply_to_dict(sig_view, DifferentialView)

        charge_fakes_view = MedianView(highv=sig_view['charge_fakes']['sys_up'], centv=sig_view['charge_fakes']['central'])

        stack = views.StackView(
            sig_view['zz'],
            charge_fakes_view,
            sig_view['wz_3l'],
            sig_view['wz'],
            sig_view['fakes'],
        )
        histo = stack.Get(variable)
        histo.Draw()
        histo.GetHistogram().GetXaxis().SetTitle(xaxis)
        histo.GetHistogram().GetYaxis().SetTitle(yaxis)

        if x_range:
            histo.GetHistogram().GetXaxis().SetRangeUser(x_range[0], x_range[1])

        data = sig_view['data'].Get(variable)

        # Add legend
        legend  = self.add_legend(histo, leftside=False, entries=4)
        #latex   = ROOT.TLatex(0.01, 0.9, "")
        #pad     = ROOT.TPad('da','fuq',0.1,0.8,0.5,0.9)
        #self.canvas.cd()
        #latexit = ''
        bkg_error = None
        if show_error:
            #correct_qcd_view = None
            #if qcd_weight_fraction == 0:
            #    fakes05 = sig_view['weighted_fakes'][1.]
            #    correct_qcd_view = MedianView(lowv=fakes05, centv=sig_view['fakes'])
            #
            #elif qcd_weight_fraction == 0.5:
            #    fakes1 = sig_view['weighted_fakes'][1.]
            #    correct_qcd_view = MedianView(highv=fakes1, centv=sig_view['fakes'])
            #
            #elif qcd_weight_fraction == 1:
            #    fakes05 = sig_view['weighted_fakes'][0.5]
            #    correct_qcd_view = MedianView(lowv=fakes05, centv=sig_view['fakes'])

            bkg_error_view = BackgroundErrorView(
                sig_view['fakes'], #correct_qcd_view, 
                views.SumView(sig_view['wz'], sig_view['wz_3l']),
                sig_view['zz'],
                charge_fakes_view,
                fake_error=0.3,
                **kwargs
            )
            bkg_error = bkg_error_view.Get(variable)
            self.keep.append(bkg_error)
            bkg_error.Draw('pe2,same')
            legend.AddEntry(bkg_error)
            #if show_chi2:
            #    chival  = get_chi_square(data, bkg_error)
            #    latexit = '#chi^{2}/#bins = %.2f / %i' % chival 
                
        data.Draw('same')
        if isinstance(maxy, (int, long, float)):
            histo.SetMaximum(maxy)
        else:
            #histo.SetMaximum(
                #1.2*max(histo.GetHistogram().GetMaximum(), data.GetMaximum()))
            histo.SetMaximum(2 * max(data.GetMaximum(), histo.GetMaximum()))
        self.keep.append(data)
        self.keep.append(histo)
        legend.Draw()
        if show_ratio:
            ratio_plot = self.add_ratio_plot(data, bkg_error, x_range, ratio_range=ratio_range)
            if fit:
                #set_trace()
                fitrange = fit.get('range', False)
                if not fitrange:
                    nbins = ratio_plot.GetNbinsX()
                    fitrange = x_range if x_range else [ ratio_plot.GetBinLowEdge(1), 
                       ratio_plot.GetBinLowEdge(nbins)+ratio_plot.GetBinWidth(nbins)]
                self.lower_pad.cd()
                function = self.fit_shape(ratio_plot, fit['model'], fitrange, fit.get('options','IRMENS'))
                toprint  = '#chi^{2} / DoF = %.2f / %i\n' % (function.GetChisquare() , function.GetNDF())
                for i in range(function.GetNpar()):
                    name  = function.GetParName(i) 
                    value = function.GetParameter(i)
                    error = function.GetParError(i)
                    toprint += '%s = %s\n' % (name, smart_float_format((value, error))) #%f #pm %f
            
                stat_box = self.make_text_box(toprint[:-1],fit.get('stat position','bottom-left'))
                stat_box.Draw()
                self.keep.append(stat_box)
                #print toprint
                self.pad.cd()

    def plot_final_f3_split(self, variable, rebin=1, xaxis='', maxy=None):
        ''' Plot the final F3 control region - with bkg. estimation '''
        sig_view = self.make_obj3_fail_cr_views(
            False, 0.5)
        sig_view = self.apply_to_dict( sig_view, RebinView, rebin ) #Rebin

        stack = views.StackView(
            sig_view['zz'],
            MedianView(highv=sig_view['charge_fakes']['sys_up'], centv=sig_view['charge_fakes']['central']),
            sig_view['wz'],
            sig_view['obj1'],
            sig_view['obj2'],
            sig_view['obj12'],
        )
        histo = stack.Get(variable)
        histo.Draw()
        histo.GetHistogram().GetXaxis().SetTitle(xaxis)

        # Add legend
        legend = self.add_legend(histo, leftside=False, entries=4)

        data = sig_view['data'].Get(variable)
        data.Draw('same')
        if maxy:
            histo.SetMaximum(maxy)
        else:
            histo.SetMaximum(
                1.2 * max(
                    histo.GetHistogram().GetMaximum(), data.GetMaximum()))
        self.keep.append(data)
        self.keep.append(histo)

        #legend.AddEntry(data)
        legend.Draw()

    def plot_qcd_contamination(self, variable, control_region, rebin):
        proj_views = self.make_qcd_proj_views(control_region, rebin)
        data = proj_views['obs'].Get(variable)
        qcd = proj_views['qcd'].Get(variable)

        data.Draw()
        qcd.Draw('same')
        self.keep.append(data)
        self.keep.append(qcd)

        legend = self.add_legend([data, qcd], leftside=False, entries=2)
        legend.Draw()
