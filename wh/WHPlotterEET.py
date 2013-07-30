'''

Plot the EET channel

Usage: python WHPlotterEET.py

'''

import glob
import logging
import os
import ROOT
import sys
import WHPlotterBase
from WHPlotterBase import make_styler, remove_name_entry, parser
import rootpy.plotting.views as views
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors
#from pudb import set_trace; set_trace()


logging.basicConfig(stream=sys.stderr, level=logging.WARNING)

class WHPlotterEET(WHPlotterBase.WHPlotterBase):
    def __init__(self):
        obj1_charge_mapper={"e1_e2_Mass":"e*1_e2_Mass","e1_t_Mass":"e*1_t_Mass"}
        obj2_charge_mapper={
            "e1_e2_Mass" : "e1_e*2_Mass",
            "e2_t_Mass"  : "e*2_t_Mass",
            "e2_t_Mass#faking_prob" : "e*2_t_Mass#faking_prob",
            "e2_t_Mass#log_prob"    : "e*2_t_Mass#log_prob"   ,
            "e2_t_Mass#LT"          : "e*2_t_Mass#LT"         ,
            "e2_t_Mass#tPt"         : "e*2_t_Mass#tPt"        ,}

        super(WHPlotterEET, self).__init__('EET', obj1_charge_mapper, obj2_charge_mapper)

if __name__ == "__main__":
    plotter = WHPlotterEET()
    sqrts   = plotter.sqrts
    plotter.defaults['show_charge_fakes'] = True
    options,NOTUSED = parser.parse_args()
    
    if not options.dry_run:
        ###########################################################################
        ##  Zmm control plots #####################################################
        ###########################################################################
        plotter.set_subdir('mc_data')
        ## # Control Z->mumu + jet region
        ## plotter.plot_mc_vs_data('os/p1p2f3', 'e1_e2_Mass', xaxis='m_{ee} (GeV)', xrange=(60, 120))
        ## plotter.add_cms_blurb(sqrts)
        ## plotter.save('mcdata-os-p1p2f3-e1e2Mass')

        ## plotter.plot_mc_vs_data('os/p1p2f3/w3', 'e1_e2_Mass')
        ## plotter.save('mcdata-os-p1p2f3-w3-e1e2Mass')

        ## plotter.plot_mc_vs_data('os/p1f2p3', 'e1_e2_Mass', xaxis='m_{ee} (GeV)', xrange=(60, 120))
        ## plotter.add_cms_blurb(sqrts)
        ## plotter.save('mcdata-os-p1f2p3-e1e2Mass')

        ## # Check PU variables
        ## plotter.plot_mc_vs_data('os/p1p2f3', 'rho')
        ## plotter.save('mcdata-os-p1p2f3-rho')

        ## plotter.plot_mc_vs_data('os/p1p2f3', 'nvtx')
        ## plotter.save('mcdata-os-p1p2f3-nvtx')

        ## # Lower stat but closer to signal region
        ## plotter.plot_mc_vs_data('os/p1p2p3', 'rho')
        ## plotter.save('mcdata-os-p1p2p3-rho')

        ## plotter.plot_mc_vs_data('os/p1p2p3', 'nvtx')
        ## plotter.save('mcdata-os-p1p2p3-nvtx')

        ## # Make Z->mumu + tau jet control

        ## #
        ## # Makes Tau fr control plot
        ## #
        ## zmm_weighted = plotter.plot('data', 'os/p1p2f3/w3/e1_e2_Mass',  'hist', styler=make_styler(2, 'hist'), xrange=(60, 120))
        ## zmm_weighted.SetTitle("Zee + fake #tau_{h} est.")
        ## zmm_weighted.legendstyle='l'
        ## zmm_weighted.GetXaxis().SetTitle("m_{ee} (GeV)")

        ## zmm_unweighted = plotter.plot('data', 'os/p1p2p3/e1_e2_Mass', 'same', styler=make_styler(1), xrange=(60, 120))
        ## zmm_unweighted.SetTitle("Zee observed")
        ## zmm_unweighted.SetTitle("Zee + fake #tau_{h} obs.")
        ## zmm_unweighted.legendstyle='pe'

        ## plotter.add_legend([zmm_weighted, zmm_unweighted])
        ## plotter.add_cms_blurb(sqrts)
        ## plotter.save('zmm-os-fr-control')

        ## #
        ## # Makes charge fr control plot
        ## #

        ## zeet_os_weighted = plotter.plot('data', 'os/p1p2f3/c1/e1_e2_Mass',  'hist', styler=make_styler(2, 'hist'), xrange=(60, 120))
        ## zeet_os_weighted.SetTitle("Ze^{#pm}e^{#mp} + fake #tau_{h} charge flip est.")
        ## zeet_os_weighted.legendstyle='l'
        ## zeet_os_weighted.GetXaxis().SetTitle("M_{ee} (GeV)")

        ## zee_ss_unweighted = plotter.plot('data', 'ss/p1p2f3/e1_e2_Mass', 'same', styler=make_styler(1), xrange=(60, 120))
        ## zee_ss_unweighted.SetTitle("Ze^{#pm}e^{#pm} + fake #tau_{h} obs.")
        ## zee_ss_unweighted.legendstyle='pe'

        ## plotter.add_legend([zeet_os_weighted, zee_ss_unweighted])
        ## plotter.add_cms_blurb(sqrts)
        ## plotter.save('zee-os-ss-charge-flip-control')

        ## plotter.plot('Zjets_M50', 'os/p1p2f3/weight')
        ## plotter.save('zee-mc-event-weights')
        ## # Check MC weights
        ## ## plotter.plot('Zjets_M50', 'os/p1p2f3/weight_nopu')
        ## ## plotter.save('zee-mc-event-weight_nopu')


        ## ###########################################################################
        ## ##  FR sideband MC-vs-Data ################################################
        ## ###########################################################################

        plotter.plot_mc_vs_data('ss/p1f2p3', 'e2_t_Mass', rebin=10, xaxis='m_{e2#tau} (GeV)', leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-ss-p1f2p3-subMass')

        plotter.plot_mc_vs_data('ss/p1f2p3', 'LT', rebin=10, xaxis='LT (GeV)', leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-ss-p1f2p3-LT')

        ###########################################################################
        ##  Signal region plots    ################################################
        ###########################################################################
        plotter.set_subdir('')
        rebin_slim = [0,20]+range(30, 81, 10)+[100, 200]

        plotter.plot_final('e2_t_Mass#LT', rebin_slim, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, project=[0, 130], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTLow')

        plotter.plot_final('e2_t_Mass#LT', rebin_slim, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTHigh')

        plotter.plot_final('e2_t_Mass#LT', rebin_slim, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, project=[80, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTCut')

        plotter.plot_final('e2_t_Mass#LT', 20, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, project=[80, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTCut-flatbin')

        plotter.plot_final('e1Pt', 10, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-e1Pt')

        plotter.plot_final('e2Pt', 10, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-e2Pt')

        plotter.plot_final('tPt', 10, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-tPt')

        plotter.plot_final('e1AbsEta', 10, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-e1AbsEta')

        plotter.plot_final('e2AbsEta', 10, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-e2AbsEta')

        plotter.plot_final('tAbsEta', 10, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-tAbsEta')

        plotter.plot_final('e2_t_Mass', 20, xaxis='m_{e_{2}#tau} (GeV)', maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass')

        plotter.plot_final('e1_t_Mass', 20, xaxis='m_{e_{1}#tau} (GeV)', maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-e1tMass')

        plotter.plot_final('e1_e2_Mass', 20, xaxis='m_{ee} (GeV)', maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-e1e2Mass')

        plotter.plot_final('tAbsEta', 5, xaxis='|#eta_#tau|', maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-tAbsEta')

        plotter.plot_final('e2RelPFIsoDB', 10, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-e2Iso')

        ## plotter.plot_final('metSignificance', 5)
        ## plotter.add_cms_blurb(sqrts)
        ## plotter.save('final-metSig')

        plotter.plot_final('LT', 10, maxy='auto', xaxis='m_{e_{1}#tau} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-LT')

        plotter.plot_final('e2_t_Mass#LT', 20, xaxis='subleading mass from projection', maxy=None, project=[0, 600], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMassProj')

        plotter.plot_final('e2_t_Mass#LT', 20, xaxis='M_{e_{2}#tau} (GeV)', maxy=None, project=[0, 70], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTLo')

        plotter.plot_final('e2_t_Mass#LT', 20, xaxis='M_{e_{2}#tau} (GeV)', maxy=None, project=[70, 600], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTHi')

        ###########################################################################
        ##  f3 region plots        ################################################
        ###########################################################################
        categories = {
            'LTCut' : [80, 650],
            'LTLow' : [0, 130],
            'LTHigh': [130, 650],
            }
                
        for label, proj_range in categories.iteritems():
            factor = 1.5 if label == 'LTHigh' else 1
            plotter.set_subdir('f3/%s' % label)
            plotter.plot_final_f3('e2_t_Mass#LT', rebin_slim, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-subMass-%s' % label)

            plotter.plot_final_f3('e2_t_Mass#LT', 200, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-subMass-%s-counting' % label, dotc=True, dotroot=True)

            #pt
            plotter.plot_final_f3("e1Pt#LT"    , int(factor*10), xaxis='p_{Te_{1}} (GeV)', maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-e1Pt-%s' % label)

            plotter.plot_final_f3("e2Pt#LT"    , int(factor*10), xaxis='p_{Te_{2}} (GeV)', maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-e2Pt-%s' % label)

            plotter.plot_final_f3("e1JetPt#LT" , int(factor*10), xaxis='p_{T Jet e_{1}} (GeV)', maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-e1JetPt-%s' % label)

            plotter.plot_final_f3("e2JetPt#LT" , int(factor*10), xaxis='p_{T Jet e_{2}} (GeV)', maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-e2JetPt-%s' % label)

            #eta
            plotter.plot_final_f3("e1AbsEta#LT", 10, xaxis='|#eta_{e_{1}}|', maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-e1AbsEta-%s' % label)

            plotter.plot_final_f3("e2AbsEta#LT", 10, xaxis='|#eta_{e_{2}}|'  , maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-e2AbsEta-%s' % label)

            #DR
            plotter.plot_final_f3("e1_t_DR#LT", 10, xaxis='#DeltaR_{e_{1}#tau}', maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-e1_t_DR-%s' % label)

            plotter.plot_final_f3("e2_t_DR#LT", 10, xaxis='#DeltaR_{e_{2}#tau}'  , maxy=None, project=proj_range, project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-e2_t_DR-%s' % label)



        plotter.set_subdir('f3')
        plotter.plot_final_f3('e2_t_Mass', 20, xaxis='m_{e_{2}#tau} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass')

        plotter.plot_final_f3('e2_t_Mass', 200, xaxis='m_{e_{2}#tau} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-counting-like', dotc=True, dotroot=True)

        plotter.plot_final_f3('e2_t_Mass', 20, xaxis='m_{e_{2}#tau} (GeV)', qcd_weight_fraction=0, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-wjetfake-subMass')

        plotter.plot_final_f3('e2_t_Mass', 20, xaxis='m_{e_{2}#tau} (GeV)', qcd_weight_fraction=1, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-qcdfake-subMass')

        plotter.plot_final_f3('e1_t_Mass', 20, xaxis='m_{e_{1}#tau} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e1tMass')

        plotter.plot_final_f3('e1_e2_Mass', 20, xaxis='m_{ee} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e1e2Mass')

        plotter.plot_final_f3('e1Pt', 10, qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e1Pt')

        plotter.plot_final_f3('e2Pt', 10, qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e2Pt')

        plotter.plot_final_f3('e2JetPt', 10, qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e2JetPt')

        plotter.plot_final_f3('tPt', 10, qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-tPt')

        plotter.plot_final_f3('e1AbsEta', 10, qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e1AbsEta')

        plotter.plot_final_f3('e2AbsEta', 10, qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e2AbsEta')

        plotter.plot_final_f3('tAbsEta', 10, qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-tAbsEta')

        plotter.plot_final_f3('LT', 5, qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-LT')

        plotter.plot_final_f3('e1_t_DR', 20, xaxis='#DeltaR_{e_{1}#tau}', x_range=[0,5])
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e1_t_DR')

        plotter.plot_final_f3('e2_t_DR', 20, xaxis='#DeltaR_{e_{2}#tau}', x_range=[0,5])
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-e2_t_DR')

        plotter.plot_final_f3('e2_t_Mass#LT', 20, xaxis='subleading mass from projection', maxy=None, project=[0, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMassProj')

        plotter.plot_final_f3('e2_t_Mass#LT', 20, xaxis='M_{e_{2}#tau} (GeV)', maxy=None, project=[0, 70], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-LTLo')

        plotter.plot_final_f3('e2_t_Mass#LT', 20, xaxis='M_{e_{2}#tau} (GeV)', maxy=None, project=[70, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-LTHi')

        ###########################################################################
        ##  charge flip region plots        #######################################
        ###########################################################################
        plotter.set_subdir('charge_flip_CR_f3')

        plotter.plot_final_f3('charge_flip_CR/e2_t_Mass', 20, xaxis='m_{e_{2}#tau} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-chargeFlip-subMass')

        plotter.plot_final_f3('charge_flip_CR/e1_t_Mass', 20, xaxis='m_{e_{1}#tau} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-chargeFlip-e1tMass')

        plotter.plot_final_f3('charge_flip_CR/e1_e2_Mass', 5, xaxis='m_{ee} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-chargeFlip-e1e2Mass')
        
        #END OF if not options.dry_run:
    ###########################################################################
    ##  Making shape file     #################################################
    ###########################################################################
    plotter.set_subdir('')
    prefixes = [options.prefix+'$'] if options.prefix else ['']
    prefixes = [i+'$' for i in options.prefixes.split(',') if i] if options.prefixes else prefixes
    for prefix in prefixes:
        shape_prefix = prefix if len(prefixes) > 1 else ''
        shape_prefix = shape_prefix.replace(':','_').replace('$','_')

        plotter.plot_final_f3(prefix+'e2_t_Mass', 20, qcd_weight_fraction=0.5, xaxis='m_{e_{2}#tau} (GeV)', show_error=True, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.canvas.SetGridx()
        plotter.canvas.SetGridy()
        plotter.save('final-%s-f3-subMass' % shape_prefix)

        plotter.plot_final(prefix+'e2_t_Mass', 20, qcd_weight_fraction=0.5, xaxis='m_{e_{2}#tau} (GeV)', maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.canvas.SetGridx()
        plotter.canvas.SetGridy()
        plotter.save('final-%s-subMass' % shape_prefix)

        plotter.plot_final(prefix+'LT', 10, qcd_weight_fraction=0.5, xaxis='m_{e_{2}#tau} (GeV)', maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.canvas.SetGridx()
        plotter.canvas.SetGridy()
        plotter.save('final-%s-LT' % shape_prefix)

        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, 'LTCut_eet_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
        shape_dir = shape_file.mkdir('eetCatHigh')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', 20, shape_dir, qcd_fraction=0.5, project=[80, 650], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatHigh_w')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', 20, shape_dir, qcd_fraction=0.0, project=[80, 650], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatHigh_q')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', 20, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')
        logging.warning('shape file %s created' % shape_file.GetName())
        shape_file.Close()


        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, '%seet_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
        
        shape_dir = shape_file.mkdir('eetCatLow')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0, 100], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatLow_w')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0, 100], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatLow_q')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0, 100], project_axis='X')
        
        shape_dir = shape_file.mkdir('eetCatHigh')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[100, 650], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatHigh_w')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[100, 650], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatHigh_q')
        plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[100, 650], project_axis='X')
        
        logging.warning('shape file %s created' % shape_file.GetName()) 
        shape_file.Close()


        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, '%seet_f3_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
        
        shape_dir = shape_file.mkdir('eetCatLow')
        plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatLow_w')
        plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatLow_q')
        plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
        
        shape_dir = shape_file.mkdir('eetCatHigh')
        plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatHigh_w')
        plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
        shape_dir = shape_file.mkdir('eetCatHigh_q')
        plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
        
        logging.warning('shape file %s created' % shape_file.GetName()) 
        shape_file.Close()
