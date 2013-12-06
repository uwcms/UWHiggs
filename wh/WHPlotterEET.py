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
from WHPlotterBase import make_styler, remove_name_entry, parser, project_x
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

rebin_slim = range(20, 81, 10)+[100, 130, 300]
categories = {
    'LTCut' : [80, 650],
    'LTLow' : [0, 100],
    'LTHigh': [100, 650],
    'Full'  : [0,  650],
}

if __name__ == "__main__":
    plotter = WHPlotterEET()
    sqrts   = plotter.sqrts
    plotter.defaults['show_charge_fakes'] = True
    options,NOTUSED = parser.parse_args()
    
    if not options.dry_run:
        if not options.no_mc_data:
            ###########################################################################
            ##  Zee control plots #####################################################
            ###########################################################################
            ## # Control Z->ee + jet region
            plotter.set_subdir('mc_data/zee')
        
            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'e1_e2_Mass', rebin=10, 
                                    xaxis='m_{#mu#mu} (GeV)', leftside=False)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zee-e1e2Mass')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'e1Pt#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zee-e1Pt')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'e2Pt#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zee-e2Pt')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'e1AbsEta#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zee-e1AbsEta')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'e2AbsEta#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zee-e2AbsEta')

            # Check PU variables
            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'rho')
            plotter.save('mcdata-Zee-rho')
        
            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'nvtx')
            plotter.save('mcdata-Zee-nvtx')

        if not options.no_wz:
            ###########################################################################
            ##  WZ control plots #####################################################
            ###########################################################################
            plotter.set_subdir('mc_data/wz_enhanced')
        
            plotter.plot_mc_vs_data('ss/tau_os/p1p2p3_enhance_wz', 'e2_t_Mass#LT', xaxis='m_{#mu#mu} (GeV)', 
                                    xrange=(0, 120), rebin=10, leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-ss-p1p2p3-enhance_wz-subMass')
        
            plotter.plot_mc_vs_data('ss/tau_os/p1p2p3_enhance_wz', 'e1_t_Mass', xaxis='m_{#mu#mu} (GeV)', 
                                    xrange=(0, 120), rebin=10, leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-ss-p1p2p3-enhance_wz-leadMass')
                    
            plotter.set_subdir('WZ_enhanced')

            plotter.plot_final_wz('e1_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-leadMass')

            plotter.plot_final_wz('e2_t_Mass#LT', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-subMass')

            plotter.plot_final_wz('e2Pt#LT', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-e2Pt')
        
            plotter.plot_final_wz('e2JetPt#LT', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-e2JetPt')

        if not options.no_signal:
            ###########################################################################
            ##  Signal region plots    ################################################
            ###########################################################################
            plotter.set_subdir('')
                    
            for label, proj_range in categories.iteritems():
                factor = 1.5 if label == 'LTHigh' else 1
                for tau_charge in ['tau_os', 'tau_ss']:
                    if tau_charge == 'tau_os':
                        plotter.set_subdir('%s' % label)
                    else:
                        plotter.set_subdir('%s_charge3' % label)

                    plotter.plot_final('e2_t_Mass#LT', rebin_slim, xaxis='m_{e_{2}#tau} (GeV)', 
                                       maxy=None, project=proj_range, project_axis='X', 
                                       differential=True, yaxis='Events / GeV', show_error=True, 
                                       tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-subMass-%s' % label)

                    plotter.plot_final('e2_t_Mass#LT', rebin_slim, xaxis='m_{e_{2}#tau} (GeV)', 
                                       maxy=None, project=proj_range, project_axis='X', 
                                       differential=True, yaxis='Events / GeV', show_error=True,
                                       tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.canvas.SetLogy(True)
                    plotter.save('final-subMass-logscale-%s' % label)

                    plotter.plot_final('e2_t_Mass#LT', 300, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, 
                                       project=proj_range, project_axis='X', tau_charge=tau_charge, 
                                       show_error=True)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-subMass-%s-counting' % label, dotc=True, dotroot=True)

                    #pt
                    plotter.plot_final("e1Pt#LT"    , int(factor*10), xaxis='p_{Te_{1}} (GeV)', 
                                       maxy=None, project=proj_range, project_axis='X', 
                                       tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-e1Pt-%s' % label)

                    plotter.plot_final("e2Pt#LT", int(factor*10), xaxis='p_{Te_{2}} (GeV)', maxy=None, 
                                       project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-e2Pt-%s' % label)

                    plotter.plot_final("tPt#LT", int(factor*10), xaxis='p_{T#tau} (GeV)', maxy=None, 
                                       project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-tPt-%s' % label)

                    #plotter.plot_final("e1JetPt#LT", int(factor*10), xaxis='p_{T Jet e_{1}} (GeV)', 
                    #                   maxy=None, project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e1JetPt-%s' % label)
                    #
                    #plotter.plot_final("e2JetPt#LT", int(factor*10), xaxis='p_{T Jet e_{2}} (GeV)', 
                    #                   maxy=None, project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e2JetPt-%s' % label)

                    #eta
                    plotter.plot_final("e1AbsEta#LT", 10, xaxis='|#eta_{e_{1}}|', maxy=None, project=proj_range, 
                                       project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-e1AbsEta-%s' % label)

                    plotter.plot_final("e2AbsEta#LT", 10, xaxis='|#eta_{e_{2}}|', maxy=None, 
                                       project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-e2AbsEta-%s' % label)

                    #DR
                    #plotter.plot_final("e1_t_DR#LT", 10, xaxis='#DeltaR_{e_{1}#tau}', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e1_t_DR-%s' % label)
                    #
                    #plotter.plot_final("e2_t_DR#LT", 10, xaxis='#DeltaR_{e_{2}#tau}', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e2_t_DR-%s' % label)

            #plotter.plot_final('e2RelPFIsoDB', 10, maxy='auto')
            #plotter.add_cms_blurb(sqrts)
            #plotter.save('final-e2Iso')

            ## plotter.plot_final('metSignificance', 5)
            ## plotter.add_cms_blurb(sqrts)
            ## plotter.save('final-metSig')

            plotter.plot_final('LT', 10, maxy='auto', xaxis='m_{e_{1}#tau} (GeV)')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-LT')

        if not options.no_f3:
            ###########################################################################
            ##  f3 region plots        ################################################
            ###########################################################################
                    
            for label, proj_range in categories.iteritems():
                factor = 1.5 if label == 'LTHigh' else 1
                for tau_charge in ['tau_os', 'tau_ss']:
                    if tau_charge == 'tau_os':
                        plotter.set_subdir('f3/%s' % label)
                    else:
                        plotter.set_subdir('f3/%s_charge3' % label)

                    plotter.plot_final_f3('e2_t_Mass#LT', rebin_slim, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, 
                                          project=proj_range, project_axis='X', differential=True, 
                                          yaxis='Events / GeV', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-subMass-%s' % label)

                    plotter.plot_final_f3('e2_t_Mass#LT', 300, xaxis='m_{e_{2}#tau} (GeV)', maxy=None, 
                                          project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-subMass-%s-counting' % label, dotc=True, dotroot=True)

                    #pt
                    plotter.plot_final_f3("e1Pt#LT", int(factor*10), xaxis='p_{Te_{1}} (GeV)', maxy=None, 
                                          project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-e1Pt-%s' % label)

                    plotter.plot_final_f3("e2Pt#LT", int(factor*10), xaxis='p_{Te_{2}} (GeV)', maxy=None, 
                                          project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-e2Pt-%s' % label)

                    plotter.plot_final_f3("tPt#LT", int(factor*10), xaxis='p_{T#tau} (GeV)', maxy=None, 
                                          project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-tPt-%s' % label)

                    #plotter.plot_final_f3("e1JetPt#LT", int(factor*10), xaxis='p_{T Jet e_{1}} (GeV)', 
                    #                      maxy=None, project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-e1JetPt-%s' % label)
                    #
                    #plotter.plot_final_f3("e2JetPt#LT" , int(factor*10), xaxis='p_{T Jet e_{2}} (GeV)', 
                    #                      maxy=None, project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-e2JetPt-%s' % label)

                    #eta
                    plotter.plot_final_f3("e1AbsEta#LT", 10, xaxis='|#eta_{e_{1}}|', maxy=None, 
                                          project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-e1AbsEta-%s' % label)

                    plotter.plot_final_f3("e2AbsEta#LT", 10, xaxis='|#eta_{e_{2}}|', maxy=None, 
                                          project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-e2AbsEta-%s' % label)

                    #DR
                    #plotter.plot_final_f3("e1_t_DR#LT", 10, xaxis='#DeltaR_{e_{1}#tau}', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-e1_t_DR-%s' % label)
                    #
                    #plotter.plot_final_f3("e2_t_DR#LT", 10, xaxis='#DeltaR_{e_{2}#tau}', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-e2_t_DR-%s' % label)



            plotter.set_subdir('f3')

            plotter.plot_final_f3('LT', 5, show_error=True)
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-LT')

            ###########################################################################
            ##  charge flip region plots        #######################################
            ###########################################################################
            #plotter.set_subdir('charge_flip_CR_f3')
            #
            #plotter.plot_final_f3('charge_flip_CR/e2_t_Mass', 20, xaxis='m_{e_{2}#tau} (GeV)', show_error=True)
            #plotter.add_cms_blurb(sqrts)
            #plotter.save('final-chargeFlip-subMass')
            #
            #plotter.plot_final_f3('charge_flip_CR/e1_t_Mass', 20, xaxis='m_{e_{1}#tau} (GeV)', show_error=True)
            #plotter.add_cms_blurb(sqrts)
            #plotter.save('final-chargeFlip-e1tMass')
            #
            #plotter.plot_final_f3('charge_flip_CR/e1_e2_Mass', 5, xaxis='m_{ee} (GeV)', show_error=True)
            #plotter.add_cms_blurb(sqrts)
            #plotter.save('final-chargeFlip-e1e2Mass')
        
        #END OF if not options.dry_run:
    if not options.no_shapes:
        ###########################################################################
        ##  Making shape file     #################################################
        ###########################################################################
        plotter.set_subdir('')
        prefixes = [options.prefix+'$'] if options.prefix else ['']
        prefixes = [i+'$' for i in options.prefixes.split(',') if i] if options.prefixes else prefixes
        for prefix in prefixes:
            shape_prefix = prefix if len(prefixes) > 1 else ''
            shape_prefix = shape_prefix.replace(':','_').replace('$','_')
            binning_7TeV = [0,40,80,120,200]

            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, 'LTCut_eet_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
            shape_dir = shape_file.mkdir('eetCatHigh')
            plotter.write_shapes(prefix+'e2_t_Mass#LT', binning_7TeV, shape_dir, qcd_fraction=0., project=[80, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHigh_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', 20, shape_dir, qcd_fraction=0.0, project=[80, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHigh_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', 20, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')

            logging.warning('shape file %s created' % shape_file.GetName())
            shape_file.Close()


            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, 'LTCut_eet_shapes_%s_different_fakes.root' % ( plotter.period) ), 'RECREATE')
            shape_dir = shape_file.mkdir('eetCatHigh')
            plotter.write_shapes(prefix+'e2_t_Mass#LT', binning_7TeV, shape_dir, qcd_fraction=0., project=[80, 650], project_axis='X', different_fakes=True)
            #shape_dir = shape_file.mkdir('eetCatHigh_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', 20, shape_dir, qcd_fraction=0.0, project=[80, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHigh_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', 20, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')

            logging.warning('shape file %s created' % shape_file.GetName())
            shape_file.Close()


            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, 'LTAll_eet_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
            shape_dir = shape_file.mkdir('eet')
            plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[0, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eet_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eet_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0, 650], project_axis='X')

            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()


            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, '%seet_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            
            shape_dir = shape_file.mkdir('eetCatLow')
            plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[0, 100], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatLow_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0, 100], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatLow_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0, 100], project_axis='X')
            
            shape_dir = shape_file.mkdir('eetCatHigh')
            plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[100, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHigh_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[100, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHigh_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[100, 650], project_axis='X')
            
            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()



            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, '%seet_shapes_%s_different_fakes.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            
            shape_dir = shape_file.mkdir('eetCatLow')
            plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[0, 100], project_axis='X', different_fakes=True)
            #shape_dir = shape_file.mkdir('eetCatLow_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0, 100], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatLow_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0, 100], project_axis='X')
            
            shape_dir = shape_file.mkdir('eetCatHigh')
            plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[100, 650], project_axis='X', different_fakes=True)
            #shape_dir = shape_file.mkdir('eetCatHigh_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[100, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHigh_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[100, 650], project_axis='X')
            
            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()


            #shape_file = ROOT.TFile(
            #    os.path.join(plotter.outputdir, '%seet_f3_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            #
            #shape_dir = shape_file.mkdir('eetCatLow')
            #plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
            ##shape_dir = shape_file.mkdir('eetCatLow_w')
            ##plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            ##shape_dir = shape_file.mkdir('eetCatLow_q')
            ##plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('eetCatHigh')
            #plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
            ##shape_dir = shape_file.mkdir('eetCatHigh_w')
            ##plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            ##shape_dir = shape_file.mkdir('eetCatHigh_q')
            ##plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #logging.warning('shape file %s created' % shape_file.GetName()) 
            #shape_file.Close()
            #
            #shape_file = ROOT.TFile(
            #    os.path.join(plotter.outputdir, '%seet_f3_all_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            #
            #shape_dir = shape_file.mkdir('eetCatLow')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatLow_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatLow_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('eetCatHigh')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHigh_w')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHigh_q')
            #plotter.write_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('eetCatLowf3')
            #plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatLowf3_w')
            #plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatLowf3_q')
            #plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('eetCatHighf3')
            #plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHighf3_w')
            #plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('eetCatHighf3_q')
            #plotter.write_f3_shapes(prefix+'e2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #logging.warning('shape file %s created' % shape_file.GetName()) 
            #shape_file.Close()
