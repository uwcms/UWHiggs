'''

Plot the EMT channel

Usage: python WHPlotterEMT.py

'''

import glob
import logging
import os
import ROOT
import sys
import WHPlotterBase
from WHPlotterBase import make_styler, parser, project_x
import rootpy.plotting.views as views
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors

logging.basicConfig(stream=sys.stderr, level=logging.WARNING)

class WHPlotterEMT(WHPlotterBase.WHPlotterBase):
    def __init__(self):
        obj2_charge_mapper={
            "subMass#LT"          : "subMass*#LT"         ,
            "subMass#tPt"         : "subMass*#tPt"        ,
            "e_t_Mass"            : "e*_t_Mass"           ,
            "e_m_Mass"            : "e*_m_Mass"           ,
            "subMass"             : "subMass*"            ,} 
        super(WHPlotterEMT, self).__init__('EMT', {}, obj2_charge_mapper)

categories = {
    'LTCut' : [80, 650],
    'LTLow' : [0, 130],
    'LTHigh': [130, 650],
    'Full'  : [0, 650],
}
rebin_slim = [20]+range(30, 81, 10)+[100,130,200]
binning_7TeV = [0,40,80,120,200]

if __name__ == "__main__":
    plotter = WHPlotterEMT()
    sqrts   = plotter.sqrts
    plotter.defaults['show_charge_fakes'] = True
    options,NOTUSED = parser.parse_args()
    if not options.dry_run:
        if not options.no_mc_data:
            ###########################################################################
            ##  Ztt control plots #####################################################
            ###########################################################################
            plotter.set_subdir('mc_data/em')

            # Control Z->tautau + jet region
            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'e_m_Mass#LT', rebin=10, xaxis='m_{e#mu} (GeV)', 
                                    leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-em-emMass')


            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'mPt#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-em-mPt')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'ePt#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-em-ePt')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'mAbsEta#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-em-mAbsEta')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'eAbsEta#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-em-eAbsEta')

        if not options.no_wz:
            ###########################################################################
            ##  WZ control plots #####################################################
            ###########################################################################
            plotter.set_subdir('mc_data/wz_enhanced')
        
            plotter.plot_mc_vs_data('ss/tau_os/p1p2p3_enhance_wz', 'e_t_Mass#LT', xaxis='m_{e#tau} (GeV)', 
                                    xrange=(0, 120), rebin=10, leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-ss-p1p2p3-enhance_wz-e_t_Mass')
        
            plotter.plot_mc_vs_data('ss/tau_os/p1p2p3_enhance_wz', 'm_t_Mass#LT', xaxis='m_{#mu#tau} (GeV)', 
                                    xrange=(0, 120), rebin=10, leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-ss-p1p2p3-enhance_wz-leadMass')
                    
            plotter.set_subdir('WZ_enhanced')

            plotter.plot_final_wz('e_t_Mass#LT', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-e_t_Mass')

            plotter.plot_final_wz('m_t_Mass#LT', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-m_t_Mass')

            plotter.plot_final_wz('ePt#LT', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-ePt')
        
            plotter.plot_final_wz('eJetPt#LT', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-eJetPt')

        if not options.no_signal:
            ###########################################################################
            ##  Signal region plots    ################################################
            ###########################################################################
            plotter.set_subdir('')

            for label, proj_range in categories.iteritems():
                for tau_charge in ['tau_os', 'tau_ss']:
                    if tau_charge == 'tau_os':
                        plotter.set_subdir('%s' % label)
                    else:
                        plotter.set_subdir('%s_charge3' % label)

                    plotter.plot_final('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, 
                                       project=proj_range, project_axis='X', differential=True, 
                                       yaxis='Events / GeV', show_error=True, tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-subMass-%s' % label)

                    plotter.plot_final('subMass#LT', binning_7TeV, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, 
                                       project=proj_range, project_axis='X', differential=True, 
                                       yaxis='Events / GeV', show_error=True, tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-subMass-7TeVBin-%s' % label)


                    #plotter.plot_final("e_m_Mass#LT", 20, xaxis='m_{e#mu} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e_m_Mass-%s' % label)
                    #
                    #plotter.plot_final("m_t_Mass#LT", 20, xaxis='m_{#mu#tau} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m_t_Mass-%s' % label)
                    #
                    #plotter.plot_final("e_t_Mass#LT", 20, xaxis='m_{e#tau} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e_t_Mass-%s' % label)

                    plotter.plot_final('subMass#LT', 300, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, 
                                       project=proj_range, project_axis='X', tau_charge=tau_charge, 
                                       show_error=True)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-subMass-%s-counting' % label, dotc=True, dotroot=True)

                    #pt
                    #plotter.plot_final("mPt#LT"    , 10, xaxis='p_{T#mu} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-mPt-%s' % label)
                    #
                    #plotter.plot_final("ePt#LT"    , 10, xaxis='p_{Te} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-ePt-%s' % label)
                    #
                    #plotter.plot_final("tPt#LT"    , 10, xaxis='p_{T#tau} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-tPt-%s' % label)

                    #plotter.plot_final("mJetPt#LT" , 10, xaxis='p_{T Jet#mu} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-mJetPt-%s' % label)
                    #
                    #plotter.plot_final("eJetPt#LT" , 10, xaxis='p_{T Jet e} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-eJetPt-%s' % label)

                    #eta
                    #plotter.plot_final("mAbsEta#LT", 10, xaxis='|#eta_{#mu}|', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-mAbsEta-%s' % label)
                    #
                    #plotter.plot_final("eAbsEta#LT", 10, xaxis='|#eta_{e}|'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-eAbsEta-%s' % label)
                    #
                    #plotter.plot_final("tAbsEta#LT", 10, xaxis='|#eta_{#tau}|'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-tAbsEta-%s' % label)

                    #DR
                    #plotter.plot_final("m_t_DR#LT", 10, xaxis='#DeltaR_{#mu#tau}', maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m_t_DR-%s' % label)
                    #
                    #plotter.plot_final("e_t_DR#LT", 10, xaxis='#DeltaR_{e#tau}'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e_t_DR-%s' % label)

                    #Jet BTag
                    #plotter.plot_final("eJetBtag#LT", 2, xaxis='e Jet Btag'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-eJetBtag-%s' % label)
                    #
                    #plotter.plot_final("mJetBtag#LT", 2, xaxis='#mu Jet Btag'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-mJetBtag-%s' % label)


            plotter.plot_final('LT', 5, xaxis='LT (GeV)', maxy=15)
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-LT')
            
        if not options.no_f3:
            ###########################################################################
            ##  F3 enhanced region plots    ###########################################
            ###########################################################################
            plotter.set_subdir('f3')
                    
            for label, proj_range in categories.iteritems():
                for tau_charge in ['tau_os', 'tau_ss']:
                    if tau_charge == 'tau_os':
                        plotter.set_subdir('f3/%s' % label)
                    else:
                        plotter.set_subdir('f3/%s_charge3' % label)

                    plotter.plot_final_f3('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, 
                                          project=proj_range, project_axis='X', differential=True, 
                                          yaxis='Events / GeV', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-subMass-%s' % label)

                    plotter.plot_final_f3('subMass#LT', binning_7TeV, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, 
                                          project=proj_range, project_axis='X', differential=True, 
                                          yaxis='Events / GeV', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-subMass-7TeVBin-%s' % label)

                    #plotter.plot_final_f3('subMass#LT', [20, 40, 120, 300], xaxis='m_{l_{2}#tau} (GeV)', maxy=None, 
                    #                      project=proj_range, project_axis='X', differential=True, 
                    #                      yaxis='Events / GeV', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-subMass-widebin-%s' % label)
                    #
                    #plotter.plot_final_f3('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-subMass-%s-notDifferential' % label)

                    #plotter.plot_final_f3("e_m_Mass#LT", 20, xaxis='m_{e#mu} (GeV)', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e_m_Mass-%s' % label)
                    #
                    #plotter.plot_final_f3("m_t_Mass#LT", 20, xaxis='m_{#mu#tau} (GeV)', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m_t_Mass-%s' % label)
                    #
                    #plotter.plot_final_f3("e_t_Mass#LT", 20, xaxis='m_{e#tau} (GeV)'  , maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-e_t_Mass-%s' % label)

                    plotter.plot_final_f3('subMass#LT', 300, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, 
                                          project=proj_range, project_axis='X', tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-subMass-%s-counting' % label, dotc=True, dotroot=True)

                    #pt
                    #plotter.plot_final_f3("mPt#LT"    , 10, xaxis='p_{T#mu} (GeV)', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-mPt-%s' % label)
                    #
                    #plotter.plot_final_f3("ePt#LT"    , 10, xaxis='p_{Te} (GeV)', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-ePt-%s' % label)
                    #
                    #plotter.plot_final_f3("tPt#LT"    , 10, xaxis='p_{T#tau} (GeV)', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-tPt-%s' % label)

                    #plotter.plot_final_f3("mJetPt#LT" , 10, xaxis='p_{T Jet#mu} (GeV)', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-mJetPt-%s' % label)
                    #
                    #plotter.plot_final_f3("eJetPt#LT" , 10, xaxis='p_{T Jet e} (GeV)', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-eJetPt-%s' % label)

                    #eta
                    #plotter.plot_final_f3("mAbsEta#LT", 10, xaxis='|#eta_{#mu}|', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-mAbsEta-%s' % label)
                    #
                    #plotter.plot_final_f3("eAbsEta#LT", 10, xaxis='|#eta_{e}|'  , maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-eAbsEta-%s' % label)
                    #
                    #plotter.plot_final_f3("tAbsEta#LT", 10, xaxis='|#eta_{#tau}|'  , maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-tAbsEta-%s' % label)

                    #DR
                    #plotter.plot_final_f3("m_t_DR#LT", 10, xaxis='#DeltaR_{#mu#tau}', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m_t_DR-%s' % label)
                    #
                    #plotter.plot_final_f3("e_t_DR#LT", 10, xaxis='#DeltaR_{e#tau}'  , maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-e_t_DR-%s' % label)

                    #Jet BTag
                    #plotter.plot_final_f3("eJetBtag#LT", 2, xaxis='e Jet Btag'  , maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-eJetBtag-%s' % label)
                    #
                    #plotter.plot_final_f3("mJetBtag#LT", 2, xaxis='#mu Jet Btag'  , maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-mJetBtag-%s' % label)

            plotter.set_subdir('f3')
            plotter.plot_final_f3('LT', 5, xaxis='LT (GeV)', show_error=True)
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-f3-LT')

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
                os.path.join(plotter.outputdir, 'LTAll_emt_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
            shape_dir = shape_file.mkdir('emt')
            plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[0, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emt_w')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emt_q')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0, 650], project_axis='X')

            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()


            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, 'LTCut_emt_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
            shape_dir = shape_file.mkdir('emt')
            plotter.write_shapes(prefix+'subMass#LT', binning_7TeV, shape_dir, qcd_fraction=0., project=[80, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_w')
            #plotter.write_shapes(prefix+'subMass#LT', 20, shape_dir, qcd_fraction=0.0, project=[80, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_q')
            #plotter.write_shapes(prefix+'subMass#LT', 20, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')

            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()


            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, 'LTCut_emt_shapes_%s_different_fakes.root' % ( plotter.period) ), 'RECREATE')
            shape_dir = shape_file.mkdir('emt')
            plotter.write_shapes(prefix+'subMass#LT', binning_7TeV, shape_dir, qcd_fraction=0., project=[80, 650], project_axis='X', different_fakes=True)
            #shape_dir = shape_file.mkdir('emtCatHigh_w')
            #plotter.write_shapes(prefix+'subMass#LT', 20, shape_dir, qcd_fraction=0.0, project=[80, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_q')
            #plotter.write_shapes(prefix+'subMass#LT', 20, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')

            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()

            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, '%semt_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            
            shape_dir = shape_file.mkdir('emtCatLow')
            plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[0, 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLow_w')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0, 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLow_q')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0, 130], project_axis='X')
            
            shape_dir = shape_file.mkdir('emtCatHigh')
            plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_w')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_q')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            
            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()



            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, '%semt_shapes_%s_different_fakes.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            
            shape_dir = shape_file.mkdir('emtCatLow')
            plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[0, 130], project_axis='X', different_fakes=True)
            #shape_dir = shape_file.mkdir('emtCatLow_w')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0, 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLow_q')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0, 130], project_axis='X')
            
            shape_dir = shape_file.mkdir('emtCatHigh')
            plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[130, 650], project_axis='X', different_fakes=True)
            #shape_dir = shape_file.mkdir('emtCatHigh_w')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_q')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            
            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()


            #shape_file = ROOT.TFile(
            #    os.path.join(plotter.outputdir, '%semt_f3_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            #
            #shape_dir = shape_file.mkdir('emtCatLow')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLow_w')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLow_q')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('emtCatHigh')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_w')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_q')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #logging.warning('shape file %s created' % shape_file.GetName()) 
            #shape_file.Close()
            #
            #shape_file = ROOT.TFile(
            #    os.path.join(plotter.outputdir, '%semt_f3_all_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            #
            #shape_dir = shape_file.mkdir('emtCatLow')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLow_w')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLow_q')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('emtCatHigh')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_w')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHigh_q')
            #plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('emtCatLowf3')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLowf3_w')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatLowf3_q')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('emtCatHighf3')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHighf3_w')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('emtCatHighf3_q')
            #plotter.write_f3_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #logging.warning('shape file %s created' % shape_file.GetName()) 
            #shape_file.Close()
