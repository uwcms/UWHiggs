'''

Plot the MMT channel

Usage: python WHPlotterMMT.py

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

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
ROOT.gStyle.SetOptTitle(0)
#ROOT.gROOT.SetStyle('Plain')

class WHPlotterMMT(WHPlotterBase.WHPlotterBase):
    def __init__(self):
        super(WHPlotterMMT, self).__init__('MMT')

rebin_slim = range(20, 81, 10)+[100, 130, 200]
binning_7TeV = [0,40,80,120,200]
categories = {
    'LTCut' : [80, 650],
    'Full'  : [0,  650],
    'LTLow' : [0,  130],
    'LTHigh': [130, 650],
}



if __name__ == "__main__":
    plotter = WHPlotterMMT()
    sqrts   = plotter.sqrts
    options,NOTUSED = parser.parse_args()
    if not options.dry_run:
        if not options.no_mc_data:
            ###########################################################################
            ##  Zmm control plots #####################################################
            ###########################################################################
            plotter.set_subdir('mc_data/zmm')
        
            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm1_m2_Mass', rebin=10, 
                                    xaxis='m_{#mu#mu} (GeV)', leftside=False)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zmm-m1m2Mass')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm1Pt#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zmm-m1Pt')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm2Pt#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zmm-m2Pt')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm1AbsEta#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zmm-m1AbsEta')

            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'm2AbsEta#LT', rebin=10, 
                                    xaxis='p_{T}^{#mu1} (GeV)', leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-Zmm-m2AbsEta')

            # Check PU variables
            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'rho')
            plotter.save('mcdata-Zmm-rho')
        
            plotter.plot_mc_vs_data('os/tau_?s/p1p2?3', 'nvtx')
            plotter.save('mcdata-Zmm-nvtx')

        if not options.no_wz:
            ###########################################################################
            ##  WZ control plots #####################################################
            ###########################################################################
            plotter.set_subdir('mc_data/wz_enhanced')
        
            plotter.plot_mc_vs_data('ss/tau_os/p1p2p3_enhance_wz', 'm2_t_Mass#LT', xaxis='m_{#mu#mu} (GeV)', 
                                    xrange=(0, 120), rebin=10, leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-ss-p1p2p3-enhance_wz-subMass')
        
            plotter.plot_mc_vs_data('ss/tau_os/p1p2p3_enhance_wz', 'm1_t_Mass', xaxis='m_{#mu#mu} (GeV)', 
                                    xrange=(0, 120), rebin=10, leftside=False, preprocess=project_x)
            plotter.add_cms_blurb(sqrts)
            plotter.save('mcdata-ss-p1p2p3-enhance_wz-leadMass')
                    
            plotter.set_subdir('WZ_enhanced')

            plotter.plot_final_wz('m1_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-leadMass')

            plotter.plot_final_wz('m2_t_Mass#LT', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-subMass')

            plotter.plot_final_wz('m2Pt#LT', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-m2Pt')
        
            plotter.plot_final_wz('m2JetPt#LT', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)',
                                  project=[0,650], project_axis='X')
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-wz-m2JetPt')

        if not options.no_signal:
            ###########################################################################
            ##  Signal region plots    ################################################
            ###########################################################################
            plotter.set_subdir('')
            #rebin_slim = range(20, 91, 10)+[110,200]
                    
            for label, proj_range in categories.iteritems():
                for tau_charge in ['tau_os', 'tau_ss']:
                    if tau_charge == 'tau_os':
                        plotter.set_subdir('%s' % label)
                    else:
                        plotter.set_subdir('%s_charge3' % label)

                    plotter.plot_final('m2_t_Mass#LT', rebin_slim, xaxis='m_{#mu_{2}#tau} (GeV)', maxy=None, 
                                       project=proj_range, project_axis='X', differential=True, 
                                       yaxis='Events / GeV', show_error=True, tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-subMass-%s' % label)

                    plotter.plot_final('m2_t_Mass#LT', binning_7TeV, xaxis='m_{#mu_{2}#tau} (GeV)', 
                                       project=proj_range, project_axis='X', 
                                       yaxis='Events', show_error=True, x_range=[0,199], 
                                       maxy=15.6, tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-subMass-7TeVBin-%s' % label, verbose=True)

                    plotter.plot_final('m2_t_Mass#LT', 300, xaxis='m_{#mu_{2}#tau} (GeV)', maxy=None, 
                                       project=proj_range, project_axis='X', show_error=True, tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-subMass-%s-counting' % label, json=True) #, dotc=True, dotroot=True)
                    
                    #pt
                    #plotter.plot_final("m1Pt#LT"    , 10, xaxis='p_{T#mu_{1}} (GeV)', 
                    #                   maxy=None, project=proj_range, project_axis='X', 
                    #                   show_error=True, tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m1Pt-%s' % label)
                    #
                    #plotter.plot_final("m2Pt#LT"    , 10, xaxis='p_{T#mu_{2}} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', 
                    #                   show_error=True, tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m2Pt-%s' % label)
                    #
                    #plotter.plot_final("tPt#LT"    , 10, xaxis='p_{T#tau} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', 
                    #                   show_error=True, tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-tPt-%s' % label)
                    
                    #plotter.plot_final("m1JetPt#LT" , 10, xaxis='p_{T Jet #mu_{1}} (GeV)', 
                    #                   maxy=None, project=proj_range, project_axis='X', 
                    #                   show_error=True, tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m1JetPt-%s' % label)
                    #
                    #plotter.plot_final("m2JetPt#LT" , 10, xaxis='p_{T Jet #mu_{2}} (GeV)', 
                    #                   maxy=None, project=proj_range, project_axis='X', 
                    #                   show_error=True, tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m2JetPt-%s' % label)
                    #
                    #plotter.plot_final("tPt#LT"    , 10, xaxis='p_{T#tau} (GeV)', maxy=None, 
                    #                   project=proj_range, project_axis='X', show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-tPt-%s' % label)
                    
                    #eta
                    #plotter.plot_final("m1AbsEta#LT", 10, xaxis='|#eta_{#mu_{1}}|', maxy=None, 
                    #                   project=proj_range, project_axis='X', show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m1AbsEta-%s' % label)
                    #
                    #plotter.plot_final("m2AbsEta#LT", 10, xaxis='|#eta_{#mu_{2}}|'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m2AbsEta-%s' % label)
                    #
                    #plotter.plot_final("tAbsEta#LT", 10, xaxis='|#eta_{#tau}|'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-tAbsEta-%s' % label)
                    
                    #DR
                    #plotter.plot_final("m1_t_DR#LT", 10, xaxis='#DeltaR_{#mu_{1}#tau}', maxy=None, 
                    #                   project=proj_range, project_axis='X', show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m1_t_DR-%s' % label)
                    #
                    #plotter.plot_final("m2_t_DR#LT", 10, xaxis='#DeltaR_{#mu_{2}#tau}'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m2_t_DR-%s' % label)
                    
                    #Jet BTag
                    
                    #from pdb import set_trace; set_trace()
                    #plotter.plot_final("m2JetBtag#LT", [-100, -40, 20, 100], xaxis='#mu_{2} Jet Btag'  , 
                    #                   maxy=None, project=proj_range, project_axis='X', show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m2JetBtag-%s' % label, dotroot=True)
                    #
                    #plotter.plot_final("m1JetBtag#LT", 1, xaxis='#mu_{1} Jet Btag'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m1JetBtag-%s' % label)
                    #
                    #plotter.plot_final("m2JetCSVBtag#LT", 1, xaxis='#mu_{2} Jet CSV Btag'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', x_range=[0,1], show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m2JetCSVBtag-%s' % label, dotroot=True)
                    #
                    #plotter.plot_final("m1JetCSVBtag#LT", 1, xaxis='#mu_{1} Jet CSV Btag'  , maxy=None, 
                    #                   project=proj_range, project_axis='X', x_range=[0,1], show_error=True, 
                    #                   tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-m1JetCSVBtag-%s' % label)

            plotter.set_subdir('')
            plotter.plot_final('LT', 10, xaxis='LT (GeV)', maxy=15)
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-LT')

        if not options.no_f3:
            ###########################################################################
            ##  F3 enhanced region plots    ###########################################
            ###########################################################################
            plotter.set_subdir('f3')
            #plotter.plot_final_f3('LT', [0, 50, 75, 100, 125, 150, 200, 300, 500], xaxis='LT (GeV)', maxy=None, 
            #                      show_ratio=True, fit = { 'model' : ("slope*x + offset", "offset[0, -1, 1], slope[0,-1,1]"), 
            #                                               'range' : [50, 500], 'options' : 'QIRMENS', 
            #                                               'stat position' : 'bottom-left'})
            #plotter.add_cms_blurb(sqrts)
            #plotter.save('final-LT-tauOS-chi2', verbose=True)

            plotter.plot_final_f3('LT', [0, 50, 75, 100, 125, 150, 200, 300, 500], xaxis='LT (GeV)', maxy=None, 
                                  show_ratio=True, fit = { 'model' : ("slope*x + offset", "offset[0, -1, 1], slope[0,-1,1]"), 
                                                           'range' : [50, 500], 'options' : 'WLQIRMENS', 
                                                           'stat position' : 'bottom-left'})
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-LT-tauOS-likelihood')

            plotter.plot_final_f3('LT', [0, 50, 130, 500], xaxis='LT (GeV)', maxy=None, show_ratio=True)
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-LT-tauOS-2categories')

            #plotter.plot_final_f3('LT', [0, 50, 75, 100, 125, 150, 200, 300, 500], xaxis='LT (GeV)', tau_charge='tau_ss', 
            #                      maxy=None, show_ratio=True, 
            #                      fit = {'range' : [50, 500], 'model' : ("slope*x + offset", "offset[0, -1, 1], slope[0,-1,1]")})
            #plotter.add_cms_blurb(sqrts)
            #plotter.save('final-LT-tauSS')

            plotter.plot_final_f3('LT', [0, 50, 130, 500], xaxis='LT (GeV)', tau_charge='tau_ss', maxy=None, show_ratio=True)
            plotter.add_cms_blurb(sqrts)
            plotter.save('final-LT-tauSS-2categories')

            for label, proj_range in categories.iteritems():
                for tau_charge in ['tau_os', 'tau_ss']:
                    if tau_charge == 'tau_os':
                        plotter.set_subdir('f3/%s' % label)
                    else:
                        plotter.set_subdir('f3/%s_charge3' % label)

                    plotter.plot_final_f3('m2_t_Mass#LT', rebin_slim, xaxis='m_{#mu_{2}#tau} (GeV)', 
                                          maxy=None, project=proj_range, project_axis='X', differential=True,
                                          show_error=True, tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-subMass-%s' % label)

                    plotter.plot_final_f3('m2_t_Mass#LT', binning_7TeV, xaxis='m_{#mu_{2}#tau} (GeV)', 
                                          maxy=None, project=proj_range, project_axis='X', differential=True,
                                          show_error=True, tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-subMass-7TeVBin-%s' % label)

                    plotter.plot_final_f3('m2_t_Mass#LT', 300, xaxis='m_{#mu_{2}#tau} (GeV)', 
                                          maxy=None, project=proj_range, project_axis='X',
                                          tau_charge=tau_charge)
                    plotter.add_cms_blurb(sqrts)
                    plotter.save('final-f3-subMass-%s-counting' % label, json=True, verbose=True)

                    #pt
                    #plotter.plot_final_f3("m1Pt#LT"    , 10, xaxis='p_{T#mu_{1}} (GeV)', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m1Pt-%s' % label)
                    #
                    #plotter.plot_final_f3("m2Pt#LT"    , 10, xaxis='p_{T#mu_{2}} (GeV)', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m2Pt-%s' % label)
                    #
                    #plotter.plot_final_f3("tPt#LT", 10, xaxis='p_{T#tau} (GeV)', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-tPt-%s' % label)
                    #
                    ##plotter.plot_final_f3("m1JetPt#LT" , 10, xaxis='p_{T Jet #mu_{1}} (GeV)', 
                    ##                      maxy=None, project=proj_range, project_axis='X', 
                    ##                      tau_charge=tau_charge)
                    ##plotter.add_cms_blurb(sqrts)
                    ##plotter.save('final-f3-m1JetPt-%s' % label)
                    ##
                    ##plotter.plot_final_f3("m2JetPt#LT" , 10, xaxis='p_{T Jet #mu_{2}} (GeV)', 
                    ##                      maxy=None, project=proj_range, project_axis='X', 
                    ##                      tau_charge=tau_charge)
                    ##plotter.add_cms_blurb(sqrts)
                    ##plotter.save('final-f3-m2JetPt-%s' % label)
                    #
                    ##eta
                    #plotter.plot_final_f3("m1AbsEta#LT", 10, xaxis='|#eta_{#mu_{1}}|', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m1AbsEta-%s' % label)
                    #
                    #plotter.plot_final_f3("m2AbsEta#LT", 10, xaxis='|#eta_{#mu_{2}}|', maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m2AbsEta-%s' % label)
                    #
                    #plotter.plot_final_f3("tAbsEta#LT", 10, xaxis='|#eta_{#tau}|'  , maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-tAbsEta-%s' % label)

                    #DR
                    #plotter.plot_final_f3("m1_t_DR#LT", 5, xaxis='#DeltaR_{#mu_{1}#tau}', 
                    #                      maxy=None, project=proj_range, project_axis='X',
                    #                      tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m1_t_DR-%s' % label)
                    #
                    #plotter.plot_final_f3("m2_t_DR#LT", 5, xaxis='#DeltaR_{#mu_{2}#tau}'  , maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m2_t_DR-%s' % label)

                    #Jet BTag
                
                    #from pdb import set_trace; set_trace()
                    #plotter.plot_final_f3("m2JetBtag#LT", [-100, -40, 20, 100], xaxis='#mu_{2} Jet Btag', 
                    #                      maxy=None, project=proj_range, project_axis='X', 
                    #                      tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m2JetBtag-%s' % label, dotroot=True)
                    #
                    #plotter.plot_final_f3("m1JetBtag#LT", 1, xaxis='#mu_{1} Jet Btag'  , maxy=None, 
                    #                      project=proj_range, project_axis='X', tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m1JetBtag-%s' % label)
                    #
                    #plotter.plot_final_f3("m2JetCSVBtag#LT", 1, xaxis='#mu_{2} Jet CSV Btag', maxy=None, 
                    #                      project=proj_range, project_axis='X', x_range=[0,1], 
                    #                      tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m2JetCSVBtag-%s' % label, dotroot=True)
                    #
                    #plotter.plot_final_f3("m1JetCSVBtag#LT", 1, xaxis='#mu_{1} Jet CSV Btag', maxy=None, 
                    #                      project=proj_range, project_axis='X', x_range=[0,1], 
                    #                      tau_charge=tau_charge)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-m1JetCSVBtag-%s' % label)

                    #TAU ISO
                    #plotter.plot_final_f3("tRawIso3Hits#LT", [0,5,10,15,20,25,30,35,40,50,60,120,200], 
                    #                      xaxis='#tau 3Hits Raw Iso', maxy=None, 
                    #                      project=proj_range, project_axis='X', 
                    #                      tau_charge=tau_charge, show_ratio=True,
                    #                      fit = { 'model' : ("slope*x + offset", "offset[0, -1, 1], slope[0,-1,1]"), 
                    #                              'options' : 'WLQIRMENS', 
                    #                              'stat position' : 'bottom-right'})
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-tRawIso3Hits-Likelihood-%s' % label)
                    #
                    #plotter.plot_final_f3("tRawIso3Hits#LT", [0,5,10,15,20,25,30,35,40,50,60,120,200], 
                    #                      xaxis='#tau 3Hits Raw Iso', maxy=None, 
                    #                      project=proj_range, project_axis='X', 
                    #                      tau_charge=tau_charge, show_ratio=True,
                    #                      fit = { 'model' : ("slope*x + offset", "offset[0, -1, 1], slope[0,-1,1]"), 
                    #                              'options' : 'QIRMENS', 
                    #                              'stat position' : 'bottom-right'})
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-tRawIso3Hits-chi2-%s' % label)
                    #
                    #plotter.plot_final_f3("tRawIso3Hits#LT", [0,5,10,15,20,25,30,35,40,50,60,120,200], 
                    #                      xaxis='#tau 3Hits Raw Iso', maxy=None, project=proj_range, 
                    #                      project_axis='X', tau_charge=tau_charge, show_ratio=True,
                    #                      wz_error=0., zz_error =0.)
                    #plotter.add_cms_blurb(sqrts)
                    #plotter.save('final-f3-tRawIso3Hits-noDibosonErr-%s' % label)

        #END OF if not options.dry_run:

    if not options.no_shapes:
        ###########################################################################
        ##  Making shape file     #################################################
        ###########################################################################
        plotter.set_subdir('')
        prefixes = [options.prefix+'$'] if options.prefix else ['']
        prefixes = [i+'$' for i in options.prefixes.split(',') if i] if options.prefixes else prefixes
        for prefix in prefixes:
            #plotter.plot_final(prefix+'m2_t_Mass', 20, xaxis='m_{#mu_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
            #plotter.add_cms_blurb(sqrts)
            #plotter.save('final-%s-qweight05-subMass' % options.prefix)
            #
            #plotter.plot_final_f3(prefix+'m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0.5, show_error=True)
            #plotter.add_cms_blurb(sqrts)
            #plotter.save('final-%s-f3-qweight05-subMass' % options.prefix)

            shape_prefix = prefix if len(prefixes) > 1 else ''
            shape_prefix = shape_prefix.replace(':','_').replace('$','_')
            binning_7TeV = [0,40,80,120,200]
            
            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, 'LTCut_mmt_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
            shape_dir = shape_file.mkdir('mmt')
            plotter.write_shapes(prefix+'m2_t_Mass#LT', binning_7TeV, shape_dir, qcd_fraction=0., project=[80, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatHigh_w')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', 20, shape_dir, qcd_fraction=-1., project=[80, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatHigh_q')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', 20, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')

            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()            
            
            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, 'LTAll_mmt_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
            shape_dir = shape_file.mkdir('mmt')
            plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[0, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmt_w')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[0, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmt_q')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0, 650], project_axis='X')
            
            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()
            
            #rebin_slim = range(20, 91, 10)+[110,200]
            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, '%smmt_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            
            shape_dir = shape_file.mkdir('mmtCatLow')
            plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatLow_w')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatLow_q')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            
            shape_dir = shape_file.mkdir('mmtCatHigh')
            plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0., project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatHigh_w')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatHigh_q')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            
            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()


            #shape_file = ROOT.TFile(
            #    os.path.join(plotter.outputdir, '%smmt_f3_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            #
            #shape_dir = shape_file.mkdir('mmtCatLow')
            #plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            ##shape_dir = shape_file.mkdir('mmtCatLow_w')
            ##plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[0., 130], project_axis='X')
            ##shape_dir = shape_file.mkdir('mmtCatLow_q')
            ##plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('mmtCatHigh')
            #plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            ##shape_dir = shape_file.mkdir('mmtCatHigh_w')
            ##plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[130, 650], project_axis='X')
            ##shape_dir = shape_file.mkdir('mmtCatHigh_q')
            ##plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #logging.warning('shape file %s created' % shape_file.GetName()) 
            #shape_file.Close()
            #
            #shape_file = ROOT.TFile(
            #    os.path.join(plotter.outputdir, '%smmt_f3_all_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            #
            #shape_dir = shape_file.mkdir('mmtCatLow')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatLow_w')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatLow_q')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('mmtCatHigh')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatHigh_w')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatHigh_q')
            #plotter.write_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('mmtCatLowf3')
            #plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatLowf3_w')
            #plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[0., 130], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatLowf3_q')
            #plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[0., 130], project_axis='X')
            #
            #shape_dir = shape_file.mkdir('mmtCatHighf3')
            #plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatHighf3_w')
            #plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=-1., project=[130, 650], project_axis='X')
            #shape_dir = shape_file.mkdir('mmtCatHighf3_q')
            #plotter.write_f3_shapes(prefix+'m2_t_Mass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[130, 650], project_axis='X')
            #
            #logging.warning('shape file %s created' % shape_file.GetName()) 
            #shape_file.Close()

        

