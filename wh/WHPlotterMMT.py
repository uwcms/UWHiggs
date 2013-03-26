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
from WHPlotterBase import make_styler
import rootpy.plotting.views as views
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

class WHPlotterMMT(WHPlotterBase.WHPlotterBase):
    def __init__(self):
        super(WHPlotterMMT, self).__init__('MMT')

if __name__ == "__main__":
    plotter = WHPlotterMMT()
    sqrts   = plotter.sqrts

    ###########################################################################
    ##  Zmm control plots #####################################################
    ###########################################################################

    # Control Z->mumu + jet region
    plotter.plot_mc_vs_data('os/p1p2f3', 'm1_m2_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(60, 120))
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-os-p1p2f3-m1m2Mass')

    plotter.plot_mc_vs_data('os/p1p2p3', 'm1_m2_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(60, 120))
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-os-p1p2p3-m1m2Mass')

    plotter.plot_mc_vs_data('ss/p1p2p3_enhance_wz', 'm2_t_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1p2p3-enhance_wz-subMass')

    plotter.plot_mc_vs_data('ss/p1p2p3_enhance_wz', 'm1_t_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1p2p3-enhance_wz-leadMass')

    plotter.plot_mc_vs_data('ss/p1f2p3_enhance_wz', 'm1_t_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-enhance_wz-leadMass')

    plotter.plot_mc_vs_data('ss/p1f2p3_enhance_wz/w2', 'm1_t_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-w2-enhance_wz-leadMass')

    plotter.plot_mc_vs_data('ss/p1p2p3_enhance_wz', 'subMTMass', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1p2p3-enhance_wz-subMTMass')

    plotter.plot_mc_vs_data('ss/p1p2p3_enhance_wz', 'm2Pt', xaxis='m_{#mu#mu} (GeV)', xrange=(0, 120), rebin=5)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1p2p3-enhance_wz-m2Pt')

    plotter.compare_shapes('Zjets_M50', 'data', 'os/p1p2f3/nvtx')
    plotter.save('z-vs-data-nvtx-shape')
    plotter.compare_shapes('Zjets_M50', 'data', 'os/p1p2f3/rho')
    plotter.save('z-vs-data-rho-shape')

    plotter.plot_mc_vs_data('os/p1p2f3/w3', 'm1_m2_Mass')
    plotter.save('mcdata-os-p1p2f3-w3-m1m2Mass')

    plotter.plot_mc_vs_data('os/p1f2p3', 'm1_m2_Mass', xaxis='m_{#mu#mu} (GeV)', xrange=(60, 120))
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-os-p1f2p3-m1m2Mass')

    # Check PU variables
    plotter.plot_mc_vs_data('os/p1p2f3', 'rho')
    plotter.save('mcdata-os-p1p2f3-rho')

    plotter.plot_mc_vs_data('os/p1p2f3', 'nvtx')
    plotter.save('mcdata-os-p1p2f3-nvtx')

    # Lower stat but closer to signal region
    plotter.plot_mc_vs_data('os/p1p2p3', 'rho')
    plotter.save('mcdata-os-p1p2p3-rho')

    plotter.plot_mc_vs_data('os/p1p2p3', 'nvtx')
    plotter.save('mcdata-os-p1p2p3-nvtx')

    # Make Z->mumu + tau jet control

    antiiso_m2JetPt = plotter.plot('data', 'ss/p1f2p3/m2JetPt',  'hist', styler=make_styler(2, 'hist'), xrange=(0, 120), rebin=10)
    antiiso_m2JetPt.SetTitle("Anti-iso CR yield")
    antiiso_m2JetPt.legendstyle='l'
    antiiso_m2JetPt.GetXaxis().SetTitle("#mu_{2} Jet Pt")
    plotter.save('data-p1f2p3-m2JetPt')

    antiiso_m1JetPt = plotter.plot('data', 'ss/f1p2p3/m1JetPt',  'hist', styler=make_styler(2, 'hist'), xrange=(0, 120), rebin=10)
    antiiso_m1JetPt.SetTitle("Anti-iso CR yield")
    antiiso_m1JetPt.legendstyle='l'
    antiiso_m1JetPt.GetXaxis().SetTitle("#mu_{1} Jet Pt")
    plotter.save('data-f1p2p3-m1JetPt')

    zmm_weighted = plotter.plot('data', 'os/p1p2f3/w3/m1_m2_Mass',  'hist', styler=make_styler(2, 'hist'), xrange=(60, 120))
    zmm_weighted.SetTitle("Z#mu#mu + fake #tau_{h} est.")
    zmm_weighted.legendstyle='l'
    zmm_weighted.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")

    zmm_unweighted = plotter.plot('data', 'os/p1p2p3/m1_m2_Mass', 'same', styler=make_styler(1), xrange=(60, 120))
    zmm_unweighted.SetTitle("Z#mu#mu observed")
    zmm_unweighted.SetTitle("Z#mu#mu + fake #tau_{h} obs.")
    zmm_unweighted.legendstyle='pe'

    plotter.add_legend([zmm_weighted, zmm_unweighted])
    plotter.add_cms_blurb(sqrts)
    plotter.save('zmm-os-fr-control')

    ## plotter.plot('data', 'os/p1p2p3/prescale', styler=make_styler(1))
    ## plotter.save('zmm-os-prescale-check')

    plotter.plot('Zjets_M50', 'os/p1p2f3/weight')
    plotter.save('zmm-mc-event-weights')
    # Check MC weights
    ## plotter.plot('Zjets_M50', 'os/p1p2f3/weight_nopu')
    ## plotter.save('zmm-mc-event-weight_nopu')

    plotter.plot('Zjets_M50', 'os/p1p2f3/nTruePU', 'nTruePU', rebin=1, xaxis='True PU')
    plotter.save('zjets-os-p1p2f3-nTruePU')


    ###########################################################################
    ##  FR sideband MC-vs-Data ################################################
    ###########################################################################

    plotter.plot_mc_vs_data('ss/p1f2p3', 'm1Pt', rebin=10, xaxis='#mu_{1} p_{T} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-m1Pt')

    plotter.plot_mc_vs_data('ss/p1f2p3', 'm2_t_Mass', rebin=10, xaxis='m_{#mu2#tau} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-subMass')

    plotter.plot_mc_vs_data('ss/p1f2p3/w2', 'm1Pt', rebin=10, xaxis='#mu_{1} p_{T}', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-w2-m1Pt')

    plotter.plot_mc_vs_data('ss/f1p2p3', 'm2_t_Mass', rebin=20, xaxis='m_{#mu2#tau} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-f1p2p3-subMass')

    plotter.plot_mc_vs_data('ss/f1p2p3/w1', 'm2_t_Mass', rebin=20, xaxis='m_{#mu2#tau} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-f1p2p3-w1-subMass')

    plotter.plot_mc_vs_data('ss/p1f2f3', 'm2AbsEta', rebin=10, xaxis='m_{#mu2#tau} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2f3-m2AbsEta')

    plotter.plot_mc_vs_data('ss/p1f2p3', 'm2AbsEta', rebin=10, xaxis='m_{#mu2#tau} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-m2AbsEta')



    ###########################################################################
    ##  Signal region plots    ################################################
    ###########################################################################

    #BEGIN - New topologicla variables
    ## sig_view = plotter.make_signal_views(20, unblinded=(not plotter.blind))
    ## vh_10x   = views.StyleView(
    ##     sig_view['signal120'],
    ##     **data_styles['VH*']
    ##     )
    ## subm_hist= vh_10x.Get('m2_t_Mass')
    ## subm_hist.SetTitle('subleading Mass')
    ## subm_hist.Draw()
    ## true_hist= vh_10x.Get('true_mass')
    ## true_hist.SetTitle('true Mass')
    ## true_hist.SetLineColor(2)
    ## true_hist.Draw('same')
    ## subm_hist.GetYaxis().SetRangeUser(0,max([subm_hist.GetMaximum(),true_hist.GetMaximum()])*1.2)
    ## legend = plotter.add_legend([subm_hist,true_hist], leftside=False)
    ## legend.Draw()
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('study-trueMass')

    plotter.plot_final('LT', 5, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto')
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-LT')

    plotter.plot_final('tToMETDPhi', 20, qcd_weight_fraction=0.5, maxy='auto', stack_higgs=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-tToMETDPhi')

    plotter.plot_final('_recoilDaught', 20, qcd_weight_fraction=0.5, maxy='auto', stack_higgs=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-recoilDaught')

    plotter.plot_final('_recoilWithMet', 20, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto')
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-recoilWithMet')

    plotter.plot_final('lepRecoil_wMET', 5,  x_range=[0,300], qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto')
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-lepRecoilwMET')

    plotter.plot_final('lepRecoil', 5,  x_range=[0,300], qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto')
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-lepRecoil')

    plotter.plot_final('metEt'    , 5, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto', x_range=[0,500])
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-metEt')

    plotter.plot_final('m2_t_Mass', 20, xaxis='m_{#mu_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-qweight05-subMass')

    plotter.plot_final_f3('m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0.5, show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-qweight05-subMass')
    #END

    plotter.plot_final('m1Pt', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-m1Pt')

    plotter.plot_final('m2Pt', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-m2Pt')

    plotter.plot_final('m2Pt', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-m2Pt')

    plotter.plot_final('m1AbsEta', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-m1AbsEta')

    plotter.plot_final('m2AbsEta', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-m2AbsEta')

    plotter.plot_final('m2AbsEta', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-m2AbsEta')

    plotter.plot_final('m2_t_Mass', 20, xaxis='m_{#mu_{2}#tau} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-subMass')

    plotter.plot_final('m2_t_Mass', 20, xaxis='m_{#mu_{2}#tau} (GeV)', qcd_weight_fraction=1)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-qweight-subMass')

    plotter.plot_final('m2_t_Mass', 20, xaxis='m_{#mu_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-qweight05-subMass')

    plotter.plot_final('m2_t_Mass', 20, xaxis='m_{#mu_{2}#tau} (GeV)', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-subMass-werror')

    plotter.plot_final('m2_t_Mass', 20, xaxis='m_{#mu_{2}#tau} (GeV)',
                       show_error=True, fake_error=0, wz_error=0, zz_error=0)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-subMass-wshapeerror')

    plotter.plot_final('m2RelPFIsoDB', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-m2Iso')

    plotter.plot_final_wz('m1_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-wz-leadMass')

    plotter.plot_final_wz('m2Pt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-wz-m2Pt')

    plotter.plot_final_wz('m2JetPt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-wz-m2JetPt')

    plotter.plot_final_f3('m1_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-leadMass')

    plotter.plot_final_f3('m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-subMass')

    plotter.plot_final_f3('m2JetBtag', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-m2JetBtag')

    plotter.plot_final_f3('m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3q-subMass')

    plotter.plot_final_f3('m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=1)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-qweight-subMass')

    plotter.plot_final_f3('m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0.5)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-qweight05-subMass')

    plotter.plot_final_f3('m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0.5, show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-qweight05-werror-subMass')

    plotter.plot_final_f3('m2_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=False, qcd_weight_fraction=0.5,
                          show_error=True, fake_error=0, wz_error=0, zz_error=0)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-qweight05-wshapeerror-subMass')

    plotter.plot_final_f3_split('m2_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-split-subMass')

    plotter.plot_final_f3('m2_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-subMass-werror')

    plotter.plot_final_f3('m2Pt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-m2Pt')

    plotter.plot_final_f3('m2AbsEta', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-m2AbsEta')

    plotter.plot_final_f3('m2AbsEta', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3q-m2AbsEta')

    plotter.plot_final_f3('m2JetPt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-f3-m2JetPt')

    ###########################################################################
    ##  Check QCD contamination in control regions ############################
    ###########################################################################

    plotter.plot_qcd_contamination('m2_t_Mass', 2, 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-qcd2-subMass')

    plotter.plot_qcd_contamination('m2_t_Mass', 1, 20)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-qcd1-subMass')

    plotter.plot_qcd_contamination('m2JetPt', 2, 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-qcd2-m2JetPt')

    plotter.plot_qcd_contamination('m1JetPt', 1, 20)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-qcd1-m1JetPt')

    ###########################################################################
    ##  Making shape file     #################################################
    ###########################################################################

    shape_file = ROOT.TFile(
        os.path.join(plotter.outputdir, 'mmt_shapes_%s.root' % plotter.period), 'RECREATE')
    shape_dir = shape_file.mkdir('mmt')
    plotter.write_shapes('m2_t_Mass', 20, shape_dir, unblinded=True, qcd_fraction=0.5)
    shape_dir = shape_file.mkdir('mmt_w')
    plotter.write_shapes('m2_t_Mass', 20, shape_dir, unblinded=True, qcd_fraction=0.0)
    shape_dir = shape_file.mkdir('mmt_q')
    plotter.write_shapes('m2_t_Mass', 20, shape_dir, unblinded=True, qcd_fraction=1.0)
    #plotter.write_cut_and_count('subMass', shape_dir, unblinded=True)
    shape_file.Close()
