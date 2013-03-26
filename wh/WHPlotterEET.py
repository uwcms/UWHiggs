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
from WHPlotterBase import make_styler
import rootpy.plotting.views as views
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors

logging.basicConfig(stream=sys.stderr, level=logging.WARNING)

class WHPlotterEET(WHPlotterBase.WHPlotterBase):
    def __init__(self):
        obj1_charge_mapper={"e1_e2_Mass":"e*1_e2_Mass","e1_t_Mass":"e*1_t_Mass"}
        obj2_charge_mapper={"e1_e2_Mass":"e1_e*2_Mass","e2_t_Mass":"e*2_t_Mass"}
        super(WHPlotterEET, self).__init__('EET', obj1_charge_mapper, obj2_charge_mapper)

if __name__ == "__main__":
    plotter = WHPlotterEET()
    sqrts   = plotter.sqrts
    plotter.defaults['show_charge_fakes'] = True

    ###########################################################################
    ##  Zmm control plots #####################################################
    ###########################################################################

    # Control Z->mumu + jet region
    plotter.plot_mc_vs_data('os/p1p2f3', 'e1_e2_Mass', xaxis='m_{ee} (GeV)', xrange=(60, 120))
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-os-p1p2f3-e1e2Mass')

    plotter.plot_mc_vs_data('os/p1p2f3/w3', 'e1_e2_Mass')
    plotter.save('mcdata-os-p1p2f3-w3-e1e2Mass')

    plotter.plot_mc_vs_data('os/p1f2p3', 'e1_e2_Mass', xaxis='m_{ee} (GeV)', xrange=(60, 120))
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-os-p1f2p3-e1e2Mass')

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

    #
    # Makes Tau fr control plot
    #
    zmm_weighted = plotter.plot('data', 'os/p1p2f3/w3/e1_e2_Mass',  'hist', styler=make_styler(2, 'hist'), xrange=(60, 120))
    zmm_weighted.SetTitle("Zee + fake #tau_{h} est.")
    zmm_weighted.legendstyle='l'
    zmm_weighted.GetXaxis().SetTitle("m_{ee} (GeV)")

    zmm_unweighted = plotter.plot('data', 'os/p1p2p3/e1_e2_Mass', 'same', styler=make_styler(1), xrange=(60, 120))
    zmm_unweighted.SetTitle("Zee observed")
    zmm_unweighted.SetTitle("Zee + fake #tau_{h} obs.")
    zmm_unweighted.legendstyle='pe'

    plotter.add_legend([zmm_weighted, zmm_unweighted])
    plotter.add_cms_blurb(sqrts)
    plotter.save('zmm-os-fr-control')

    #
    # Makes charge fr control plot
    #

    zeet_os_weighted = plotter.plot('data', 'os/p1p2f3/c1/e1_e2_Mass',  'hist', styler=make_styler(2, 'hist'), xrange=(60, 120))
    zeet_os_weighted.SetTitle("Ze^{#pm}e^{#mp} + fake #tau_{h} charge flip est.")
    zeet_os_weighted.legendstyle='l'
    zeet_os_weighted.GetXaxis().SetTitle("M_{ee} (GeV)")

    zee_ss_unweighted = plotter.plot('data', 'ss/p1p2f3/e1_e2_Mass', 'same', styler=make_styler(1), xrange=(60, 120))
    zee_ss_unweighted.SetTitle("Ze^{#pm}e^{#pm} + fake #tau_{h} obs.")
    zee_ss_unweighted.legendstyle='pe'

    plotter.add_legend([zeet_os_weighted, zee_ss_unweighted])
    plotter.add_cms_blurb(sqrts)
    plotter.save('zee-os-ss-charge-flip-control')

    plotter.plot('Zjets_M50', 'os/p1p2f3/weight')
    plotter.save('zee-mc-event-weights')
    # Check MC weights
    ## plotter.plot('Zjets_M50', 'os/p1p2f3/weight_nopu')
    ## plotter.save('zee-mc-event-weight_nopu')


    ###########################################################################
    ##  FR sideband MC-vs-Data ################################################
    ###########################################################################

    plotter.plot_mc_vs_data('ss/p1f2p3', 'e1Pt', rebin=10, xaxis='e_{1} p_{T} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-e1Pt')

    plotter.plot_mc_vs_data('ss/p1f2p3', 'e2_t_Mass', rebin=10, xaxis='m_{e2#tau} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-subMass')

    plotter.plot_mc_vs_data('ss/p1f2p3/w2', 'e1Pt', rebin=10, xaxis='e_{1} p_{T}', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1f2p3-w2-e1Pt')

    plotter.plot_mc_vs_data('ss/f1p2p3', 'e2_t_Mass', rebin=20, xaxis='m_{e2#tau} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-f1p2p3-subMass')

    plotter.plot_mc_vs_data('ss/f1p2p3/w1', 'e2_t_Mass', rebin=20, xaxis='m_{e2#tau} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-f1p2p3-w1-subMass')

    plotter.plot_mc_vs_data('ss/p1p2f3', 'e1_e2_Mass', rebin=10, xaxis='m_{ee} (GeV)', leftside=False)
    plotter.add_cms_blurb(sqrts)
    plotter.save('mcdata-ss-p1p2f3-e1e2Mass')



    ###########################################################################
    ##  Signal region plots    ################################################
    ###########################################################################

    #BEGIN - New topologicla variables
    ## sig_view = plotter.make_signal_views(20, unblinded=(not plotter.blind))
    ## vh_10x   = views.StyleView(
    ##     sig_view['signal120'], 
    ##     **data_styles['VH*']
    ##     )
    ## subm_hist= vh_10x.Get('e2_t_Mass')
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

    plotter.plot_final('LT', 5, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-LT')

    plotter.plot_final('_recoilDaught', 3, qcd_weight_fraction=0.5, maxy='auto', stack_higgs=False, x_range=[0,300], show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-recoilDaught')

    rebin=3
    plotter.plot_final_f3('_recoilDaught', rebin, qcd_weight_fraction=0.5, maxy='auto', x_range=[0,300], show_error=True)
    sig_view = plotter.make_signal_views(  rebin, unblinded=(not plotter.blind))
    vh_10x   = views.TitleView(
        views.StyleView(
            views.ScaleView(sig_view['signal120'], 120),
            **data_styles['VH*']
            ),
            "(20#times) m_{H} = 125"
        )
    sign_hist= vh_10x.Get('_recoilDaught')
    sign_hist.Draw('same')
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-recoilDaught')

    rebin=20
    plotter.plot_final_f3('tToMETDPhi',  rebin, qcd_weight_fraction=0.5, maxy='auto', show_error=True)
    sig_view = plotter.make_signal_views(rebin, unblinded=(not plotter.blind))
    vh_10x   = views.TitleView(
        views.StyleView(
            views.ScaleView(sig_view['signal120'], 120),
            **data_styles['VH*']
            ),
            "(20#times) m_{H} = 125"
        )
    sign_hist= vh_10x.Get('tToMETDPhi')
    sign_hist.Draw('same')
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-tToMETDPhi')

    plotter.plot_final('e1_e2_Mass', 10, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-e1_e2_Mass')

    plotter.plot_final('mass', 10, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-mass')

    plotter.plot_final('pt_ratio', 10, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-pt_ratio')

    plotter.plot_final('metEt'    , 5, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto', x_range=[0,500], show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-metEt')

    plotter.plot_final_f3('e1Pt',  10, qcd_weight_fraction=0.5, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-e1Pt')

    plotter.plot_final_f3('e2Pt',  10, qcd_weight_fraction=0.5, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-e2Pt')

    plotter.plot_final_f3('e1_e2_Mass', 2, qcd_weight_fraction=0.5, show_error=True, maxy=80, x_range=[60,120])
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-e1_e2_Mass')

    plotter.plot_final_f3('e1eta_on_z_peak', 10, qcd_weight_fraction=0.5, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-e1eta_on_z_peak')

    plotter.plot_final_f3('e1pt_on_z_peak', 10, qcd_weight_fraction=0.5, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-e1pt_on_z_peak')

    plotter.plot_final_f3('e2eta_on_z_peak', 10, qcd_weight_fraction=0.5, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-e2eta_on_z_peak')

    plotter.plot_final_f3('e2pt_on_z_peak', 10, qcd_weight_fraction=0.5, maxy='auto', show_error=True)
    plotter.add_cms_blurb(sqrts)
    plotter.save('study-f3-e2pt_on_z_peak')    
    ## #END

    plotter.plot_final('e1Pt', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-e1Pt')

    plotter.plot_final('e2Pt', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-e2Pt')

    plotter.plot_final('tPt', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-tPt')

    #plotter.plot_final('e1AbsEta', 10)
    #plotter.add_cms_blurb(sqrts)
    #plotter.save('final-e1AbsEta')

    #plotter.plot_final('e2AbsEta', 10)
    #plotter.add_cms_blurb(sqrts)
    #plotter.save('final-e2AbsEta')

    plotter.plot_final('tAbsEta', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-tAbsEta')

    plotter.plot_final('e2_t_Mass', 20, xaxis='m_{e_{2}#tau} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-subMass')

    plotter.plot_final('e1_t_Mass', 20, xaxis='m_{e_{1}#tau} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-e1tMass')

    plotter.plot_final('e1_e2_Mass', 20, xaxis='m_{ee} (GeV)')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-e1e2Mass')

    plotter.plot_final('tAbsEta', 5, xaxis='|#eta_#tau|')
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-tAbsEta')

    plotter.plot_final('e2RelPFIsoDB', 10)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-e2Iso')

    ## plotter.plot_final('metSignificance', 5)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('final-metSig')

    plotter.plot_final('LT', 5)
    plotter.add_cms_blurb(sqrts)
    plotter.save('final-LT')


    ###########################################################################
    ##  Making shape file     #################################################
    ###########################################################################

    shape_file = ROOT.TFile(
        os.path.join(plotter.outputdir, 'eet_shapes_%s.root' % plotter.period), 'RECREATE')
    shape_dir = shape_file.mkdir('eet')
    plotter.write_shapes('e2_t_Mass', 20, shape_dir)
    shape_file.Close()


