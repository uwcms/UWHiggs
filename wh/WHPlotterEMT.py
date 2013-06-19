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
from WHPlotterBase import make_styler, parser
import rootpy.plotting.views as views
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors

logging.basicConfig(stream=sys.stderr, level=logging.WARNING)

class WHPlotterEMT(WHPlotterBase.WHPlotterBase):
    def __init__(self):
        super(WHPlotterEMT, self).__init__('EMT')

if __name__ == "__main__":
    plotter = WHPlotterEMT()
    sqrts   = plotter.sqrts
    plotter.defaults['show_charge_fakes'] = True
    options,NOTUSED = parser.parse_args()
    if not options.dry_run:
        ###########################################################################
        ##  Zmm control plots #####################################################
        ###########################################################################

        # Control Z->tautau + jet region
        plotter.plot_mc_vs_data('os/p1p2f3', 'e_m_Mass', rebin=10, xaxis='m_{e#mu} (GeV)', leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-os-p1p2f3-emMass')

        plotter.plot_mc_vs_data('os/p1p2f3', 'nTruePU', rebin=1, xaxis='True PU')
        plotter.add_cms_blurb(sqrts)
        plotter.save('mcdata-os-p1p2f3-nTruePU')

        plotter.plot('Zjets_M50', 'os/p1p2f3/nTruePU', 'nTruePU', rebin=1, xaxis='True PU')
        plotter.save('zjets-os-p1p2f3-nTruePU')


        ## plotter.plot_mc_vs_data('os/p1p2f3', 'bCSVVeto', rebin=1, xaxis='bveto')
        ## plotter.add_cms_blurb(sqrts)
        ## plotter.save('mcdata-os-p1p2f3-bveto')

        plotter.plot_mc_vs_data('os/p1p2f3/w3', 'e_m_Mass', 10)
        plotter.save('mcdata-os-p1p2f3-w3-emMass')

        plotter.plot_mc_vs_data('os/p1f2p3', 'e_m_Mass', 10)
        plotter.save('mcdata-os-p1f2p3-emMass')

        plotter.plot_mc_vs_data('os/f1p2p3', 'e_m_Mass', 10)
        plotter.save('mcdata-os-p1f2p3-emMass')

        # Check PU variables
        #plotter.plot_mc_vs_data('os/p1p2f3', 'rho')
        #plotter.save('mcdata-os-p1p2f3-rho')

        #plotter.plot_mc_vs_data('os/p1p2f3', 'nvtx')
        #plotter.save('mcdata-os-p1p2f3-nvtx')

        # Lower stat but closer to signal region
        #plotter.plot_mc_vs_data('os/p1p2p3', 'rho')
        #plotter.save('mcdata-os-p1p2p3-rho')

        #plotter.plot_mc_vs_data('os/p1p2p3', 'nvtx')
        #plotter.save('mcdata-os-p1p2p3-nvtx')

        # Make Z->mumu + tau jet control

        weighted = plotter.plot('data', 'os/p1p2f3/w3/e_m_Mass',  'hist', rebin=20, styler=make_styler(2, 'hist'), xaxis='m_{e#mu} (GeV)')
        unweighted = plotter.plot('data', 'os/p1p2p3/e_m_Mass', 'same', rebin=20, styler=make_styler(1), xaxis='m_{e#mu} (GeV)')
        weighted.SetTitle('e^{+}#mu^{-} + fake #tau_{h} est.')
        weighted.legendstyle = 'l'
        unweighted.SetTitle('e^{+}#mu^{-} + fake #tau_{h} obs.')
        unweighted.legendstyle = 'pe'
        plotter.add_legend([weighted, unweighted])
        plotter.add_cms_blurb(sqrts)
        plotter.save('ztt-os-fr-control')

        #plotter.plot('data', 'os/p1p2p3/prescale', styler=make_styler(1))
        #plotter.save('ztt-os-prescale-check')

        #plotter.plot('Zjets_M50', 'os/p1p2f3/weight')
        #plotter.save('ztt-mc-event-weights')
        ## Check MC weights
        #plotter.plot('Zjets_M50', 'os/p1p2f3/weight_nopu')
        #plotter.save('ztt-mc-event-weight_nopu')


        ###########################################################################
        ##  FR sideband MC-vs-Data ################################################
        ###########################################################################

        plotter.plot_mc_vs_data('ss/p1f2p3', 'mPt', 5, '#mu_{1} p_{T}', leftside=False)
        plotter.save('mcdata-ss-p1f2p3-mPt')

        plotter.plot_mc_vs_data('ss/p1f2p3', 'subMass', 20, 'Subleading mass (GeV)', leftside=False)
        plotter.save('mcdata-ss-p1f2p3-subMass')

        plotter.plot_mc_vs_data('ss/p1p2f3', 'subMass', 20, 'Subleading mass (GeV)', leftside=False)
        plotter.save('mcdata-ss-p1p2f3-subMass')

        plotter.plot_mc_vs_data('ss/f1p2p3', 'subMass', 20, 'Subleading mass (GeV)', leftside=False)
        plotter.save('mcdata-ss-f1p2p3-subMass')

        plotter.plot_mc_vs_data('ss/p1f2p3/w2', 'mPt', 5, '#mu_{1} p_{T}', leftside=False)
        plotter.save('mcdata-ss-p1f2p3-w2-mPt')

        plotter.plot_mc_vs_data('ss/p1f2p3', 'ePt', 5, 'Electron p_{T}', leftside=False)
        plotter.save('mcdata-ss-p1f2p3-ePt')

        plotter.plot_mc_vs_data('ss/p1f2p3/w2', 'ePt', 5, 'Electron p_{T}', leftside=False)
        plotter.save('mcdata-ss-p1f2p3-w2-ePt')

        plotter.plot_mc_vs_data('ss/f1p2p3', 'ePt', 5, 'Electron p_{T}', leftside=False)
        plotter.save('mcdata-ss-f1p2p3-ePt')

        plotter.plot_mc_vs_data('ss/f1p2p3/w1', 'ePt', 5, 'Electron p_{T}', leftside=False)
        plotter.save('mcdata-ss-f1p2p3-w2-ePt')

        ###########################################################################
        ##  Signal region plots    ################################################
        ###########################################################################

        #BEGIN - New topologicla variables
        plotter.plot_final('LT', 5, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto')
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-LT')
        ## sig_view = plotter.make_signal_views(20, unblinded=(not plotter.blind))
        ## vh_10x   = views.StyleView(
        ##     sig_view['signal120'], 
        ##     **data_styles['VH*']
        ##     )
        ## subm_hist= vh_10x.Get('subMass')
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

        plotter.plot_final('tToMETDPhi', 20, qcd_weight_fraction=0.5, maxy='auto', stack_higgs=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-tToMETDPhi')

        ##   plotter.plot_final('_recoilDaught', 20, qcd_weight_fraction=0.5, maxy='auto', stack_higgs=False)
        ##   plotter.add_cms_blurb(sqrts)
        ##   plotter.save('study-recoilDaught')
        ##
        ##   plotter.plot_final('_recoilWithMet', 20, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto')
        ##   plotter.add_cms_blurb(sqrts)
        ##   plotter.save('study-recoilWithMet')
        ##
        ##   plotter.plot_final('lepRecoil_wMET', 5,  x_range=[0,300], qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto')
        ##   plotter.add_cms_blurb(sqrts)
        ##   plotter.save('study-lepRecoilwMET')
        ##
        ##   plotter.plot_final('lepRecoil', 5,  x_range=[0,300], qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto')
        ##   plotter.add_cms_blurb(sqrts)
        ##   plotter.save('study-lepRecoil')

        plotter.plot_final('type1_pfMetEt'    , 5, qcd_weight_fraction=0.5, stack_higgs=False, maxy='auto', x_range=[0,500])
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-metEt')

        plotter.plot_final('subMass', 20, xaxis='m_{#l_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-subMass-qweight05')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{1}#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-f3-qweight05-werror-subMass')
        #END

        plotter.plot_final('mPt', 10)
        plotter.save('final-mPt')

        plotter.plot_final('ePt', 10)
        plotter.save('final-ePt')

        plotter.plot_final('tPt', 10)
        plotter.save('final-tPt')

        plotter.plot_final('mAbsEta', 10)
        plotter.save('final-mAbsEta')

        plotter.plot_final('eAbsEta', 10)
        plotter.save('final-eAbsEta')

        plotter.plot_final('tAbsEta', 10)
        plotter.save('final-tAbsEta')

        plotter.plot_final('subMass', 20, xaxis='m_{#l_{2}#tau} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass')

        plotter.plot_final('subMass', 20, xaxis='m_{#l_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-qweight05')
        plotter.canvas.SetLogy(True)
        plotter.save('final-subMass-qweight05-logscale')


        plotter.plot_final('subMass', 20, xaxis='m_{#l_{2}#tau} (GeV)', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-werror')

        # Shape only
        plotter.plot_final('subMass', 20, xaxis='m_{#l_{2}#tau} (GeV)', show_error=True,
                           fake_error=0, wz_error=0, zz_error=0)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-wshapeerror')

        ## plotter.plot_final('metSig', 5)
        ## plotter.save('final-metSig')
        plotter.plot_final('tLeadDR', 10)
        plotter.save('final-tLeadDR')
        plotter.plot_final('tSubDR', 10)
        plotter.save('final-tSubDR')

        plotter.plot_final('e_t_Mass', 10)
        plotter.save('final-etMass')

        plotter.plot_final_wz('e_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-wz-etMass')

        plotter.plot_final_wz('mPt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', maxy=20)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-wz-mPt')

        #plotter.plot_final_wz('mJetPt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-wz-mJetPt')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{1}#tau_{#mu}} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{1}#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-qweight05-subMass')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{1}#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-qweight05-werror-subMass')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{1}#tau_{#mu}} (GeV)', show_error=True,
                           fake_error=0, wz_error=0, zz_error=0)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-wshaperror-subMass')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{1}#tau_{#mu}} (GeV)', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-werror-subMass')

        plotter.plot_final_f3('e_t_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-etMass')

        plotter.plot_final_f3('e_m_Mass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-emMass')

        #plotter.plot_final_f3('subMass', 20, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', qcd_correction=True)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-f3q-subMass')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{#l_{1}#tau_{#mu}} (GeV)', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-werror')

        plotter.plot_final_f3('mPt', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-mPt')

        plotter.plot_final_f3('ePt', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-ePt')

        ## plotter.plot_final_f3('eChargeIdMedium', 1, xaxis='Charge ID Med', maxy=None)
        ## plotter.add_cms_blurb(sqrts)
        ## plotter.save('final-f3-eChargeIdMedium')

        plotter.plot_final_f3('eChargeIdTight', 1, xaxis='Charge ID Tight', maxy=None)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-eChargeIdTight')

        #plotter.plot_final_f3('mJetPt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-f3-mJetPt')
        #END OF if not options.dry_run:
    ###########################################################################
    ##  Making shape file     #################################################
    ###########################################################################
    prefixes = [options.prefix+'$'] if options.prefix else ['']
    prefixes = [i+'$' for i in options.prefixes.split(',') if i] if options.prefixes else prefixes
    for prefix in prefixes:
        shape_prefix = prefix if len(prefixes) > 1 else ''
        shape_prefix = shape_prefix.replace(':','_').replace('$','_')
        
        plotter.plot_final(prefix+'LT', 5, xaxis='LT (GeV)', qcd_weight_fraction=0.5)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-%s-LT-qweight05' % shape_prefix)

        plotter.plot_final(prefix+'subMass', 20, xaxis='m_{#l_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-%s-subMass-qweight05' % shape_prefix)

        plotter.plot_final_f3(prefix+'subMass', 20, xaxis='m_{l_{1}#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-%s-f3-qweight05-werror-subMass' % shape_prefix)

        print "saving shape file %s" % ('%semt_shapes_%s.root' % (shape_prefix, plotter.period))
        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, '%semt_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
        shape_dir = shape_file.mkdir('emt')
        plotter.write_shapes(prefix+'subMass', 20, shape_dir, qcd_fraction=0.5)
        shape_dir = shape_file.mkdir('emt_w')
        plotter.write_shapes(prefix+'subMass', 20, shape_dir, qcd_fraction=0.0)
        shape_dir = shape_file.mkdir('emt_q')
        plotter.write_shapes(prefix+'subMass', 20, shape_dir, qcd_fraction=1.0)
        #plotter.write_cut_and_count('subMass', shape_dir, unblinded=True)
        shape_file.Close()
