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
        obj2_charge_mapper={
            "subMass#LT"          : "subMass*#LT"         ,
            "subMass#tPt"         : "subMass*#tPt"        ,
            "e_t_Mass"            : "e*_t_Mass"           ,
            "e_m_Mass"            : "e*_m_Mass"           ,
            "subMass"             : "subMass*"            ,} 
        super(WHPlotterEMT, self).__init__('EMT', {}, obj2_charge_mapper)

if __name__ == "__main__":
    plotter = WHPlotterEMT()
    sqrts   = plotter.sqrts
    plotter.defaults['show_charge_fakes'] = True
    options,NOTUSED = parser.parse_args()
    if not options.dry_run:
        ###########################################################################
        ##  Zmm control plots #####################################################
        ###########################################################################
        plotter.set_subdir('mc_data')

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
        plotter.set_subdir('')
        rebin_slim = [0,20]+range(30, 81, 10)+[100,130,200]

        plotter.plot_final('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, project=[0, 130], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTLow')

        plotter.plot_final('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTHigh')

        plotter.plot_final('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, project=[80, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTCut')

        plotter.plot_final('subMass#LT', 20, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, project=[80, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTCut-flatbin')

        plotter.plot_final('LT', 5, xaxis='LT (GeV)', maxy=15)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-LT')

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

        plotter.plot_final('subMass', 20, xaxis='m_{l_{2}#tau} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass')

        plotter.plot_final('subMass', 20, xaxis='m_{l_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-qweight05')
        plotter.canvas.SetLogy(True)
        plotter.save('final-subMass-qweight05-logscale')


        plotter.plot_final('subMass', 20, xaxis='m_{l_{2}#tau} (GeV)', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-werror')

        # Shape only
        plotter.plot_final('subMass', 20, xaxis='m_{l_{2}#tau} (GeV)', show_error=True,
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

        plotter.plot_final('subMass#LT', 20, xaxis='subleading mass from projection', maxy=None, project=[0, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMassProj')

        plotter.plot_final('subMass#LT', 20, xaxis='M_{l_{2}#tau} (GeV)', maxy=None, project=[0, 90], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTLo')

        plotter.plot_final('subMass#LT', 20, xaxis='M_{l_{2}#tau} (GeV)', maxy=None, project=[90, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-subMass-LTHi')

        ###########################################################################
        ##  WZ enhanced region plots    ###########################################
        ###########################################################################
        plotter.set_subdir('WZ_enhanced')

        plotter.plot_final_wz('e_t_Mass', 10, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-wz-etMass')

        plotter.plot_final_wz('mPt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)', maxy=20)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-wz-mPt')

        #plotter.plot_final_wz('mJetPt', 5, xaxis='m_{#mu_{1}#tau_{#mu}} (GeV)')
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-wz-mJetPt')

        ###########################################################################
        ##  F3 enhanced region plots    ###########################################
        ###########################################################################
        plotter.set_subdir('f3')

        plotter.plot_final_f3('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, project=[0, 130], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-LTLow')

        plotter.plot_final_f3('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, project=[130, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-LTHigh')

        plotter.plot_final_f3('subMass#LT', rebin_slim, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, project=[80, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-LTCut')

        plotter.plot_final_f3('subMass#LT', 20, xaxis='m_{l_{2}#tau} (GeV)', maxy=None, project=[80, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-LTCut-flatbin')

        plotter.plot_final_f3('LT', 5, xaxis='LT (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-LT')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{2}#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass')

        plotter.plot_final_f3('subMass', 200, xaxis='m_{l_{2}#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-counting-like')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{2}#tau_{#mu}} (GeV)', qcd_weight_fraction=1, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-qcdfake-subMass')

        plotter.plot_final_f3('subMass', 20, xaxis='m_{l_{2}#tau_{#mu}} (GeV)', qcd_weight_fraction=0, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-wjetfake-subMass')

        plotter.plot_final_f3('e_t_Mass', 20, xaxis='m_{e#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-etMass')

        plotter.plot_final_f3('e_m_Mass', 20, xaxis='m_{e#mu} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-emMass')

        #plotter.plot_final_f3('m_t_Mass', 20, xaxis='m_{#mu#tau} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        #plotter.add_cms_blurb(sqrts)
        #plotter.save('final-f3-emMass')

        plotter.plot_final_f3('mPt', 10, xaxis='p_{T#mu} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-mPt')

        plotter.plot_final_f3('ePt', 10, xaxis='p_{Te} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-ePt')

        plotter.plot_final_f3('tPt', 10, xaxis='p_{T#tau} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-tPt')

        plotter.plot_final_f3('tSubDR', 10, xaxis='#DeltaR_{l_{2}#tau}', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-tSubDR')

        plotter.plot_final_f3('tLeadDR', 10, xaxis='#DeltaR_{l_{1}#tau}', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-tLeadDR')

        plotter.plot_final_f3('e_t_DR', 10, xaxis='#DeltaR_{e#tau}', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-etDR')

        plotter.plot_final_f3('m_t_DR', 10, xaxis='#DeltaR_{#mu#tau}', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-mtDR')

        plotter.plot_final_f3('subPt', 10, xaxis='p_{Tl_{2}} (GeV)', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subPt')

        plotter.plot_final_f3('subJetPt', 10, xaxis='p_{TJetl_{2}} (GeV)', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subJetPt')

        plotter.plot_final_f3('leadPt', 10, xaxis='p_{Tl_{1}} (GeV)', show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-leadPt')

        plotter.plot_final_f3('eChargeIdTight', 1, xaxis='Charge ID Tight', maxy=None)
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-eChargeIdTight')

        plotter.plot_final_f3('subMass#LT', 20, xaxis='subleading mass from projection', maxy=None, project=[0, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMassProj')

        plotter.plot_final_f3('subMass#LT', 20, xaxis='M_{l_{2}#tau} (GeV)', maxy=None, project=[0, 90], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-LTLo')

        plotter.plot_final_f3('subMass#LT', 20, xaxis='M_{l_{2}#tau} (GeV)', maxy=None, project=[90, 650], project_axis='X')
        plotter.add_cms_blurb(sqrts)
        plotter.save('final-f3-subMass-LTHi')

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
        
        plotter.plot_final(prefix+'LT', 5, xaxis='LT (GeV)', qcd_weight_fraction=0.5)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-%s-LT-qweight05' % shape_prefix)

        plotter.plot_final(prefix+'subMass', 20, xaxis='m_{#l_{2}#tau} (GeV)', qcd_weight_fraction=0.5)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-%s-subMass-qweight05' % shape_prefix)

        plotter.plot_final_f3(prefix+'subMass', 20, xaxis='m_{l_{1}#tau_{#mu}} (GeV)', qcd_weight_fraction=0.5, show_error=True)
        plotter.add_cms_blurb(sqrts)
        plotter.save('study-%s-f3-qweight05-werror-subMass' % shape_prefix)

        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, 'LTCut_emt_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
        shape_dir = shape_file.mkdir('emtCatHigh')
        plotter.write_shapes(prefix+'subMass#LT', 20, shape_dir, qcd_fraction=0.5, project=[80, 650], project_axis='X')
        shape_dir = shape_file.mkdir('emtCatHigh_w')
        plotter.write_shapes(prefix+'subMass#LT', 20, shape_dir, qcd_fraction=0.0, project=[80, 650], project_axis='X')
        shape_dir = shape_file.mkdir('emtCatHigh_q')
        plotter.write_shapes(prefix+'subMass#LT', 20, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')
        shape_file.Close()

        rebin = rebin_slim
        shape_file = ROOT.TFile(
            os.path.join(plotter.outputdir, 'LTCut_RebinSmart_emt_shapes_%s.root' % ( plotter.period) ), 'RECREATE')
        shape_dir = shape_file.mkdir('emtCatHigh')
        plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=0.5, project=[80, 650], project_axis='X')
        shape_dir = shape_file.mkdir('emtCatHigh_w')
        plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=0.0, project=[80, 650], project_axis='X')
        shape_dir = shape_file.mkdir('emtCatHigh_q')
        plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=1.0, project=[80, 650], project_axis='X')
        shape_file.Close()

        low = 0
        rebin_fat  = range(0, 80, 20)+[110,200]
        for mid in range(70, 140, 10):
            for name, rebin in zip(['slim', 'fat'],[rebin_slim, rebin_fat]):
                shape_prefix = "Rebin%s_LT_CAT_%i_%i_inf_" % (name, low, mid)
                
                shape_file = ROOT.TFile(
                    os.path.join(plotter.outputdir, '%semt_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
        
                shape_dir = shape_file.mkdir('emtCatLow')
                plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=0.5, project=[low, mid], project_axis='X')
                shape_dir = shape_file.mkdir('emtCatLow_w')
                plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=0.0, project=[low, mid], project_axis='X')
                shape_dir = shape_file.mkdir('emtCatLow_q')
                plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=1.0, project=[low, mid], project_axis='X')
        
                shape_dir = shape_file.mkdir('emtCatHigh')
                plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=0.5, project=[mid, 650], project_axis='X')
                shape_dir = shape_file.mkdir('emtCatHigh_w')
                plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=0.0, project=[mid, 650], project_axis='X')
                shape_dir = shape_file.mkdir('emtCatHigh_q')
                plotter.write_shapes(prefix+'subMass#LT', rebin, shape_dir, qcd_fraction=1.0, project=[mid, 650], project_axis='X')
        
                logging.warning('shape file %s created' % shape_file.GetName()) 
                shape_file.Close()


        for mid in [70, 80]:
            shape_prefix = "Rebin%s_LT_CAT_%i_%i_inf_" % ('smart', low, mid)
            shape_file = ROOT.TFile(
                os.path.join(plotter.outputdir, '%semt_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
            
            shape_dir = shape_file.mkdir('emtCatLow')
            plotter.write_shapes(prefix+'subMass#LT', rebin_fat, shape_dir, qcd_fraction=0.5, project=[low, mid], project_axis='X')
            shape_dir = shape_file.mkdir('emtCatLow_w')
            plotter.write_shapes(prefix+'subMass#LT', rebin_fat, shape_dir, qcd_fraction=0.0, project=[low, mid], project_axis='X')
            shape_dir = shape_file.mkdir('emtCatLow_q')
            plotter.write_shapes(prefix+'subMass#LT', rebin_fat, shape_dir, qcd_fraction=1.0, project=[low, mid], project_axis='X')
        
            shape_dir = shape_file.mkdir('emtCatHigh')
            plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.5, project=[mid, 650], project_axis='X')
            shape_dir = shape_file.mkdir('emtCatHigh_w')
            plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=0.0, project=[mid, 650], project_axis='X')
            shape_dir = shape_file.mkdir('emtCatHigh_q')
            plotter.write_shapes(prefix+'subMass#LT', rebin_slim, shape_dir, qcd_fraction=1.0, project=[mid, 650], project_axis='X')
        
            logging.warning('shape file %s created' % shape_file.GetName()) 
            shape_file.Close()


        #low = 10
        #for mid in range(25, 42, 5):
        #    shape_prefix = "TPT_CAT_%i_%i_inf_" % (20, mid)
        #    
        #    shape_file = ROOT.TFile(
        #        os.path.join(plotter.outputdir, '%semt_shapes_%s.root' % (shape_prefix, plotter.period) ), 'RECREATE')
        #    
        #    shape_dir = shape_file.mkdir('emtCatLow')
        #    plotter.write_shapes(prefix+'subMass#tPt', 20, shape_dir, qcd_fraction=0.5, project=[low, mid], project_axis='Y')
        #    shape_dir = shape_file.mkdir('emtCatLow_w')
        #    plotter.write_shapes(prefix+'subMass#tPt', 20, shape_dir, qcd_fraction=0.0, project=[low, mid], project_axis='Y')
        #    shape_dir = shape_file.mkdir('emtCatLow_q')
        #    plotter.write_shapes(prefix+'subMass#tPt', 20, shape_dir, qcd_fraction=1.0, project=[low, mid], project_axis='Y')
        #
        #    shape_dir = shape_file.mkdir('emtCatHigh')
        #    plotter.write_shapes(prefix+'subMass#tPt', 20, shape_dir, qcd_fraction=0.5, project=[mid, 200], project_axis='Y')
        #    shape_dir = shape_file.mkdir('emtCatHigh_w')
        #    plotter.write_shapes(prefix+'subMass#tPt', 20, shape_dir, qcd_fraction=0.0, project=[mid, 200], project_axis='Y')
        #    shape_dir = shape_file.mkdir('emtCatHigh_q')
        #    plotter.write_shapes(prefix+'subMass#tPt', 20, shape_dir, qcd_fraction=1.0, project=[mid, 200], project_axis='Y')
        #
        #    logging.warning('shape file %s created' % shape_file.GetName()) 
        #    shape_file.Close()
