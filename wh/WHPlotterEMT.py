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
from WHPlotterBase import make_styler

logging.basicConfig(stream=sys.stderr, level=logging.INFO)

class WHPlotterEMT(WHPlotterBase.WHPlotterBase):
    def __init__(self, files, lumifiles, outputdir):
        super(WHPlotterEMT, self).__init__(files, lumifiles, outputdir)

if __name__ == "__main__":
    jobid = os.environ['jobid']

    print "Plotting EMT for %s" % jobid

    # Figure out if we are 7 or 8 TeV
    period = '7TeV' if '7TeV' in jobid else '8TeV'

    sqrts = 7 if '7TeV' in jobid else 8

    samples = [
        'Zjets_M50',
        'WplusJets_madgraph',
        'WZJetsTo3LNu*',
        'ZZ*',
        'WW*',
        'VH*',
        'WH*',
        'TTplusJets_madgraph',
        "data_MuEG*",
    ]

    files = []
    lumifiles = []

    for x in samples:
        files.extend(glob.glob('results/%s/WHAnalyzeEMT/%s.root' % (jobid, x)))
        lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))

    outputdir = 'results/%s/plots/emt' % jobid
    plotter = WHPlotterEMT(files, lumifiles, outputdir)

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

    ###########################################################################
    ##  Making shape file     #################################################
    ###########################################################################

    shape_file = ROOT.TFile(
        os.path.join(outputdir, 'emt_shapes_%s.root' % period), 'RECREATE')
    shape_dir = shape_file.mkdir('emt')
    plotter.write_shapes('subMass', 20, shape_dir, unblinded=True, qcd_fraction=0.5)
    shape_dir = shape_file.mkdir('emt_w')
    plotter.write_shapes('subMass', 20, shape_dir, unblinded=True, qcd_fraction=0.0)
    shape_dir = shape_file.mkdir('emt_q')
    plotter.write_shapes('subMass', 20, shape_dir, unblinded=True, qcd_fraction=1.0)
    #plotter.write_cut_and_count('subMass', shape_dir, unblinded=True)
    shape_file.Close()
