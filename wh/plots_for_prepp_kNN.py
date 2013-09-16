from FinalStateAnalysis.Utilities.rootbindings import ROOT
import os
from rootpy.utils import asrootpy
from FinalStateAnalysis.PlotTools.RebinView  import RebinView
import rootpy.plotting as plotting
import logging
import sys

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

def round_to_ints(histo):
    new = histo.Clone()
    new.Reset()
    for bin in range(histo.GetNbinsX()+1):
        nentries = ROOT.TMath.Nint(histo.GetBinContent(bin)) \
                   if histo.GetBinContent(bin) >= 0 else 0
        centerx  = histo.GetXaxis().GetBinCenter(bin)
        for _ in range(nentries):
            new.Fill(centerx)
    return new


ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

jobid = os.environ['jobid']

rebinPt = (
    [10,12,15,20,25,30,35,40,45,50,60,70,100], #,150,200],
    [10,12,15,20,25,30,35,40,45,50,60,70,100], #range(0,50,5)+range(50,110,10),
    )

rebinJet = (
    [0,1,2,3,4,6,9,12],
    [0,1,2,3,4,6,9,12],
    #range(13),
    )

axes = {
    'muonPt'        : 'muon p_{T} (GeV)',
    'muonJetPt'     : 'muon p_{T}^{Jet} (GeV)',
    'numJets20'     : '# Jets',
    'electronPt'    : 'electron p_{T} (GeV)', 
    'electronJetPt' : 'electron p_{T}^{Jet} (GeV)',
}

sources = {
    ##MMT
    'm_mmt_subleading_kNN' : {
        'vars'  : {'muonPt':rebinPt, 'muonJetPt':rebinPt, 'numJets20':rebinJet},
        'wjets' : 'results/%s/fakerate_fits/m_wjets_pt10_h2taucuts020_muonInfo_k100.kNN.root' % jobid,
        'qcd'   : 'results/%s/fakerate_fits/m_qcd_pt10_h2taucuts020_muonInfo_k100.kNN.root'   % jobid,
        },
    'm_mmt_leading_kNN' : {
        'vars'  : {'muonPt':rebinPt, 'muonJetPt':rebinPt, 'numJets20':rebinJet},
        'wjets' : 'results/%s/fakerate_fits/m_wjets_pt10_h2taucuts_muonInfo_k100.kNN.root' % jobid,
        'qcd'   : 'results/%s/fakerate_fits/m_qcd_pt10_h2taucuts_muonInfo_k100.kNN.root'   % jobid,
        },

    #EMT
    'm_emt_kNN' : {
        'vars'  : {'muonPt':rebinPt, 'muonJetPt':rebinPt, 'numJets20':rebinJet},
        'wjets' : 'results/%s/fakerate_fits/m_Mwjets_pt10_h2taucuts_muonInfo_k100.kNN.root' % jobid,
        'qcd'   : 'results/%s/fakerate_fits/m_Mqcd_pt10_h2taucuts_muonInfo_k100.kNN.root'   % jobid,
        },
    'e_emt_kNN' : {
        'vars'  : {'electronPt':rebinPt, 'electronJetPt':rebinPt, 'numJets20':rebinJet},
        'wjets' : 'results/%s/fakerate_fits/e_wjets_pt10_eid12Medium_h2taucuts_electronInfo_k100.kNN.root' % jobid,
        'qcd'   : 'results/%s/fakerate_fits/e_qcd_pt10_eid12Medium_h2taucuts_electronInfo_k100.kNN.root'   % jobid,
        },

    #EET
    'e_eet_subleading_kNN' : {
        'vars'  : {'electronPt':rebinPt, 'electronJetPt':rebinPt, 'numJets20':rebinJet},
        'wjets' : 'results/%s/fakerate_fits/ee_wjetsNoZmass_pt10_eid12Medium_h2taucuts020_electronInfo_k100.kNN.root' % jobid,
        'qcd'   : 'results/%s/fakerate_fits/ee_qcd_pt10_eid12Medium_h2taucuts020_electronInfo_k100.kNN.root'   % jobid,
        },
    'e_eet_leading_kNN' : {
        'vars'  : {'electronPt':rebinPt, 'electronJetPt':rebinPt, 'numJets20':rebinJet},
        'wjets' : 'results/%s/fakerate_fits/ee_wjetsNoZmass_pt10_eid12Tight_h2taucuts_electronInfo_k100.kNN.root' % jobid,
        'qcd'   : 'results/%s/fakerate_fits/ee_qcd_pt10_eid12Tight_h2taucuts_electronInfo_k100.kNN.root'   % jobid,
        },
}

canvas = plotting.Canvas(800, 800, name='adsf', title='asdf')
canvas.SetLogy(True)

#from pdb import set_trace; set_trace()
for output, info in sources.iteritems():
    logging.info("analyzing %s" % output)
    wjets_name = info['wjets']
    qcd_name   = info['qcd']
    variables  = info['vars']

    wjets_file = ROOT.TFile.Open(wjets_name)
    logging.debug('filename: %s instance: %s' % (wjets_name, wjets_file.__repr__()))
    qcd_file   = ROOT.TFile.Open(qcd_name)
    logging.debug('filename: %s instance: %s' % (qcd_name, qcd_file.__repr__()))

    for var, rebins in variables.iteritems():
        logging.debug('variable %s' % var)
        #####
        #WJets fake rate
        #####
        wjets_view = RebinView(wjets_file, rebins[0])
        wjets_eff = asrootpy( 
            ROOT.TGraphAsymmErrors( 
                round_to_ints( wjets_view.Get('%s_pass' % var)), 
                round_to_ints( wjets_view.Get('%s_all' % var))
            ) 
        )
        wjets_eff.SetTitle('WJets fake rate')
        wjets_eff.markerstyle = 20
        wjets_eff.markercolor = ROOT.EColor.kBlue

        #####
        #QCD fake rate
        #####
        qcd_view = RebinView(qcd_file, rebins[0])
        qcd_eff = asrootpy( 
            ROOT.TGraphAsymmErrors( 
                round_to_ints(qcd_view.Get('%s_pass' % var)), 
                round_to_ints(qcd_view.Get('%s_all' % var))
            ) 
        )
        qcd_eff.SetTitle('QCD fake rate')
        qcd_eff.markerstyle = 20
        qcd_eff.markercolor = ROOT.EColor.kRed


        #####
        #WJets kNN output
        #####
        wjets_view     = RebinView(wjets_file, rebins[1])
        wjets_estimate = asrootpy( 
            wjets_view.Get('%s_estimate' % var)
        )
        wjets_estimate_denom = asrootpy( 
            wjets_view.Get('%s_estimate_all' % var)
        )
        wjets_estimate.Divide(wjets_estimate_denom)
        wjets_estimate.SetTitle('WJets kNN output')
        wjets_estimate.linecolor = ROOT.kBlue 
        wjets_estimate.linewidth = 2
        wjets_estimate.legendstyle = 'l'
        wjets_estimate.fillstyle = 0
        wjets_estimate.drawstyle = 'hist'

        #####
        #QCD kNN output
        #####
        qcd_view = RebinView(qcd_file, rebins[1])
        qcd_estimate = asrootpy( 
            qcd_view.Get('%s_estimate' % var)
        )
        qcd_estimate_denom = asrootpy( 
            qcd_view.Get('%s_estimate_all' % var)
        )
        qcd_estimate.Divide(qcd_estimate_denom)
        qcd_estimate.SetTitle('QCD kNN output')
        qcd_estimate.linecolor = ROOT.kRed
        qcd_estimate.linewidth = 2
        qcd_estimate.fillstyle = 0
        qcd_estimate.legendstyle = 'l'
        qcd_estimate.drawstyle = 'hist'
        
        canvas.SetGrid()

        wjets_estimate.GetXaxis().SetTitle(axes[var])
        wjets_estimate.GetYaxis().SetTitle('fake rate')
        wjets_estimate.GetYaxis().SetRangeUser(10**-2,1.)
        for axis in [wjets_estimate.GetXaxis(), wjets_estimate.GetYaxis()]:
            axis.SetLabelSize( axis.GetLabelSize()*0.6 )
            axis.SetTitleSize( axis.GetTitleSize()*0.6 )
            #axis.SetTitleOffset( axis.GetTitleOffset()*0.6 )

        wjets_estimate.Draw()
        qcd_estimate.Draw('same')

        wjets_eff.Draw('P SAME')
        qcd_eff.Draw('P SAME')

        #####
        #Legend
        #####
        legend = plotting.Legend(4, rightmargin=0.02, topmargin=0.05, leftmargin=0.45)
        legend.AddEntry(wjets_estimate)
        legend.AddEntry(qcd_estimate)
        legend.AddEntry(wjets_eff)
        legend.AddEntry(qcd_eff)
        legend.SetTextSize(0.5*legend.GetTextSize())
        legend.SetEntrySeparation(0.0)
        legend.SetMargin(0.35)
        legend.Draw()
        
        canvas.Update()
        canvas.SetLogy(True)
        print "saving results/%s/fakerate_fits/%s_%s.pdf(png)" % (jobid, output, var)
        canvas.SaveAs('results/%s/fakerate_fits/%s_%s.pdf' % (jobid, output, var))
        canvas.SaveAs('results/%s/fakerate_fits/%s_%s.png' % (jobid, output, var))


    wjets_file.Close()
    qcd_file.Close()  
