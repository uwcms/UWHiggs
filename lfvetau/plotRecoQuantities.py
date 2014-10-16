#from mauro plotters
import os
from sys import argv, stdout, stderr
import ROOT
import sys
from FinalStateAnalysis.PlotTools.MegaBase import make_dirs
import glob
import logging
import sys

#jobid = os.environ['jobid']
jobid = 'MCntuples_3March' 
channel = 'et'

print "\nPlotting %s for %s\n" % (channel, jobid)

mcsamples = [
        'ggHiggsToETau',
        'vbfHiggsToETau',
        'GluGluToHToTauTau_M-125_8TeV-powheg-pythia6',
        'VBF_HToTauTau_M-125_8TeV-powheg-pythia6',
        #                'TTJets*',
        #                'T_t*',
        #                'Tbar_t', 
        #                'Wplus*',
        #                'WWJets*',
        #                'WZJets*',
        #                'ZZJets*',
        #                'Z*jets_M50'
]
        
files = []
lumifiles = []

for x in mcsamples:
        files.extend(glob.glob('results/%s/LFVHETauAnalyzer/%s.root' % (jobid, x)))
        lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))
        outputdir = 'plots/%s/LFVHETauAnalyzer/%s' % (jobid, channel)
        base_out_dir = outputdir
        if not os.path.exists(outputdir):
                os.makedirs(outputdir)
                
#       blinder = None
#        blind   = 'blind' not in os.environ or os.environ['blind'] == 'YES'
#        print '\n\nRunning Blind: %s\n\n' % blind


#if blind:
        # Don't look at the SS all pass region
#        blinder = lambda x: BlindView(x, "ss/tau_os/p1p2p3/.*")
        
        
        
import rootpy.plotting.views as views
        
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

canvas = ROOT.TCanvas("canvas","canvas",800,800)
LFVStack = ROOT.THStack("stack","")

file0 = ROOT.TFile('results/%s/LFVHETauAnalyzer/ggHiggsToETau.root'  % (jobid) )


dir = file0.Get('os/gg/ept0')
hlist = dir.GetListOfKeys()
iter = ROOT.TIter(hlist)
for i in iter:
        print i.GetName()
        histo0 = file0.Get('%s/%s' % ('os/gg/ept0', i.GetName()))
        if histo0.Integral() != 0 :
                histo0.Scale(1./histo0.Integral())
        histo0.Draw("E")
        histo0.SetLineWidth(2)
        histo0.SetLineColor(1)
        variable=''
        if histo0.GetName()[4:7] == "Phi" : variable = "#phi"
        if histo0.GetName()[4:7] == "Eta" : variable = "#eta"
        if histo0.GetName()[4:10] == "Energy" : 
                variable = "energy (GeV)"
                canvas.SetLogy(1)
        if histo0.GetName()[4:6] == "Pt" : 
                variable = "p_{T} (GeV)"
                canvas.SetLogy(1)
        if histo0.GetName()[9:17] == "DeltaPhi" : 
                variable = "#Delta#phi"
                p = 'e#tau'
                canvas.SetLogy(1)
#        axislabel = p+" "+variable
        axislabel = variable
        if variable != '' : histo0.GetXaxis().SetTitle(axislabel)
        maxy = histo0.GetBinContent(histo0.GetMaximumBin())
        legend = ROOT.TLegend(0.7,0.85,0.88,0.75)
        legend.SetFillColor(0)
        legend.AddEntry(histo0, "ggHigsToETau")

        for n, sample in enumerate(mcsamples) :
                if n == 0 : continue 
                file = ROOT.TFile("results/%s/LFVHETauAnalyzer/%s.root"  %(jobid, sample))
                histo=file.Get('%s/%s' % ('os/gg/ept0', i.GetName()))

                
                if histo.Integral() != 0  : 
                        histo.Scale(1./histo.Integral())
                        
                p = histo.GetName()[0:1]
                if p == 't' : p = '#tau'
                histo.Draw("ESAME") 
                
                histo.SetLineColor(n)
                histo.SetLineWidth(2)

                canvas.SetLogy(0)
                if maxy < histo.GetBinContent(histo.GetMaximumBin()) : maxy = histo.GetBinContent(histo.GetMaximumBin())
                legend.AddEntry(histo, "mcsamples")
        
        if canvas.GetLogy()==0 :
                histo0.GetYaxis().SetRangeUser(0, maxy*1.3)
                
        canvas.Update()
        canvas.SaveAs(outputdir+'/os_ept0_'+i.GetName()+'.png')
 
