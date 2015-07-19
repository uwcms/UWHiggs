#from lfvmutau plotter
import os
from sys import argv, stdout, stderr
import ROOT
import sys
from FinalStateAnalysis.PlotTools.MegaBase import make_dirs

ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

shape_norm = True
if shape_norm == False:
	ynormlabel = "Normalized to Data "
else:
	ynormlabel = "Normalized to 1 "

canvas = ROOT.TCanvas("canvas","canvas",800,800)
LFVStack = ROOT.THStack("stack","")

#lfvfilelist = ['results/MCntuples25Feb/LFVHAnalyzeGEN/ggHiggsToETau.root', 'results/MCntuples25Feb/LFVHAnalyzeGEN/vbfHiggsToETau.root']
#smfilelist = ['results/MCntuples25Feb/LFVHAnalyzeGEN/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root', 'results/MCntuples25Feb/LFVHAnalyzeGEN/VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root']

lfvfilelist = ['results/MCntuples_otherCh2/LFVHAnalyzeGENEMu/ggHiggsToMuTau.root', 'results/MCntuples_otherCh2/LFVHAnalyzeGENEMu/vbfHiggsToMuTau.root','results/MCntuples_otherCh2/LFVHAnalyzeGENMuTau/ggHiggsToMuTau.root', 'results/MCntuples_otherCh2/LFVHAnalyzeGENMuTau/vbfHiggsToMuTau.root']
smfilelist = ['results/MCntuples_otherCh2/LFVHAnalyzeGENEMu/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root', 'results/MCntuples_otherCh2/LFVHAnalyzeGENEMu/VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root','results/MCntuples_otherCh2/LFVHAnalyzeGENMuTau/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6.root', 'results/MCntuples_otherCh2/LFVHAnalyzeGENMuTau/VBF_HToTauTau_M-125_8TeV-powheg-pythia6.root']

for n, file in enumerate(lfvfilelist):

        ETaufile = ROOT.TFile(file)
        SMHTauTaufile = ROOT.TFile(smfilelist[n])

        gendir = ETaufile.Get('gen')
        hlist = gendir.GetListOfKeys()

        iter = ROOT.TIter(hlist)

        filepath = 'plots/'+ETaufile.GetName()[8 : len(ETaufile.GetName()) -5]
        print filepath
        startingdir = os.getcwd()
        print startingdir
        dirs = filepath.split('/')
        print dirs
        thechannel = 'e #tau_{h}'
        for d in dirs:
                currentdir = os.getcwd()
                dirlist = os.listdir(currentdir)
                
                if d in dirlist:
                        os.chdir(d)
                else: 
                        os.makedirs(d)
                        os.chdir(d)
		if d == "LFVHAnalyzeGENEMu" : 
                        thechannel='#mu #tau_{e}'
                        thesmchannel = "#tau_{#mu}#tau_{e}"
		if d == "LFVHAnalyzeGENMuTau" : 
                        thechannel='#mu #tau_{h}'
                        thesmchannel = "#tau_{#mu}#tau_{h}"
                
        os.chdir(startingdir)
                
        for i in iter:
                lfv_histo = ETaufile.Get('gen/'+i.GetName())
                sm_histo = SMHTauTaufile.Get('gen/'+i.GetName())
                
                if lfv_histo.Integral() != 0 and sm_histo.Integral() != 0  : 
                        lfv_histo.Scale(1./lfv_histo.Integral())
                        sm_histo.Scale(1./sm_histo.Integral())

                        p = sm_histo.GetName()[0:1]
                        if p == 't' : p = '#tau'
                        if p == 'h' :
                                lfv_histo.Draw("E") 
                        else:

                                sm_histo.Draw("E")
                                lfv_histo.Draw("ESAME") 
                        
                        lfv_histo.SetLineWidth(2)
                        sm_histo.SetLineColor(2)
                        sm_histo.SetLineWidth(2)
                        
                               
                        canvas.SetLogy(0)
                        variable=''
                        if sm_histo.GetName()[4:7] == "Phi" : variable = "#phi"
                        if sm_histo.GetName()[4:7] == "Eta" : variable = "#eta"
                        if sm_histo.GetName()[4:10] == "Energy" : 
                                variable = "energy (GeV)"
                                canvas.SetLogy(1)
                        if sm_histo.GetName()[4:6] == "Pt" : 
                                variable = "p_{T} (GeV)"
                                canvas.SetLogy(1)
                        if sm_histo.GetName()[9:17] == "DeltaPhi" : 
                                variable = "#Delta#phi"
                                p = 'e#tau'
                                canvas.SetLogy(1)
                        if sm_histo.GetName() == "higgsPt" : 
                                variable = "p_{T} (GeV)"
                                canvas.SetLogy(1)
                        axislabel = p+" "+variable
                        if variable != '' : lfv_histo.GetXaxis().SetTitle(axislabel)

                        maxlfv = lfv_histo.GetBinContent(lfv_histo.GetMaximumBin())
                        maxsm  = sm_histo.GetBinContent(sm_histo.GetMaximumBin())
                        if canvas.GetLogy()==0 :
                                if maxlfv > maxsm :
                                        lfv_histo.GetYaxis().SetRangeUser(0, maxlfv*1.3)
                                else :
                                        lfv_histo.GetYaxis().SetRangeUser(0, maxsm*1.3)
                                
                        if sm_histo.GetName()[4:13] == 'DecayMode': 
                                legend = ROOT.TLegend(0.6,0.75,0.8,0.85)
                        else :
                                legend = ROOT.TLegend(0.38,0.15,0.62,0.25)
                        
                        legend.SetFillColor(0)
                        if  thechannel == 'e #tau_{h}' :
                                legend.AddEntry(sm_histo, "H #rightarrow #tau_{e}#tau_{h} ")
                                legend.AddEntry(lfv_histo, "H #rightarrow e#tau_{h} ")
                        else:
                                legend.AddEntry(sm_histo, "H #rightarrow " +thesmchannel)
                                legend.AddEntry(lfv_histo, "H #rightarrow "+ thechannel)
                               
                        if p != 'h' :
                                legend.Draw()

                        canvas.Update()
                        canvas.SaveAs(filepath+'/gen_'+i.GetName()+'.png')
 
