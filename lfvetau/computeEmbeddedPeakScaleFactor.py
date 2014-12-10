from FinalStateAnalysis.Utilities.rootbindings import ROOT
import math
import logging
import sys

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch()

file_dataA = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/ZetauEmbedded_Run2012A.root')
file_dataB = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/ZetauEmbedded_Run2012B.root')
file_dataC = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/ZetauEmbedded_Run2012C.root')
file_dataD = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/ZetauEmbedded_Run2012D.root')

dataA=file_dataA.Get('os/gg/ept30/h_collmass_pfmet')
dataB=file_dataB.Get('os/gg/ept30/h_collmass_pfmet')
dataC=file_dataC.Get('os/gg/ept30/h_collmass_pfmet')
dataD=file_dataD.Get('os/gg/ept30/h_collmass_pfmet')

data=dataC.Clone()
data.Add(dataB)
data.Add(dataA)
data.Add(dataD)

c= ROOT.TCanvas("c","c", 800, 1000)
c.Draw()
c.SetGridx(1)
c.SetGridy(1)


njets=[0,1,2,3,4]

file_MC0 = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/Z0jets_M50_skimmedTT.root')
file_MC1 = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/Z1jets_M50_skimmedTT.root')
file_MC2 = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/Z2jets_M50_skimmedTT.root')
file_MC3 = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/Z3jets_M50_skimmedTT.root')
file_MC4 = ROOT.TFile('results/newNtuple_5Nov/LFVHETauAnalyzerMVA/Z4jets_M50_skimmedTT.root')

mc0 = file_MC0.Get('os/gg/ept30/h_collmass_pfmet')
mc1 = file_MC1.Get('os/gg/ept30/h_collmass_pfmet')
mc2 = file_MC2.Get('os/gg/ept30/h_collmass_pfmet')
mc3 = file_MC3.Get('os/gg/ept30/h_collmass_pfmet')
mc4 = file_MC4.Get('os/gg/ept30/h_collmass_pfmet')

datalumi = 0
#for i in ['A', 'B', 'C', 'D'] :
for i in ['C', 'D'] :
    lumifile = 'inputs/newNtuple_5Nov/data_SingleElectron_Run2012%s_22Jan2013_v1.lumicalc.sum' %(i)
    f = open(lumifile, 'r')
    lumistring = f.readline()
    datalumi += float(lumistring)


mc_lumi=[]

for i in range(0,5):
    
    lumifile = 'inputs/newNtuple_5Nov/Z%sjets_M50_skimmedTT.lumicalc.sum' %(str(int(i)))
    f = open(lumifile, 'r')
    lumistring = f.readline()
    mc_lumi.append(float(lumistring))

mc=mc0.Clone()
mc.Scale(datalumi/mc_lumi[0])
mc.Add(mc1, datalumi/mc_lumi[1])
mc.Add(mc2, datalumi/mc_lumi[2])
mc.Add(mc3, datalumi/mc_lumi[3])
mc.Add(mc4, datalumi/mc_lumi[4])
mc.SetMarkerStyle(20)
mc.SetMarkerColor(4)
mc.SetLineColor(4)

data.Scale(mc.Integral()/data.Integral()  )

max_data = data.GetBinContent(data.GetMaximumBin())
max_mc = mc.GetBinContent(mc.GetMaximumBin())
mc.Draw()
data.SetMarkerStyle(20)
data.SetMarkerColor(2)
data.SetLineColor(2)
data.Draw()

m          = ROOT.RooRealVar('m', 'Collinear Mass', 90,40,250)
os_data    = ROOT.RooDataHist('z_data','z_data',ROOT.RooArgList(m),data)
os_mc      = ROOT.RooDataHist('z_mc','z_mc',ROOT.RooArgList(m),mc)
mean       = ROOT.RooRealVar('mean', 'mean', 90,70,110)
sigmaL     = ROOT.RooRealVar('sigmaL' ,'sigmaL' ,10 ,0 ,100)
sigmaR     = ROOT.RooRealVar('sigmaR' ,'sigmaR' ,10 ,0 ,100)
alphaL     = ROOT.RooRealVar('alphaL' ,'alphaL' ,1  ,0 ,30 )
alphaR     = ROOT.RooRealVar('alphaR' ,'alphaR' ,1  ,0 ,30 )
command    = ROOT.RooCruijff("os_func", "os_func",m, mean, sigmaL, sigmaR, alphaL, alphaR)


frame = m.frame(ROOT.RooFit.Title("Z mass peak"))

fit_result_mc = command.fitTo(
    os_mc,
    ROOT.RooFit.Save(True))

os_mc.plotOn(
    frame,ROOT.RooFit.LineColor(ROOT.EColor.kBlue),ROOT.RooFit.MarkerColor(ROOT.EColor.kBlue)
)
mcpeak = mean.getVal()
mcpeakerr = mean.getError()

command.plotOn(frame,#ROOT.RooFit.VisualizeError(fit_result_mc,1),
               ROOT.RooFit.LineColor(ROOT.EColor.kBlue))

leg = ROOT.TLegend(0.5,0.8,0.85,0.7)


print '____________________'

fit_result = command.fitTo(
    os_data,
    ROOT.RooFit.Save(True))
os_data.plotOn(frame,ROOT.RooFit.LineColor(ROOT.EColor.kRed),ROOT.RooFit.MarkerColor(ROOT.EColor.kRed) )
command.plotOn(frame,ROOT.RooFit.LineColor(ROOT.EColor.kRed))
datapeak = mean.getVal()
datapeakerr= mean.getError()
print 'datapeak : %f +/- %f, mcpeak: %f +/- %f' %(datapeak,datapeakerr, mcpeak, mcpeakerr )

leg.AddEntry(mc, "DY->#tau#tau (MC)", "lp")
leg.AddEntry(data, "DY->#tau#tau (embedded)", "lp")

pt = ROOT.TPaveText(0.5, 0.5, 0.85, 0.65)
pt.AddText("MC peak: %.2f #pm %.2f GeV" %(mcpeak, mcpeakerr))
pt.AddText("Emb peak: %.2f #pm %.2f GeV" %(datapeak, datapeakerr))

frame.addObject(pt)

#text = ROOT.TText(124.5,1400, 'MC peak: %.2f #pm %.2f GeV' %(mcpeak, mcpeakerr))
#text2 = ROOT.TText(120, 1300, 'Emb peak: %.2f #pm %.2f GeV'%(datapeak, datapeakerr))
#text.SetTextSize(0.03)
#text.SetTextColor(4)
#text2.SetTextSize(0.03)
#frame.addObject(text)
#text2.SetTextColor(2)
#frame.addObject(text2)
leg.SetFillColor(0)
frame.addObject(leg)
frame.Draw()
#pt.AddText(text)
#pt.AddText(text2)
#Leg.Draw()
#leg.SetFillColor(0)
pt.SetFillColor(0)
pt.SetShadowColor(0)
pt.Paint()

c.Update()
c.SaveAs('Zfit.pdf')
c.SaveAs('Zfit.png')
