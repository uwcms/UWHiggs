#from lfvmutau plotter
import os
from sys import argv, stdout, stderr
import ROOT
import sys
from FinalStateAnalysis.PlotTools.MegaBase import make_dirs
from math import sqrt

ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

shape_norm = True
if shape_norm == False:
	ynormlabel = "Normalized to Data "
else:
	ynormlabel = "Normalized to 1 "

#jobid = os.environ['jobid']
jobid = 'MCntuples_14April'
import inspect

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

lfvfilelist = []
legendlist = []    

if jobid!='MCntuples_14April': 
        lfvfilelist = ['ggHiggsToETau','vbfHiggsToETau',
                       ['Zjets_M50'],
                       ['WplusJets_madgraph'],
                       ['WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola', 'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola', 'WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola', 'ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola','ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola'],
                       ['TTJetsFullLepMGDecays', 'TTJetsSemiLepMGDecays'],
                       ['T_t-channel_TuneZ2star_8TeV-powheg-tauola','T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola','Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola','Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'], ['GluGluToHToTauTau_M-125_8TeV-powheg-pythia6', 'VBF_HToTauTau_M-125_8TeV-powheg-pythia6']
               ]
        
        legendlist = ['LFV GG Higgs', 'LFV VBF Higgs', 'DY + jets', 'W + jets', 'EWK Dibosons', 'ttbar', 'single Top', 'SM Higgs']
        
else:
        lfvfilelist = ['ggHiggsToETau','vbfHiggsToETau',
                       ['Zjets_M50'],
                       #['Zjets_M50_skimmedLL', 'Z1jets_M50_skimmedLL', 'Z2jets_M50_skimmedLL', 'Z3jets_M50_skimmedLL'],
                       #['Zjets_M50_skimmedTT', 'Z1jets_M50_skimmedTT', 'Z2jets_M50_skimmedTT', 'Z3jets_M50_skimmedTT','Z4jets_M50_skimmedTT'],
                       ['WplusJets_madgraph_skimmed','Wplus4Jets_madgraph','Wplus2Jets_madgraph_tapas', 'Wplus2Jets_madgraph', 'Wplus1Jets_madgraph_tapas', 'Wplus3Jets_madgraph', 'Wplus1Jets_madgraph'],
                       ['WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola', 'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola', 'WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola', 'ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola','ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola'],
                       #['TTJetsFullLepMGDecays', 'TTJetsSemiLepMGDecays'],
                       ['TTJetsSemiLepMGDecays'],
                       ['T_t-channel_TuneZ2star_8TeV-powheg-tauola','T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola','Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola','Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'], ['GluGluToHToTauTau_M-125_8TeV-powheg-pythia6', 'VBF_HToTauTau_M-125_8TeV-powheg-pythia6']
               ]
        
        legendlist = ['LFV GG Higgs', 'LFV VBF Higgs', 'Z/#gamma* #rightarrow ll + jets', #'Z/#gamma* #rightarrow #tau#tau + jets',
                      'W + jets', 'EWK Dibosons', 'ttbar', 'single Top', 'SM Higgs']



sourcepath = 'results/'+jobid+'/LFVHETauAnalyzer/'
sign = ['os', 'ss']
process = ['gg']#['gg','vbf']
ptcut = [0]#[0, 40]
jetN = [0]#[0,1,2,3]
legend = ROOT.TLegend(0.58,0.15,0.82,0.35)
legend.SetFillColor(0)

                        
histolist = ['tPt', 'tPhi', 'tEta', 
'ePt', 'ePhi', 'eEta', 
'et_DeltaPhi', 'et_DeltaR', 
'tPFMET_DeltaPhi', 'tPFMET_Mt', 'tMVAMET_DeltaPhi', 'tMVAMET_Mt', 
'ePFMET_DeltaPhi','ePFMET_Mt', 'eMVAMET_DeltaPhi', 'eMVAMET_Mt', 
'jetN_20', 'jetN_30',
'h_collmass_pfmet',  'h_collmass_mvamet', 'h_vismass'
]
rebins = [5,5,2,5,5,2,1, 1, 2, 5,  2, 5, 2, 5, 2, 5,1,1,1,1,1]
axistitle = ['#tau p_{T} (GeV)','#tau #phi','#tau #eta', 
'e p_{T} (GeV)','e #phi','e #eta',
'e-#tau #Delta#phi','e-#tau #DeltaR',
'#tau-PFMET #Delta#phi','#tau-PFMET M_{T} (GeV) ','#tau-MVAMET #Delta#phi','#tau-MVAMET M_{T} (GeV)',
'e-PFMET #Delta#phi','e-PFMET M_{T} (GeV)','e-MVAMET #Delta#phi','e-MVAMET #M_{T} (GeV)',
'Number of jets (p_{T}>20)','Number of jets (p_{T}>30)', 
'M_{e#tau}coll (GeV)','M_{e#tau}coll (GeV)','M_{e#tau} vis (GeV)'
]

scalefactors=[]
for myfile in lfvfilelist:
        if type(myfile) is list:
                for subfile in myfile :
                        lumifile = 'inputs/'+jobid+'/'+subfile[:(len(subfile))]+'.lumicalc.sum'
                        f = open(lumifile, 'r')
                        lumistring = f.readline()
                        #print lumistring
                        intlumi = float(lumistring)
                        f.close()
                        scalefactor = 19800./intlumi
                        scalefactors.append(scalefactor)
        else :
                lumifile = 'inputs/'+jobid+'/'+myfile[:(len(myfile))]+'.lumicalc.sum'
                f = open(lumifile, 'r')
                lumistring = f.readline()
                #print lumistring
                intlumi = float(lumistring)
                f.close()
                scalefactor = 19800./intlumi
                scalefactors.append(scalefactor)
canvas = ROOT.TCanvas("canvas","canvas",800,800)
canvas.Draw()

for s in sign :
        outfile = ROOT.TFile("outfile"+s+".root", "RECREATE")
        #print lineno()
        for pr in process:
                #print lineno()
                for c in ptcut:
                        for jn in jetN:
                                #print lineno()
                                filepath = 'plots/'+jobid+'/LFVHETauAnalyzer/et/'+s+'/'+pr+'/ept'+str(c)+'/'+str(jn)
                        
                                startingdir = os.getcwd()
                                dirs = filepath.split('/')
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

                                
                                for nhisto, i  in enumerate(histolist):
                                        #outfile = ROOT.TFile(filepath+'/'+i+".root","RECREATE")
                                        outfile.cd()
                                        legend.Clear()
                                        LFVStack = ROOT.THStack("stack"+i,"")
                                        histos = [] #fare vector di histo e clonare quello che leggo 
                                        ETaufiles = []
                                        nfile=-1
                                        newhistos =[]
                                        for n, file in enumerate(lfvfilelist):
                                                nfile+=1
                                                if type(file) is not list and file.endswith('HiggsToETau') : continue
                                                #print n, " " , file 
                                                nn = n-2
                                                if type(file) is list :
                                                
                                                        file0=ROOT.TFile(sourcepath+file[0]+'.root')
                                                        h0 = file0.Get(s+'/'+pr+'/ept'+str(c)+'/'+str(jn)+'/'+i).Clone()
                                                        h0.Scale(scalefactors[nfile])
                                                        newhisto=h0.Clone()
                                                        outfile.cd()
                                                        newhistos.append(file0.Get(s+'/'+pr+'/ept'+str(c)+'/'+str(jn)+'/'+i).Clone())
                                                        newhistos[nn].Scale(scalefactors[nfile])
                                                        newhistos[nn].SetName('h'+str(n))
                                                        for isub, subfile in enumerate(file) :
                                                                #print sourcepath+subfile+'.root'
                                                                if isub ==0 : continue
                                                                nfile +=1
                                                                newFile = ROOT.TFile(sourcepath+subfile+'.root')
                                                                outfile.cd()
                                                                newhistos[nn].Add(newFile.Get(s+'/'+pr+'/ept'+str(c)+'/'+str(jn)+'/'+i).Clone(),scalefactors[nfile])
                                                                
                                                                #print lineno() , nfile
                                                                newFile.Close()
                                                
                                                        histos.append(newhistos[nn])
                                                        histos[nn].Rebin(rebins[nhisto])
                                                        histos[nn].SetLineColor(n+1)
                                                        histos[nn].SetFillColor(n+1)
                                                        histos[nn].SetMarkerColor(n+1)
                                                        # LFVStack.Add(histos[nn])
                                                        #legend.AddEntry(histos[nn], legendlist[n])
                                                        file0.Close()
                                                
                                                #else:
                                               
                                                        #ETaufiles.append(ROOT.TFile(sourcepath+file+'.root'))
                                                        #histos.append(ETaufiles[nn].Get(s+'/'+pr+'/ept'+str(c)+'/'+str(jn)+'/'+i).Clone())
                                                        #histos[nn].SetName('h'+str(n))
                                                        #histos[nn].Scale(scalefactors[nfile])
                                                        #histos[nn].Rebin(rebins[nhisto])
                                                        #histos[nn].SetLineWidth(1)
                                                        #histos[nn].SetLineColor(n+1)
                                                        #histos[nn].SetFillColor(n+1)
                                                        #histos[nn].SetMarkerColor(n+1)
                                                
                                                #LFVStack.Add(histos[nn])
                                                #legend.AddEntry(histos[nn], legendlist[n])
                                        
                                        legend.Clear()
                                        LFVStack.Add(histos[4])
                                        LFVStack.Add(histos[5])
                                        LFVStack.Add(histos[2])
                                        LFVStack.Add(histos[3])
                                        LFVStack.Add(histos[1])
                                        LFVStack.Add(histos[0])
                                        legend.AddEntry(histos[4],legendlist[6])
                                        legend.AddEntry(histos[5],legendlist[7])
                                        legend.AddEntry(histos[2],legendlist[4])
                                        legend.AddEntry(histos[3],legendlist[5])
                                        legend.AddEntry(histos[0],legendlist[2])
                                        legend.AddEntry(histos[1],legendlist[3])
                                        # print histos
                                        canvas.cd()
                                        canvas.Clear()
                                        #LFVStack.Print()
                                        #newStack=LFVStack.Clone()
                                        #newStack.RecursiveRemove()
                                        #newStacl.Add(LFVStack.GetHists() )
                                
                                        LFVStack.Draw('hist')
                                        
                                        LFVStack.GetXaxis().SetTitle(axistitle[nhisto])

                                        ggLFVHfile = ROOT.TFile(sourcepath+lfvfilelist[0]+'.root')
                                        vbfLFVHfile = ROOT.TFile(sourcepath+lfvfilelist[1]+'.root')
                                        ggLFVHhisto= ggLFVHfile.Get(s+'/'+pr+'/ept'+str(c)+'/'+str(jn)+'/'+i).Clone()
                                        vbfLFVHhisto= vbfLFVHfile.Get(s+'/'+pr+'/ept'+str(c)+'/'+str(jn)+'/'+i).Clone()
                                        ggLFVHhisto.Rebin(rebins[nhisto])
                                        vbfLFVHhisto.Rebin(rebins[nhisto])
 
                                        ggLFVHhisto.SetLineColor(1)
                                        ggLFVHhisto.SetLineStyle(1)
                                        ggLFVHhisto.SetLineWidth(2)
                                        vbfLFVHhisto.SetLineColor(1)
                                        vbfLFVHhisto.SetLineStyle(2)
                                        vbfLFVHhisto.SetLineWidth(2)
                                        outfile.cd()
                                        ibin=0
                                        snhisto = ggLFVHhisto.Clone()
                                        
                                        snhisto.SetName("significance_vs_"+i)
                                        snhisto.SetTitle("significance_vs_"+i)
                                        
                                        while ibin<=LFVStack.GetXaxis().GetNbins():
                                                ibin+=1
                                                bkg = 0
                                                for h in histos:
                                                        bkg+=h.GetBinContent(ibin)
                                                sig = ggLFVHhisto.GetBinContent(ibin) + vbfLFVHhisto.GetBinContent(ibin)
                                                
                                                if sig+bkg> 0:
                                                        sn = sig/sqrt(sig+bkg)
                                                        snhisto.SetBinContent(ibin, sn)
                                                        err = sqrt(sig/(2*sqrt(sig+bkg)))
                                                        snhisto.SetBinError(ibin, err)     
                                                        #print bkg, sig, sn, err
                                                snhisto.GetXaxis().SetTitle(axistitle[nhisto])
                                                snhisto.GetYaxis().SetTitle("S/#sqrt{S+B}")

                                        ggLFVHhisto.Scale(10*scalefactors[0])
                                        vbfLFVHhisto.Scale(10*scalefactors[1])
                                        mymax = max (ggLFVHhisto.GetBinContent(ggLFVHhisto.GetMaximumBin()), vbfLFVHhisto.GetBinContent(vbfLFVHhisto.GetMaximumBin()))
                                        if mymax > LFVStack.GetMaximum() :  LFVStack.SetMaximum( 1.2*mymax )

                                        ggLFVHhisto.Draw('SAMEHIST')
                                        vbfLFVHhisto.Draw('SAMEHIST')
                                        legend.AddEntry(ggLFVHhisto, 'LFV GG Higgs x 10')
                                        legend.AddEntry(vbfLFVHhisto, 'LFV VBF Higgs x 10')                                
                                        legend.Draw() 
                                        canvas.SetLogy(0)
                                        if i.endswith('Pt') : canvas.SetLogy(1)
                                        canvas.Update()
                                        canvas.SaveAs(filepath+'/mc_'+i+'.png')
                                        #outfile.cd()
                                        
                                        #canvas.Clear()
                                        snhisto.Draw()
                                        canvas.SaveAs(filepath+'/mc_significance_vs_'+i+'.png')
                                        #canvas.Write()

                                                                                
                                
        outfile.Close()
        os.system("rm outfile"+s+".root")
 
