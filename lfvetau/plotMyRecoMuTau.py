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


import inspect

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno


lfvfilelist = ['ggHiggsToMuTau', 'vbfHiggsToMuTau',
'Zjets_M50',
'WplusJets_madgraph',
#'Wplus4Jets_madgraph.root',
#'Wplus2Jets_madgraph_tapas.root',
#'Wplus2Jets_madgraph.root',
#'Wplus1Jets_madgraph_tapas.root',
#'Wplus3Jets_madgraph.root',
#'Wplus1Jets_madgraph.root',
['WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola', 'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola', 'WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola', 'ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola','ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola'],
['TTJetsFullLepMGDecays', 'TTJetsSemiLepMGDecays'],
['T_t-channel_TuneZ2star_8TeV-powheg-tauola','T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola','Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola','Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'], ['GluGluToHToTauTau_M-125_8TeV-powheg-pythia6', 'VBF_HToTauTau_M-125_8TeV-powheg-pythia6']
]
 


legendlist = ['LFV GG Higgs', 'LFV VBF Higgs', 'DY + jets', 'W + jets', 'EWK Dibosons', 'ttbar', 'single Top', 'SM Higgs'
]
#print 'lfvfilelist ', len(lfvfilelist), ' legendlist ' , len(legendlist)
sourcepath = 'results/MCntuples_8April/LFVHMuTauAnalyzer/'
sign = ['os', 'ss']
process = ['gg','vbf']
ptcut = [0, 40]
jetN = [0,1,2,3]
legend = ROOT.TLegend(0.58,0.15,0.82,0.35)
legend.SetFillColor(0)

                        
histolist = ['tPt', 'tPhi', 'tEta', 
'mPt', 'mPhi', 'mEta', 
'mt_DeltaPhi', 'mt_DeltaR', 
'tPFMET_DeltaPhi', 'tPFMET_Mt', 'tMVAMET_DeltaPhi', 'tMVAMET_Mt', 
'mPFMET_DeltaPhi','mPFMET_Mt', 'mMVAMET_DeltaPhi', 'mMVAMET_Mt', 
'jetN_20', 'jetN_30',
'h_collmass_pfmet',  'h_collmass_mvamet', 'h_vismass'
]
rebins = [5,5,2,5,5,2,1, 1, 2, 5,  2, 5, 2, 5, 2, 5,1,1,1,1,1]
axistitle = ['#tau p_{T} (GeV)','#tau #phi','#tau #eta', 
'#mu p_{T} (GeV)','#mu #phi','#mu #eta',
'#mu-#tau #Delta#phi','#mu-#tau #DeltaR',
'#tau-PFMET #Delta#phi','#tau-PFMET M_{T} (GeV) ','#tau-MVAMET #Delta#phi','#tau-MVAMET M_{T} (GeV)',
'#mu-PFMET #Delta#phi','#mu-PFMET M_{T} (GeV)','#mu-MVAMET #Delta#phi','#mu-MVAMET #M_{T} (GeV)',
'Number of jets (p_{T}>20)','Number of jets (p_{T}>30)', 
'M_{#mu#tau}coll (GeV)','M_{#mu#tau}coll (GeV)','M_{#mu#tau} vis (GeV)'
]

scalefactors=[]
for myfile in lfvfilelist:
        if type(myfile) is list:
                for subfile in myfile :
                        lumifile = 'inputs/MCntuples_8April/'+subfile[:(len(subfile))]+'.lumicalc.sum'
                        f = open(lumifile, 'r')
                        lumistring = f.readline()
                        #print lumistring
                        intlumi = float(lumistring)
                        f.close()
                        scalefactor = 19800./intlumi
                        scalefactors.append(scalefactor)
        else :
                lumifile = 'inputs/MCntuples_8April/'+myfile[:(len(myfile))]+'.lumicalc.sum'
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
        #print lineno()
        for pr in process:
                #print lineno()
                for c in ptcut:
                        for jn in jetN:
                                #print lineno()
                                filepath = 'plots/MCntuples_8April/LFVHETauAnalyzer/mt/'+s+'/'+pr+'/mpt'+str(c)+'/'+str(jn)
                        
                                startingdir = os.getcwd()
                                dirs = filepath.split('/')
                                thechannel = '#mu #tau_{h}'
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
                                        outfile = ROOT.TFile(filepath+'/'+i+".root","RECREATE")
                                        legend.Clear()
                                        LFVStack = ROOT.THStack("stack"+i,"")
                                        histos = [] #fare vector di histo e clonare quello che leggo 
                                        ETaufiles = []
                                        nfile=-1
                                        newhistos =[]
                                        for n, file in enumerate(lfvfilelist):
                                                nfile+=1
                                                if type(file) is not list and file.endswith('HiggsToMuTau') : continue
                                                #print n, " " , file 
                                                nn = n-2
                                                if type(file) is list :
                                                
                                                        file0=ROOT.TFile(sourcepath+file[0]+'.root')
                                                        h0 = file0.Get(s+'/'+pr+'/mpt'+str(c)+'/'+str(jn)+'/'+i).Clone()
                                                        h0.Scale(scalefactors[nfile])
                                                        newhisto=h0.Clone()
                                                        outfile.cd()
                                                        newhistos.append(file0.Get(s+'/'+pr+'/mpt'+str(c)+'/'+str(jn)+'/'+i).Clone())
                                                        newhistos[n-4].Scale(scalefactors[nfile])
                                                
                                                        newhistos[n-4].SetName('h'+str(n))
                                                        for isub, subfile in enumerate(file) :
                                                                #print sourcepath+subfile+'.root'
                                                                if isub ==0 : continue
                                                                nfile +=1
                                                                newFile = ROOT.TFile(sourcepath+subfile+'.root')
                                                                outfile.cd()
                                                                newhistos[n-4].Add(newFile.Get(s+'/'+pr+'/mpt'+str(c)+'/'+str(jn)+'/'+i).Clone(),scalefactors[nfile])
                                                                #newhistos[n-4].Scale()
                                                                #print lineno() , nfile
                                                                newFile.Close()
                                                
                                                        histos.append(newhistos[n-4])
                                                        histos[nn].Rebin(rebins[nhisto])
                                                        histos[nn].SetLineColor(n+1)
                                                        histos[nn].SetFillColor(n+1)
                                                        histos[nn].SetMarkerColor(n+1)
                                                        # LFVStack.Add(histos[nn])
                                                        #legend.AddEntry(histos[nn], legendlist[n])
                                                        file0.Close()
                                                
                                                else:
                                               
                                                        ETaufiles.append(ROOT.TFile(sourcepath+file+'.root'))
                                                        
                                                        histos.append(ETaufiles[nn].Get(s+'/'+pr+'/mpt'+str(c)+'/'+str(jn)+'/'+i).Clone())
                                                        histos[nn].Scale(scalefactors[nfile])
                                                        histos[nn].SetName('h'+str(n))
                                                        
                                                        histos[nn].Rebin(rebins[nhisto])
                                                        histos[nn].SetLineWidth(1)
                                                        histos[nn].SetLineColor(n+1)
                                                        histos[nn].SetFillColor(n+1)
                                                        histos[nn].SetMarkerColor(n+1)
                                                
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
                                        ggLFVHhisto= ggLFVHfile.Get(s+'/'+pr+'/mpt'+str(c)+'/'+str(jn)+'/'+i).Clone()
                                        vbfLFVHhisto= vbfLFVHfile.Get(s+'/'+pr+'/mpt'+str(c)+'/'+str(jn)+'/'+i).Clone()
                                        ggLFVHhisto.Rebin(rebins[nhisto])
                                        vbfLFVHhisto.Rebin(rebins[nhisto])

                                        ggLFVHhisto.SetLineColor(1)
                                        ggLFVHhisto.SetLineStyle(1)
                                        ggLFVHhisto.SetLineWidth(2)
                                        vbfLFVHhisto.SetLineColor(1)
                                        vbfLFVHhisto.SetLineStyle(2)
                                        vbfLFVHhisto.SetLineWidth(2)
                                        ggLFVHhisto.Scale(10*scalefactors[0])
                                        vbfLFVHhisto.Scale(10*scalefactors[1])
                                        mymax = max(ggLFVHhisto.GetBinContent(ggLFVHhisto.GetMaximumBin()), vbfLFVHhisto.GetBinContent(vbfLFVHhisto.GetMaximumBin()))
                                        if mymax > LFVStack.GetMaximum() :  
                                                LFVStack.SetMaximum(1.2*mymax )
                                                print mymax, LFVStack.GetMaximum() 
                                        #LFVStack.Draw('hist')
                                        ggLFVHhisto.Draw('SAMEHIST')
                                        vbfLFVHhisto.Draw('SAMEHIST')
                                        legend.AddEntry(ggLFVHhisto, 'LFV GG Higgs x 10')
                                        legend.AddEntry(vbfLFVHhisto, 'LFV VBF Higgs x 10')                                
                                        legend.Draw() 
                                        canvas.SetLogy(0)
                                        if i.endswith('Pt') : canvas.SetLogy(1)
                                        canvas.Update()
                                        canvas.SaveAs(filepath+'/mc_'+i+'.png')
                                        outfile.cd()
                                        canvas.Write()
                                        outfile.Close()
                                        
                                        os.system("rm "+filepath+'/'+i+".root")
