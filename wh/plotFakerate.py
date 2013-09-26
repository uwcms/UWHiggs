#!/usr/bin/env python
"""
Use: ./plotFakeRate.py myFakeRate bool where myFakeRate MM, EM, EE, MMT, MMM etc.
bool = True if you want to save the images of  plots
"""
print __doc__
import os
import sys
import shutil
import rootpy
import subprocess 

#rootpy.log.basic_config_colorized()
from rootpy.io import open as ropen, DoesNotExist
from rootpy.io import File as File
from rootpy.plotting import Hist, Hist2D
from glob import glob
from rootpy.plotting import Hist, Hist2D, Hist3D, HistStack, Legend, Canvas
from rootpy.plotting.views import NormalizeView, SumView, SubdirectoryView
from rootpy.interactive import wait
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptTitle(0)
from collections import deque

def GetContent(dir):
    tempList = dir.GetListOfKeys()
    retList = []
    for it in range(0,tempList.GetSize()):
       retList.append(tempList.At(it).ReadObj())
    return retList

def FindTrees( directory, dirName, objectList ):
    dirContent = GetContent(directory)
    for entry in dirContent:
        if type(entry) is ROOT.TDirectory or type(entry) is ROOT.TDirectoryFile:
            subdirName = os.path.join(dirName,entry.GetName())
            FindTrees(entry, subdirName,objectList)
        elif entry.InheritsFrom('TTree'):
            pathname = os.path.join(dirName,entry.GetName())
            objectList.append(pathname)
def FindTH1(directory, dirName, objectList):
    dirContent = GetContent(directory)
    for entry in dirContent:
        if type(entry) is ROOT.TDirectory or type(entry) is ROOT.TDirectoryFile:
            subdirName = os.path.join(dirName,entry.GetName())
 ##           FindDir(entry, subdirName,objectList)
 ##           print subdirName
        elif entry.InheritsFrom('TH1'):
            pathname = os.path.join(dirName,entry.GetName())
            objectList.append(pathname)

def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1", "True")

if len(sys.argv) < 2 or '-h' in sys.argv or '--help' in sys.argv:
    print __doc__
    sys.exit(1)

channel = sys.argv[-2]
saveplot = str2bool(sys.argv[-1])
print sys.argv
mydirpath = 'results/2013-May-14-8TeV/FakeRates'
mydir = ''.join([mydirpath,channel,'/'])
print 'Analysing Fakerates in dir ', mydir

files =  glob(mydir+'alldata.root')
##outputfiles = glob(mydir+'output.root')
##print outputfiles
print files
nfile =0

tfile=ropen(files[0])
dir=tfile.Get('qcd')
subdirlistQCD = GetContent(dir) #[ i for i in GetContent(dir) if i.GetName().startswith('pt10_') or i.GetName().startswith('pt20_')]  

total_view = SumView( *[ropen(i) for i in files] )
qcd_view   = SubdirectoryView(total_view, 'qcd')
if channel=='EE':
    wjet_view = SubdirectoryView(total_view, 'wjetsNoZmass')
else:
    wjet_view  = SubdirectoryView(total_view, 'wjets')

for m,n in enumerate (subdirlistQCD):
    
    thlistQCD = []

    FindTH1( subdirlistQCD[m], '', thlistQCD)

    print thlistQCD
    for i in thlistQCD:
        if not (i.find('pt10_') or i.find('pt20_')):
            continue
        histoQCD = qcd_view.Get('/'.join([subdirlistQCD[m].GetName(),i]))
        histoWJets = wjet_view.Get('/'.join([subdirlistQCD[m].GetName(),i]))
        #print '/'.join([subdirlistQCD[m].GetName(),i])
        if (histoQCD.Integral() ==0 or histoWJets.Integral()==0): continue
        if (i.find('JetQGLikelihoodID') != -1 or i.find('JetQGMVAID')!=-1):
            print "QCD %s, number of events %.0f, underflow %.1f, overflow %.1f" % \
                (i, histoQCD.Integral(), histoQCD.GetBinContent(0), histoQCD.GetBinContent(histoQCD.GetNbinsX()+1))
            print "WJets %s, number of events %.0f, underflow %.1f, overflow %.1f" % \
                (i, histoWJets.Integral(), histoWJets.GetBinContent(0), histoWJets.GetBinContent(histoWJets.GetNbinsX()+1))

        histoQCD.Scale(1./histoQCD.Integral())
        histoWJets.Scale(1./histoWJets.Integral())
        canvas = Canvas ()
        canvas.Draw()
        histoQCD.Draw()
        histoQCD.SetMarkerStyle(20)
        histoWJets.SetLineColor(2)
        histoWJets.SetMarkerColor(2)
        histoWJets.SetMarkerStyle(20)
        histoWJets.Draw('SAME')
        if (histoWJets.GetBinContent(histoWJets.GetMaximumBin())> histoQCD.GetBinContent(histoQCD.GetMaximumBin())):
            histoQCD.GetYaxis().SetRangeUser(0,histoWJets.GetBinContent(histoWJets.GetMaximumBin())*1.1)
        leg = Legend(2)
        leg.AddEntry(histoWJets)#, label='wjets')
        leg.AddEntry(histoQCD)#, label='QCD')
        leg.Draw()
        
        if  (i.find('Jetmult') != -1):
            ibin = 5
            bestbin =0
            bestsnratio=0.
            while ibin <= histoQCD.GetNbinsX():
                snratio=histoWJets.Integral(1,ibin)/ \
                         (histoQCD.Integral(1, ibin)+histoWJets.Integral(1,ibin)) 
                ibin=ibin+1
                if snratio > bestsnratio:
                    bestsnratio=snratio
                    bestbin=ibin-1
            print "Best sn ratio of %s %.2f at %f " % (i, bestsnratio, histoWJets.GetXaxis().GetBinCenter(bestbin))
                
        if  (i.find('JetptD') != -1):
            ibin = histoQCD.GetNbinsX()-10
            bestsnratio=0.
            bestbin=0
            while ibin > 0:
                snratio=histoWJets.Integral(ibin, histoQCD.GetNbinsX())/ \
                         (histoQCD.Integral(ibin,histoQCD.GetNbinsX())+histoWJets.Integral(ibin, histoQCD.GetNbinsX()))
                #print snratio
                ibin=ibin-1
                if snratio > bestsnratio:
                    bestsnratio=snratio
                    bestbin=ibin+1
            print "Best sn ratio of %s %.2f %.2f" % (i, bestsnratio, histoWJets.GetXaxis().GetBinCenter(bestbin))
                


        if not saveplot: 
            continue
        filelist = ['JetArea', 'EtaEtaMoment', 'EtaPhiMoment', 'EtaPhiSpread', 'PhiPhiMoment','JetWidth', 'JetWidthCorr', 'JetEtaPhiMomentOverSpread','JetDet', 'JetWidthOverSpread','JetEtaEtaPhiPhiDiff','JetWidthCorrOverSpread','JetPhiPhiMomentOverSpread','JetEtaEtaMomentOverSpread', 'MET', 'MET_Z', 'MET_NoZ', 'JetptD','Jetaxis', 'Jetmult', 'JetQG']
        
        ##filelist = ['JetptD','Jetaxis', 'Jetmult', 'JetQG']
        ofile=ropen(mydir+'outputfile.root', "RECREATE")
        ofile.cd()
        for j in filelist:
            #print i.find(j)
            if i.find(j)!=-1:
                canvas.SaveAs(''.join([mydir,subdirlistQCD[m].GetName(),'_',i,'.png']))
                canvas.Write()

        ofile.close()
##                mystr = ''.join(['cp ', mydir,subdirlistQCD[m].GetName(),'_',i,'.png', '  /afs/hep.wisc.edu/home/taroni/public_html/fakerateStudies', channel,'/'])
##                print mystr
##                proc = subprocess.Popen(mystr)
##                stdout, stderr = proc.communicate()
##                exitc  = proc.wait()

