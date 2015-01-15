import ROOT
import sys
import os
import math

try:
    massrange = sys.argv[1].split(',')
    print massrange
except:
        
    print 'please give me the mass range in the format 50,300 '

ROOT.gROOT.SetBatch()        # don't pop up canvases

ROOT.gROOT.SetStyle('Plain') # white background
ROOT.gStyle.SetOptStat(0)    
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)


   
samples = [
    "zjetsother",
    "diboson", 
    "SMVBF126", 
    "singlet",
    "SMGG126",
    "ttbar",
    "ztautau",
    "wplusjets",
    "LFVGG",
    "LFVVBF"
]

names = {
    "zjetsother" : "zjetsother",
    "diboson"    : "diboson"   ,
    "SMVBF126"   : "SMH"       ,
    "SMGG126"    : "SMH"       ,
    "singlet"    : "singlet"   ,
    "ttbar"      : "ttbar"     ,
    "ztautau"    : "ztautau"   ,
    "wplusjets"  : "wplusjets"     ,
    "LFVGG"      : "LFVH"      ,
    "LFVVBF"     : "LFVH"      
}


mymapper ={
    "wplusjets"  : ["W $+$ jets "],
    "ztautau"    : ["$ Z \\rightarrow \\tau\\tau$ "],
    "diboson"    : ["EWK Diboson "],
    "zjetsother" : ["$Z \\rightarrow ee, \\mu\\mu$"],
    "ttbar"      : ["$t\\bar{t}$ "],
    "singlet"    : ["$t$, $\\bar{t}$ "],
    "SMH"        : ["SM Higgs Background "],
    "tot"        : ["Sum of Backgrounds "],
    "LFVH"       : ["LFV Higgs Signal "],
}

dirname = ['gg0etau', 'boostetau', 'vbfetau']
     
tPt = [30, 35, 40, 45, 50, 60, 70, 80, 90]
ePt = [30, 35, 40, 45, 50, 60, 70, 80, 90]
tMtToPfMet = [5, 10, 20, 30, 35, 40, 50, 60, 70, 80, 90]

dphi = [0.50, 1.0, 1.50, 2.20, 2.40, 2.70, 3.00]

vbf_mass = [200, 300, 400, 450, 500, 550, 600, 700, 800, 900]
vbf_deta = [2.0, 2.5, 3.0, 3.5, 4.0]




for jet in range (0, 3) :
    file1 = ROOT.TFile.Open('plots/newNtuple_5Nov/lfvet/shapes.%s.root' %(str(jet)))
    cuts= [
        ['selected', 'tPt30', 'tPt35', 'tPt40', 'tPt45','tPt50','tPt60','tPt70','tPt80','tPt90',
         'ePt30', 'ePt35', 'ePt40', 'ePt45','ePt50','ePt60','ePt70','ePt80','ePt90',
         'dphi0.50','dphi1.00','dphi1.50','dphi2.20','dphi2.40','dphi2.70', 'dphi3.00',
         'tMtToPfMet5','tMtToPfMet10','tMtToPfMet20','tMtToPfMet30','tMtToPfMet35','tMtToPfMet40', 'tMtToPfMet50', 'tMtToPfMet60', 'tMtToPfMet70' , 'tMtToPfMet80', 'tMtToPfMet90' ],
        ['selected', 'tPt30', 'tPt35', 'tPt40', 'tPt45','tPt50','tPt60','tPt70','tPt80','tPt90',
         'ePt30', 'ePt35', 'ePt40', 'ePt45', 'ePt50','ePt60','ePt70','ePt80','ePt90',
         'tMtToPfMet5','tMtToPfMet10','tMtToPfMet20','tMtToPfMet30','tMtToPfMet35','tMtToPfMet40' , 'tMtToPfMet50', 'tMtToPfMet60', 'tMtToPfMet70' , 'tMtToPfMet80', 'tMtToPfMet90'],
        ['selected', 'tPt30', 'tPt35', 'tPt40', 'tPt45','tPt50','tPt60','tPt70','tPt80','tPt90',
         'ePt30', 'ePt35', 'ePt40', 'ePt45', 'ePt50','ePt60','ePt70','ePt80','ePt90',
         'tMtToPfMet5','tMtToPfMet10','tMtToPfMet20','tMtToPfMet30','tMtToPfMet35','tMtToPfMet40', 'tMtToPfMet50', 'tMtToPfMet60', 'tMtToPfMet70', 'tMtToPfMet80', 'tMtToPfMet90',
         'vbf_mass200', 'vbf_mass300', 'vbf_mass400', 'vbf_mass450', 'vbf_mass500', 'vbf_mass550',  'vbf_mass600',  'vbf_mass700',  'vbf_mass800',  'vbf_mass900', 
         'vbf_deta2.0','vbf_deta2.5','vbf_deta3.0','vbf_deta3.5', 'vbf_deta4.0' ]
    ]
    
    
    tPtcount =0
    ePtcount =0
    dphicount=0
    tMtcount =0 
    vbfmasscount =0
    vbfdetacount =0

    tPt_g  = ROOT.TGraphErrors(len(tPt))
    ePt_g  = ROOT.TGraphErrors(len(tPt))
    dphi_g = ROOT.TGraphErrors(len(dphi))
    tMt_g  = ROOT.TGraphErrors(len(tMtToPfMet))
    vbfMass_g = ROOT.TGraphErrors(len(vbf_mass))
    vbfdeta_g = ROOT.TGraphErrors(len(vbf_deta))

    tPt_ratio  = ROOT.TGraphErrors(len(tPt))
    ePt_ratio  = ROOT.TGraphErrors(len(tPt))
    dphi_ratio = ROOT.TGraphErrors(len(dphi))
    tMt_ratio  = ROOT.TGraphErrors(len(tMtToPfMet))
    vbfMass_ratio = ROOT.TGraphErrors(len(vbf_mass))
    vbfdeta_ratio = ROOT.TGraphErrors(len(vbf_deta))


    for cut in cuts[jet]:
        mydir = file1.Get("%s_%s" %(dirname[jet], cut))
        #mylist=mydir.GetListOfKeys()
        totbackground =0.
        totbackgrounderr2=0.
        sig_int =0.
        sig_err2 =0. 

        


        for sample in samples:
            print dirname[jet], cut, mydir.GetName(), sample
            histo = mydir.Get(sample)
            integral = 0 
            err2=0
        
            for i in range(histo.GetXaxis().FindBin(float(massrange[0])), histo.GetXaxis().FindBin(float(massrange[1]))+1):
                integral += histo.GetBinContent(i) 
                err2 += histo.GetBinError(i)*histo.GetBinError(i)
            
            if not 'LFV' in sample:
                totbackground += integral 
                totbackgrounderr2 += err2
        
            else:                 
                sig_int =integral
                sig_err2 =err2

                
        den = (sig_int+totbackground)
        #print "_".join([dirname[jet], cut]), sig_int/(sig_int+totbackground) , math.sqrt((sig_err2+totbackgrounderr2)* pow(sig_int,2) / pow(den,4))
        if den==0 : 
            if 'tPt' in cut :  tPtcount+=1
            if 'ePt' in cut : ePtcount+=1
            if 'dphi' in cut : dphicount+=1 
            if 'tMtToPfMet' in cut : tMtcount+=1
            if 'vbf_mass' in cut :vbfmasscount+=1
            if 'vbf_deta' in cut : vbfdetacount+=1       
            continue

        error = math.sqrt((totbackgrounderr2 * pow(sig_int,2) +sig_err2*(4*pow(den,2) - pow(sig_int,2) ))/4*pow(den,3))
        if 'tPt' in cut : 
            tPt_g.SetPoint(tPtcount, tPt[tPtcount], sig_int/math.sqrt(sig_int+totbackground))
            tPt_g.SetPointError(tPtcount, 0, math.sqrt((sig_err2+totbackgrounderr2)* pow(sig_int,2) / pow(den,4)))
            tPt_ratio.SetPoint(tPtcount, tPt[tPtcount], sig_int/error)
            tPt_ratio.SetPointError(tPtcount, 0,  math.sqrt(sig_err2)/error)
            tPtcount+=1
        if 'ePt' in cut : 
            print ePtcount, len(ePt), ePt[ePtcount]
            ePt_g.SetPoint(ePtcount, ePt[ePtcount], sig_int/math.sqrt(sig_int+totbackground))
            ePt_g.SetPointError(ePtcount, 0, math.sqrt((sig_err2+totbackgrounderr2)* pow(sig_int,2) / pow(den,4)))
            ePt_ratio.SetPoint(ePtcount, ePt[ePtcount], sig_int/error)
            ePt_ratio.SetPointError(ePtcount, 0,  math.sqrt(sig_err2)/error)
            ePtcount+=1
        if 'dphi' in cut : 
            dphi_g.SetPoint(dphicount, dphi[dphicount], sig_int/math.sqrt(sig_int+totbackground))
            dphi_g.SetPointError(dphicount, 0, math.sqrt((sig_err2+totbackgrounderr2)* pow(sig_int,2) / pow(den,4)))
            dphi_ratio.SetPoint(dphicount, dphi[dphicount], sig_int/error)
            dphi_ratio.SetPointError(dphicount, 0,  math.sqrt(sig_err2)/error)
            dphicount+=1 
        if 'tMtToPfMet' in cut : 
            tMt_g.SetPoint(tMtcount, tMtToPfMet[tMtcount], sig_int/math.sqrt(sig_int+totbackground))
            tMt_g.SetPointError(tMtcount, 0, math.sqrt((sig_err2+totbackgrounderr2)* pow(sig_int,2) / pow(den,4)))
            tMt_ratio.SetPoint(tMtcount, tMtToPfMet[tMtcount], sig_int/error)
            tMt_ratio.SetPointError(tMtcount, 0,  math.sqrt(sig_err2)/error)
            tMtcount+=1
        if 'vbf_mass' in cut : 
            vbfMass_g.SetPoint(vbfmasscount, vbf_mass[vbfmasscount], sig_int/math.sqrt(sig_int+totbackground))
            vbfMass_g.SetPointError(vbfmasscount, 0, math.sqrt((sig_err2+totbackgrounderr2)* pow(sig_int,2) / pow(den,4)))
            vbfMass_ratio.SetPoint(vbfmasscount, vbf_mass[vbfmasscount], sig_int/error)
            vbfMass_ratio.SetPointError(vbfmasscount, 0, math.sqrt(sig_err2)/error)
            vbfmasscount+=1
        if 'vbf_deta' in cut : 
            vbfdeta_g.SetPoint(vbfdetacount, vbf_deta[vbfdetacount], sig_int/math.sqrt(sig_int+totbackground))
            vbfdeta_g.SetPointError(vbfdetacount, 0, math.sqrt((sig_err2+totbackgrounderr2)* pow(sig_int,2) / pow(den,4)))
            vbfdeta_ratio.SetPoint(vbfdetacount, vbf_deta[vbfdetacount], sig_int/error)
            vbfdeta_ratio.SetPointError(vbfdetacount, 0, math.sqrt(sig_err2)/error)
            vbfdetacount+=1             
    
    c=ROOT.TCanvas("c", "c", 800, 800) 
    massInterval = massrange[0]+'_'+massrange[1]
    c.Draw()
    c.Divide(1,2)
    
    c.SetGridx(1)
    c.SetGridy(1)
    if jet==0:
        c.cd(1)
        tPt_g.SetTitle("0 jet, #tau p_{T} threshold")
        tPt_g.Draw("AP")
        tPt_g.GetXaxis().SetTitle("#tau p_{T} (GeV)") 
        tPt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        tPt_g.SetMarkerStyle(20)
        ##tPt_g.SetMarkerSize(0.3)
        c.cd(2)
        tPt_ratio.Draw("AP")
        tPt_ratio.GetXaxis().SetTitle("#tau p_{T} (GeV)") 
        tPt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        tPt_ratio.SetMarkerStyle(20)
        
        c.Update()
        c.SaveAs("0jet_tPt_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        ePt_g.SetTitle("0 jet, e p_{T} threshold")
        ePt_g.Draw("AP")
        ePt_g.GetXaxis().SetTitle("e p_{T} (GeV)") 
        ePt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        ePt_g.SetMarkerStyle(20)
        c.cd(2)
        ePt_ratio.Draw("AP")
        ePt_ratio.GetXaxis().SetTitle("e p_{T} (GeV)") 
        ePt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        ePt_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("0jet_ePt_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        dphi_g.SetTitle("0 jet, #Delta#Phi")
        dphi_g.Draw("AP")
        dphi_g.GetXaxis().SetTitle("e - #tau #Delta#Phi")
        dphi_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        dphi_g.SetMarkerStyle(20)
        c.cd(2)
        dphi_ratio.Draw("AP")
        dphi_ratio.GetXaxis().SetTitle("e - #tau #Delta#Phi") 
        dphi_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        dphi_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("0jet_dphi_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        tMt_g.SetTitle("0 jet, #tau Met M_{T}")
        tMt_g.Draw("AP")
        tMt_g.GetXaxis().SetTitle("#tau Met M_{T}")
        tMt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        tMt_g.SetMarkerStyle(20)
        c.cd(2)
        tMt_ratio.Draw("AP")
        tMt_ratio.GetXaxis().SetTitle("#tau Met M_{T}") 
        tMt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        tMt_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("0jet_tMt_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
    if jet==1:
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        tPt_g.SetTitle("1 jet, #tau p_{T} threshold")
        tPt_g.Draw("AP")
        tPt_g.GetXaxis().SetTitle("#tau p_{T} (GeV)") 
        tPt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        tPt_g.SetMarkerStyle(20)
        c.cd(2)
        tPt_ratio.Draw("AP")
        tPt_ratio.GetXaxis().SetTitle("#tau p_{T} (GeV)") 
        tPt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        tPt_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("1jet_tPt_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        ePt_g.SetTitle("1 jet, e p_{T} threshold")
        ePt_g.Draw("AP")
        ePt_g.GetXaxis().SetTitle("e p_{T} (GeV)") 
        ePt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        ePt_g.SetMarkerStyle(20)
        c.cd(2)
        ePt_ratio.Draw("AP")
        ePt_ratio.GetXaxis().SetTitle("e p_{T} (GeV)") 
        ePt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        ePt_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("1jet_ePt_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        tMt_g.SetTitle("1 jet, #tau Met M_{T}")
        tMt_g.Draw("AP")
        tMt_g.GetXaxis().SetTitle("#tau Met M_{T}")
        tMt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        tMt_g.SetMarkerStyle(20)
        c.cd(2)
        tMt_ratio.Draw("AP")
        tMt_ratio.GetXaxis().SetTitle("#tau Met M_{T}") 
        tMt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        tMt_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("1jet_tMt_%s.pdf" %massInterval)
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
    if jet==2:
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        tPt_g.SetTitle("2 jet, #tau p_{T} threshold")
        tPt_g.Draw("AP")
        tPt_g.GetXaxis().SetTitle("#tau p_{T} (GeV)") 
        tPt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        tPt_g.SetMarkerStyle(20)
        c.cd(2)
        tPt_ratio.Draw("AP")
        tPt_ratio.GetXaxis().SetTitle("#tau p_{T} (GeV)") 
        tPt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        tPt_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("2jet_tPt_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        ePt_g.SetTitle("2 jet, e p_{T} threshold")
        ePt_g.Draw("AP")
        ePt_g.GetXaxis().SetTitle("e p_{T} (GeV)") 
        ePt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        ePt_g.SetMarkerStyle(20)
        c.cd(2)
        ePt_ratio.Draw("AP")
        ePt_ratio.GetXaxis().SetTitle("e p_{T} (GeV)") 
        ePt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        ePt_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("2jet_ePt_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        tMt_g.SetTitle("2 jet, #tau Met M_{T}")
        tMt_g.Draw("AP")
        tMt_g.GetXaxis().SetTitle("#tau Met M_{T}")
        tMt_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        tMt_g.SetMarkerStyle(20)
        c.cd(2)
        tMt_ratio.Draw("AP")
        tMt_ratio.GetXaxis().SetTitle("#tau Met M_{T}") 
        tMt_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        tMt_ratio.SetMarkerStyle(20)
        c.Update()
        c.SaveAs("2jet_tMt_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        vbfMass_g.SetTitle("2 jet, vbf Mass")
        vbfMass_g.Draw("AP")
        vbfMass_g.GetXaxis().SetTitle("M_{jj}")
        vbfMass_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        vbfMass_g.SetMarkerStyle(20)
        c.cd(2)
        vbfMass_ratio.Draw("AP")
        vbfMass_ratio.GetXaxis().SetTitle("M_{jj}")
        vbfMass_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        vbfMass_ratio.SetMarkerStyle(20)

        c.Update()
        c.SaveAs("2jet_vbfMass_%s.pdf" %massInterval) 
        c.Clear()
        c.Divide(1,2)
        c.cd(1)
        vbfdeta_g.SetTitle("2 jet, #Delta#eta_{jj}")
        vbfdeta_g.Draw("AP")
        vbfdeta_g.GetXaxis().SetTitle("#Delta#eta_{jj}")
        vbfdeta_g.GetYaxis().SetTitle("S/#sqrt(S+B)") 
        vbfdeta_g.SetMarkerStyle(20)
        c.cd(2)
        vbfdeta_ratio.Draw("AP")
        vbfdeta_ratio.GetXaxis().SetTitle("2 jet, #Delta#eta_{jj}")
        vbfdeta_ratio.GetYaxis().SetTitle("S/#sigma(S/#sqrt(S+B))") 
        vbfdeta_ratio.SetMarkerStyle(20)

        c.Update()
        c.SaveAs("2jet_vbfDeta_%s.pdf" %massInterval) 
        c.Clear()
        

        file1.Close()
        
