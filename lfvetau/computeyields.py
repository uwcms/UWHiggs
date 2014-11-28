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

samples = [
    "zjetsother",
    "diboson", 
    "SMVBF126", 
    "singlet",
    "SMGG126",
    "ttbar",
    "ztautau",
    "fakes",
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
    "fakes"      : "fakes"     ,
    "LFVGG"      : "LFVH"      ,
    "LFVVBF"     : "LFVH"      
}


mymapper ={
    "fakes"      : ["Fakes "],
    "ztautau"    : ["$ Z \\rightarrow \\tau\\tau$ "],
    "diboson"    : ["EWK Diboson "],
    "zjetsother" : ["$Z \\rightarrow ee, \\mu\\mu$"],
    "ttbar"      : ["$t\\bar{t}$ "],
    "singlet"    : ["$t$, $\\bar{t}$ "],
    "SMH"        : ["SM Higgs Background "],
    "tot"        : ["Sum of Backgrounds "],
    "LFVH"       : ["LFV Higgs Signal "],
}

     
     
    
for jet in range (0, 3) :

    #norm_file = 'plots/newNtuple_5Nov/lfvet/unc.%s.vals' %(str(jet))
    #f = open(norm_file, 'r')
    #lines= f.readlines()
    #for line in lines :
    #    norm= line.split()
    #    print jet, norm, norm[1].split(',')


    #print jet
    file1 = ROOT.TFile.Open('plots/newNtuple_5Nov/lfvet/shapes.%s.root' %(str(jet)))
    #print 'plots/newNtuple_5Nov/lfvet/shapes.%s.root' %(str(jet))
    #print file1.GetListOfKeys()[0].GetName()
    mydir = file1.Get(file1.GetListOfKeys()[0].GetName())
    mylist=mydir.GetListOfKeys()
    totbackground =0.
    totbackgrounderr2=0.
    for sample in samples:
        #print sample
        histo = mydir.Get(sample)
        sublist = filter(lambda x :  sample in x.GetName()  and sample != x.GetName(), mylist) 
       
        integral = 0 
        err2=0

        for i in range(histo.GetXaxis().FindBin(float(massrange[0])), histo.GetXaxis().FindBin(float(massrange[1]))+1):
            integral += histo.GetBinContent(i) 
            err2 += histo.GetBinError(i)*histo.GetBinError(i)
        if integral != 0:
            print sample, integral, math.sqrt(err2), math.sqrt(err2)/integral
        else:
            print sample, integral, math.sqrt(err2)

        for histname in sublist:
            #print 'histo name', histname.GetName()
            systHistInt=0.
            h = mydir.Get(histname.GetName()).Clone()
            #print histo.Integral(), h.Integral()
            for i in range(histo.GetXaxis().FindBin(float(massrange[0])), histo.GetXaxis().FindBin(float(massrange[1]))+1):
                systHistInt+=pow(h.GetBinContent(i) - histo.GetBinContent(i), 2)
                #print sample, h.GetBinContent(i) ,  histo.GetBinContent(i), h.GetBinContent(i) - histo.GetBinContent(i)
                #print histname.GetName(), systHistInt
            #print systHistInt, integral
            err2+= systHistInt
            

                
        #print histo.FindBin(float(massrange[0])), histo.FindBin(float(massrange[1])), histo.GetXaxis().GetNbins()
        
        if integral != 0:
            print sample, integral, math.sqrt(err2), math.sqrt(err2)/integral
        else:
            print sample, integral, math.sqrt(err2),

        if not 'LFV' in sample:
            totbackground += integral 
            totbackgrounderr2 += err2
            #print totbackground, totbackgrounderr2
            
        try:
            print mymapper[names[sample]][1+3*jet],mymapper[names[sample]][3+3*jet]
            mymapper[names[sample]][1+3*jet]=float(mymapper[names[sample]][1+3*jet]) + integral
            mymapper[names[sample]][3+3*jet]=math.sqrt(pow(float(mymapper[names[sample]][3+3*jet]),2) + err2)

        except:
            mymapper[names[sample]].extend([ integral, "$\pm$" ,  math.sqrt(err2)])

    print 'tot', totbackground, math.sqrt(totbackgrounderr2)
    mylist = ["%.2f" % totbackground, "$\pm$" , "%.2f" % math.sqrt(totbackgrounderr2)]
    mymapper["tot"].extend(mylist)

        #error=ROOT.Double()
        #print histo.IntegralAndError(histo.FindBin(float(massrange[0])), histo.FindBin(float(massrange[1])), error), error

for k,v in mymapper.iteritems():
    for n,obj in enumerate(v) :
        if isinstance(obj, float): mymapper[k][n]=str('%.2f' %obj)

outfile = open('yields_MassRange_%s_%s.tex' %(str(massrange[0]),str(massrange[1])) ,'w') 
outfile.write('\\begin{table}[h] \n')
outfile.write('\\begin{tabular}{|c|rcl|rcl|rcl|}\n')

for k,v in mymapper.iteritems():
    print ' & '.join(v)+'\\\\\\hline\n'
    outfile.write(' & '.join(v)+'\\\\\n')

outfile.write('\\end{tabular}\n')
outfile.write('\\end{table}\n')

outfile.close()
