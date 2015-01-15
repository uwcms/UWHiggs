import ROOT 
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.SetBatch()

file_data = ROOT.TFile('results/newNtuple_3Dec/efakerate_fits/e_os_eLoose_eTigh_e3AbsEta.corrected_inputs.root')

data_uncorrected = file_data.Get('numerator_uncorr')
den_data_uncorrected = file_data.Get('denominator_uncorr')
data_corrected = file_data.Get('numerator')
den_data_corrected = file_data.Get('denominator')

data_uncorrected.Sumw2()
data_corrected.Sumw2()

data_uncorrected.Divide(den_data_uncorrected)
data_corrected.Divide(den_data_corrected)

c= ROOT.TCanvas("c","c", 800, 1000)
c.Draw()
c.SetGridx(1)
c.SetGridy(1)

pad1 = ROOT.TPad("pad1", "pad1", 0, 0.2, 1, 1)
pad1.Draw()
pad1.cd()
pad1.SetGridx(1)
pad1.SetGridy(1)
data_uncorrected.Draw('PE')
data_uncorrected.SetMarkerStyle(23)
data_uncorrected.GetXaxis().SetRangeUser(0, 2.2)
data_uncorrected.GetYaxis().SetRangeUser(0, 1.2)
data_uncorrected.SetMarkerColor(4)
data_uncorrected.SetLineColor(4)
data_uncorrected.GetYaxis().SetTitle("fakerate")
data_uncorrected.GetYaxis().SetTitleOffset(1.3)
data_uncorrected.GetXaxis().SetTitle("e |#eta|")
data_corrected.Draw('SAMEPE')
data_corrected.SetMarkerStyle(22)

file_MC  = ROOT.TFile('results/newNtuple_3Dec/efakerate_fits_MC/e_os_eLoose_eTigh_e3AbsEta_plot.root')

cMC = file_MC.Get('asdf')
dataMC = cMC.GetPrimitive('hxy_xy_data')

pad1.cd()
dataMC.Draw('SAMEP')
#dataMC.GetYaxis().SetRangeUser(0,1)


dataMC.SetMarkerStyle(20)
dataMC.SetMarkerColor(2)
dataMC.SetLineColor(2)

leg=ROOT.TLegend(0.3,0.85,0.7,0.7)
leg.SetFillColor(0)
leg.AddEntry(data_corrected, 'data, bkg corrected', 'lp')
leg.AddEntry(data_uncorrected, 'data', 'lp')
leg.AddEntry(dataMC, 'Z+jets MC', 'lp')
leg.Draw()

c.cd()

pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1,0.2)
pad2.Draw()
pad2.cd()
pad2.SetGridx(1)
pad2.SetGridy(1)

denom=ROOT.TH1F("denom", "denom", data_corrected.GetXaxis().GetNbins(), data_corrected.GetXaxis().GetXmin(), data_corrected.GetXaxis().GetXmax())
for bin in range(1, dataMC.GetN()+1):
    denom.Fill( dataMC.GetX()[bin-1],dataMC.GetY()[bin-1])
    x = dataMC.GetX()[bin-1]
    err = math.sqrt(math.pow(data_corrected.GetBinError(bin)/dataMC.GetY()[bin-1],2)+math.pow(dataMC.GetErrorY(bin-1)/dataMC.GetY()[bin-1],2)*math.pow(data_corrected.GetBinContent(bin)/dataMC.GetY()[bin-1],2))

    print x, data_corrected.GetBinContent(bin)/dataMC.GetY()[bin-1] , err, data_corrected.GetBinError(bin), dataMC.GetErrorY(bin-1) 
    denom.SetBinError(denom.FindBin(x), err)

ratio= data_corrected.Clone()
ratio.Divide(denom)
ratio.SetMarkerStyle(22)
ratio.Draw('PE1')
ratio.GetXaxis().SetRangeUser(0,2.2)
ratio.GetYaxis().SetTitle('data/mc')
ratio.GetYaxis().SetRangeUser(0,5)
ratio.GetYaxis().SetTitleSize(0.1)
ratio.GetYaxis().SetLabelSize(.1)
ratio.GetXaxis().SetTitleSize(.1)
ratio.GetXaxis().SetLabelSize(.1)
ratio.GetYaxis().SetNdivisions(504, False)
ratio.GetYaxis().SetTitleOffset(0.3)

myline = ROOT.TLine(0.05, 1, 2.29, 1)
myline.SetLineColor(2)
myline.SetLineWidth(2)
myline.SetLineStyle(2)
myline.Draw()


c.Update()
c.SaveAs('eFakerateComparison_Eta.pdf')


file_data.Close()
file_MC.Close()

file_data = ROOT.TFile('results/newNtuple_3Dec/efakerate_fits/e_os_eLoose_eTigh_e3Pt.corrected_inputs.root')

data_uncorrected = file_data.Get('numerator_uncorr')
den_data_uncorrected = file_data.Get('denominator_uncorr')
data_corrected = file_data.Get('numerator')
den_data_corrected = file_data.Get('denominator')

data_uncorrected.Sumw2()
data_corrected.Sumw2()

data_uncorrected.Divide(den_data_uncorrected)
data_corrected.Divide(den_data_corrected)

c= ROOT.TCanvas("c","c", 800, 1000)
c.Draw()
c.SetGridx(1)
c.SetGridy(1)
pad1 = ROOT.TPad("pad1", "pad1", 0, 0.2, 1, 1)
pad1.Draw()
pad1.cd()
pad1.SetGridx(1)
pad1.SetGridy(1)
pad1.cd()
data_uncorrected.Draw('SAMEPE')
data_uncorrected.GetXaxis().SetRangeUser(0, 2.3)
data_uncorrected.GetYaxis().SetRangeUser(0, 1.2)
data_uncorrected.GetYaxis().SetTitle("fakerate")
data_uncorrected.GetYaxis().SetTitleOffset(1.3)
data_uncorrected.GetXaxis().SetTitle("e p_{T} (GeV)")
data_uncorrected.SetMarkerStyle(23)
data_uncorrected.SetMarkerColor(4)
data_uncorrected.SetLineColor(4)
data_corrected.Draw('PE')
data_corrected.SetMarkerStyle(22)

file_MC  = ROOT.TFile('results/newNtuple_3Dec/efakerate_fits_MC/e_os_eLoose_eTigh_e3Pt_plot.root')

cMC = file_MC.Get('asdf')
dataMC = cMC.GetPrimitive('hxy_xy_data')

pad1.cd()
dataMC.Draw('SAMEP')
#dataMC.GetYaxis().SetRangeUser(0,1)


dataMC.SetMarkerStyle(20)
dataMC.SetMarkerColor(2)
dataMC.SetLineColor(2)

leg=ROOT.TLegend(0.3,0.85,0.7,0.7)
leg.SetFillColor(0)
leg.AddEntry(data_corrected, 'data, bkg corrected', 'lp')
leg.AddEntry(data_uncorrected, 'data', 'lp')
leg.AddEntry(dataMC, 'Z+jets MC', 'lp')
leg.Draw()
c.cd()

pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1,0.2)
pad2.Draw()
pad2.cd()
pad2.SetGridx(1)
pad2.SetGridy(1)

halfbinwidths=[5,5,5,5,5,10,10,15,24.7]

denom=ROOT.TGraphErrors(dataMC.GetN())

for ibin in range(1, dataMC.GetN()+1):
    
    denom.SetPoint(ibin-1, dataMC.GetX()[ibin-1],data_corrected.GetBinContent(ibin)/dataMC.GetY()[ibin-1])
    x = dataMC.GetX()[ibin-1]
    err = math.sqrt(math.pow(data_corrected.GetBinError(ibin)/dataMC.GetY()[ibin-1],2)+math.pow(dataMC.GetErrorY(ibin-1)/dataMC.GetY()[ibin-1],2)*math.pow(data_corrected.GetBinContent(ibin)/dataMC.GetY()[ibin-1],2))
    print x, data_corrected.GetBinContent(ibin)/dataMC.GetY()[ibin-1] , err, data_corrected.GetBinError(ibin), dataMC.GetErrorY(ibin-1)
    denom.SetPointError(ibin-1, halfbinwidths[ibin-1], err)



ratio=  denom
#ratio.Divide(denom)
ratio.SetMarkerStyle(22)
ratio.Draw('APE')
ratio.GetXaxis().SetRangeUser(30,200)
ratio.GetYaxis().SetTitle('data/mc')
ratio.GetYaxis().SetRangeUser(0,4)
ratio.GetYaxis().SetTitleSize(0.15)
ratio.GetYaxis().SetLabelSize(.15)
ratio.GetXaxis().SetTitleSize(.1)
ratio.GetXaxis().SetLabelSize(.1)
ratio.GetYaxis().SetNdivisions(504, False)

line = ROOT.TLine(35, 1, 190, 1)
line.SetLineColor(2)
line.SetLineWidth(2)
line.SetLineStyle(2)
line.Draw()

ratio.GetYaxis().SetTitleOffset(0.3)
pad2.Update()

c.Update()
c.SaveAs('eFakerateComparison_Pt.pdf')


file_data.Close()
file_MC.Close()

