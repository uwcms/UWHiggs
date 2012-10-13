'''

Compare W and QCD fake rates in dimuon events.

Author: Evan K. Friis, UW

'''

import os
import sys
import rootpy.io as io
import rootpy.plotting.views as views
import rootpy.plotting as plotting
import ROOT
from FinalStateAnalysis.PlotTools.RebinView import RebinView

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "Usage: WVsQCDComp.py file1 [file2] [file3]..."
    if not os.path.exists("fake_study_plots"):
        os.makedirs("fake_study_plots")

    data_files = [io.open(x) for x in sys.argv[1:]]

    data = views.SumView(*data_files)
    qcd_view = views.SubdirectoryView(data, "qcd")
    wjets_view = views.SubdirectoryView(data, "wjets")
    ROOT.gStyle.SetPaintTextFormat("0.1f")

    qcd_denom_2d = qcd_view.Get("pt10/muonJetVsLeptonPt")
    qcd_num_2d = qcd_view.Get("pt10/pfidiso02/muonJetVsLeptonPt")
    qcd_eff_2d = qcd_num_2d.Clone()
    qcd_eff_2d.Divide(qcd_denom_2d)
    qcd_eff_2d.Scale(100)

    triptych = ROOT.TCanvas("asdf", "asdf", 1200, 600)
    triptych.Divide(2)

    wjets_denom_2d = wjets_view.Get("pt10/muonJetVsLeptonPt")
    wjets_num_2d = wjets_view.Get("pt10/pfidiso02/muonJetVsLeptonPt")
    wjets_eff_2d = wjets_num_2d.Clone()
    wjets_eff_2d.Divide(wjets_denom_2d)
    wjets_eff_2d.Scale(100)

    triptych.cd(1)
    wjets_denom_2d.Draw("text")
    #import pdb; pdb.set_trace()
    triptych.cd(2)
    qcd_denom_2d.Draw("text")
    triptych.Draw()
    triptych.SaveAs("fake_study_plots/2d_denoms.png")

    triptych.cd(1)
    wjets_num_2d.Draw("text")
    triptych.cd(2)
    qcd_num_2d.Draw("text")
    triptych.SaveAs("fake_study_plots/2d_nums.png")

    triptych.cd(1)
    wjets_eff_2d.Draw("text")
    triptych.cd(2)
    qcd_eff_2d.Draw("text")
    triptych.SaveAs("fake_study_plots/2d_effs.png")

    canvas = plotting.Canvas(name='adsf', title='asdf')
    canvas.SetLogy(True)

    jet_pt_binning = [10,12,15,20,25,30,40,50,70,100,150,200]
    mu_pt_binning = [10, 12, 15, 18, 21,25,30,40,50]
    jet_pt_binning = 1
    mu_pt_binning = 1

    wjets_jet_pt_view = RebinView(wjets_view, jet_pt_binning)
    qcd_jet_pt_view = RebinView(qcd_view, jet_pt_binning)


    wjets_mu_pt_view = RebinView(wjets_view, mu_pt_binning)
    qcd_mu_pt_view = RebinView(qcd_view, mu_pt_binning)

    def round_to_ints(x):
        new = x.Clone()
        new.Reset()
        for bin in range(x.GetNbinsX()+1):
            nentries = ROOT.TMath.Nint(x.GetBinContent(bin))
            center = x.GetBinLowEdge(bin) + 0.5*x.GetBinWidth(bin)
            for _ in range(nentries):
                new.Fill(center)
        return new

    def get_efficiency(view, denompath, numpath, color):
        num = round_to_ints(view.Get(numpath))
        denom = round_to_ints(view.Get(denompath))
        graph = ROOT.TGraphAsymmErrors(num, denom)
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(color)
        return graph

    # Compare efficiency of jets and muons for W and QCD
    wjets_jet_eff = get_efficiency(wjets_jet_pt_view,
                                   "pt10/muonJetPt", "pt10/pfidiso02/muonJetPt",
                                   ROOT.EColor.kRed)
    qcd_jet_eff = get_efficiency(qcd_jet_pt_view,
                                 "pt10/muonJetPt", "pt10/pfidiso02/muonJetPt",
                                 ROOT.EColor.kBlue)

    wjets_jet_eff.Draw("ape")
    wjets_jet_eff.GetHistogram().SetMinimum(1e-3)
    wjets_jet_eff.GetHistogram().SetMaximum(1)
    qcd_jet_eff.Draw('pe')


    canvas.SetLogy(True)
    canvas.SaveAs("fake_study_plots/jet_eff_comp.png")
    canvas.SetLogy(False)
    canvas.SaveAs("fake_study_plots/jet_eff_comp_lin.png")

    wjets_mu_eff = get_efficiency(wjets_mu_pt_view,
                                  "pt10/muonPt", "pt10/pfidiso02/muonPt",
                                  ROOT.EColor.kRed)
    qcd_mu_eff = get_efficiency(qcd_mu_pt_view,
                                "pt10/muonPt", "pt10/pfidiso02/muonPt",
                                ROOT.EColor.kBlue)

    wjets_mu_eff.Fit("expo")
    wjets_mu_eff.Draw("ape")
    wjets_mu_eff.GetHistogram().SetMinimum(1e-3)
    wjets_mu_eff.GetHistogram().SetMaximum(1)
    qcd_mu_eff.Fit("expo")
    qcd_mu_eff.Draw('pe')

    canvas.SetLogy(True)
    canvas.SaveAs("fake_study_plots/mu_eff_comp.png")
    canvas.SetLogy(False)
    canvas.SaveAs("fake_study_plots/mu_eff_comp_lin.png")
    canvas.SetLogy(True)

    # Now check mu fake rate if we reweight by jet pt
    w2d = wjets_view.Get("pt10/muonJetVsLeptonPt")
    w2d_num = wjets_view.Get("pt10/pfidiso02/muonJetVsLeptonPt")
    q2d = qcd_view.Get("pt10/muonJetVsLeptonPt")
    q2d_num = qcd_view.Get("pt10/pfidiso02/muonJetVsLeptonPt")

    print "Initial raw fake rates"
    print "W: %0.2f" % (100.*w2d_num.Integral()/w2d.Integral())
    print "Q: %0.2f" % (100.*q2d_num.Integral()/q2d.Integral())

    w2d_jet_pt_proj = w2d.ProjectionX()
    q2d_jet_pt_proj = q2d.ProjectionX()
    w2d_jet_pt_proj.Scale(1.0/w2d_jet_pt_proj.Integral())
    q2d_jet_pt_proj.Scale(1.0/q2d_jet_pt_proj.Integral())

    q2d_reweight = q2d.Clone()
    q2d_num_reweight = q2d_num.Clone()
    for xbin in range(1, q2d_reweight.GetNbinsX()+1):
        q_val = q2d_jet_pt_proj.GetBinContent(xbin)
        w_val = w2d_jet_pt_proj.GetBinContent(xbin)
        xbin_weight = w_val/q_val if q_val else 0
        for ybin in range(1, q2d_reweight.GetNbinsY()+1):
            current_denom = q2d_reweight.GetBinContent(xbin, ybin)
            q2d_reweight.SetBinContent(xbin, ybin, current_denom*xbin_weight)
            current_num = q2d_num.GetBinContent(xbin, ybin)
            q2d_num_reweight.SetBinContent(xbin, ybin, current_num*xbin_weight)

    # These should be the same now
    q2d_reweight_jetpt_xcheck = q2d_reweight.ProjectionX("_px_reweight_check")
    q2d_reweight_jetpt_xcheck.Scale(1.0/q2d_reweight_jetpt_xcheck.Integral())
    w2d_jet_pt_proj.Draw()
    q2d_reweight_jetpt_xcheck.SetMarkerColor(ROOT.EColor.kRed)
    q2d_reweight_jetpt_xcheck.Draw('same, pe')
    canvas.SaveAs("fake_study_plots/xcheck.png")

    print "Reweighted Q fake rates"
    print "Q: %0.2f" % (100.*q2d_num_reweight.Integral()/q2d_reweight.Integral())

    def norm(x):
        x.Scale(1.0/x.Integral())
    # Now compare muonPt distribution
    w2d_mu_pt_proj = w2d.ProjectionY()
    norm(w2d_mu_pt_proj)
    w2d_mu_pt_proj.Draw()
    q2d_mu_pt_proj = q2d.ProjectionY()
    norm(q2d_mu_pt_proj)
    q2d_mu_pt_proj.SetMarkerColor(ROOT.EColor.kRed)
    q2d_mu_pt_proj.Draw("same")
    q2d_rew_mu_pt_proj = q2d_reweight.ProjectionY()
    norm(q2d_rew_mu_pt_proj)
    q2d_rew_mu_pt_proj.SetMarkerColor(ROOT.EColor.kBlue)
    q2d_rew_mu_pt_proj.Draw("same")

    canvas.SaveAs("fake_study_plots/mu_pt_comp.png")

    def eff(x, y, color, rebin=1):
        x = RebinView.rebin(x, mu_pt_binning)
        y = RebinView.rebin(y, mu_pt_binning)
        x = round_to_ints(x)
        y = round_to_ints(y)
        graph = ROOT.TGraphAsymmErrors(y, x)
        graph.SetMarkerColor(color)
        return graph

    fit_func = ROOT.TF1("expoconst", "[0]*TMath::Exp(-[1]*x)+[2]", 3)

    # Now compare mu pt effs
    print "wjets mu eff"
    w2d_mu_eff_denom = w2d.ProjectionY("w_d")
    w2d_mu_eff_num = w2d_num.ProjectionY("w_n")
    w2d_mu_eff_denom.Draw()
    w2d_mu_eff_num.SetMarkerColor(ROOT.EColor.kRed)
    w2d_mu_eff_num.Draw("same")
    canvas.SaveAs("fake_study_plots/debug.png")
    w2d_mu_eff = eff(w2d_mu_eff_denom, w2d_mu_eff_num, ROOT.EColor.kBlack, 2)
    w2d_mu_eff.Draw('ape')
    w2d_mu_eff.Fit('expo')
    w2d_mu_eff.GetHistogram().SetMaximum(1)
    w2d_mu_eff.GetHistogram().SetMinimum(1e-3)

    print "qcd mu eff"
    q2d_mu_eff_denom = q2d.ProjectionY()
    q2d_mu_eff_num = q2d_num.ProjectionY()
    q2d_mu_eff = eff(q2d_mu_eff_denom, q2d_mu_eff_num, ROOT.EColor.kRed, 2)
    q2d_mu_eff.Fit("expo")
    q2d_mu_eff.Draw('pe')

    print "qcd rew mu eff"
    q2d_rew_mu_eff_denom = q2d_reweight.ProjectionY()
    q2d_rew_mu_eff_num = q2d_num_reweight.ProjectionY()
    q2d_rew_mu_eff = eff(q2d_rew_mu_eff_denom, q2d_rew_mu_eff_num, ROOT.EColor.kBlue, 2)
    q2d_rew_mu_eff.Fit("expo")
    q2d_rew_mu_eff.Draw('pe')



    canvas.SaveAs("fake_study_plots/reweighted_mu_eff.png")


    # Compare mu denominators
    def norm(x):
        x.Scale(1.0/x.Integral())
        x = RebinView.rebin(x, mu_pt_binning)
        return x

    w2d_mu_eff_denom = norm(w2d_mu_eff_denom)
    q2d_mu_eff_denom = norm(q2d_mu_eff_denom)
    q2d_rew_mu_eff_denom = norm(q2d_rew_mu_eff_denom)

    def color(x, color):
        x.SetMarkerColor(getattr(ROOT.EColor, "k"+color))
    color(w2d_mu_eff_denom, "Black")
    color(q2d_mu_eff_denom, "Red")
    color(q2d_rew_mu_eff_denom, "Blue")

    w2d_mu_eff_denom.Draw()
    q2d_mu_eff_denom.Draw('same')
    q2d_rew_mu_eff_denom.Draw('same')

    canvas.SaveAs("fake_study_plots/reweighted_mu_denoms.png")
    canvas.SetLogy(False)
    canvas.SaveAs("fake_study_plots/reweighted_mu_denoms_lin.png")


    # Compare mu numerators

    w2d_mu_eff_num= norm(w2d_mu_eff_num)
    q2d_mu_eff_num= norm(q2d_mu_eff_num)
    q2d_rew_mu_eff_num=norm(q2d_rew_mu_eff_num)

    def color(x, color):
        x.SetMarkerColor(getattr(ROOT.EColor, "k"+color))
    color(w2d_mu_eff_num, "Black")
    color(q2d_mu_eff_num, "Red")
    color(q2d_rew_mu_eff_num, "Blue")

    w2d_mu_eff_num.Draw()
    w2d_mu_eff_num.SetMaximum(2*w2d_mu_eff_num.GetMaximum())
    q2d_mu_eff_num.Draw('same')
    q2d_rew_mu_eff_num.Draw('same')

    canvas.SetLogy(True)
    canvas.SaveAs("fake_study_plots/reweighted_mu_nums.png")
    canvas.SetLogy(False)
    canvas.SaveAs("fake_study_plots/reweighted_mu_nums_lin.png")

    # We can try reweighting the other direction - make the mu PT match, then
    # show jet pt fake rates

    w2d_mu_pt_proj = w2d.ProjectionY()
    q2d_mu_pt_proj = q2d.ProjectionY()
    w2d_mu_pt_proj.Scale(1.0/w2d_mu_pt_proj.Integral())
    q2d_mu_pt_proj.Scale(1.0/q2d_mu_pt_proj.Integral())

    q2d_mu_reweight = q2d.Clone()
    q2d_num_mu_reweight = q2d_num.Clone()
    for ybin in range(1, q2d_mu_reweight.GetNbinsY()+1):
        q_val = q2d_mu_pt_proj.GetBinContent(ybin)
        w_val = w2d_mu_pt_proj.GetBinContent(ybin)
        ybin_weight = w_val/q_val if q_val else 0
        for xbin in range(1, q2d_mu_reweight.GetNbinsX()+1):
            current_denom = q2d_mu_reweight.GetBinContent(xbin, ybin)
            q2d_mu_reweight.SetBinContent(xbin, ybin, current_denom*ybin_weight)
            current_num = q2d_num.GetBinContent(xbin, ybin)
            q2d_num_mu_reweight.SetBinContent(xbin, ybin, current_num*ybin_weight)

    # Now compare mu pt effs
    print "wjets jet eff"
    w2d_jet_eff_denom = w2d.ProjectionX("w_d_x")
    w2d_jet_eff_num = w2d_num.ProjectionX("w_n_x")
    w2d_jet_eff_denom.Draw()
    w2d_jet_eff_num.SetMarkerColor(ROOT.EColor.kRed)
    w2d_jet_eff_num.Draw("same")
    canvas.SaveAs("fake_study_plots/debug.png")
    w2d_jet_eff = eff(w2d_jet_eff_denom, w2d_jet_eff_num, ROOT.EColor.kBlack, 2)
    w2d_jet_eff.Draw('ape')
    w2d_jet_eff.GetHistogram().SetMaximum(1)
    w2d_jet_eff.GetHistogram().SetMinimum(1e-3)

    print "qcd jet eff"
    q2d_jet_eff_denom = q2d.ProjectionX()
    q2d_jet_eff_num = q2d_num.ProjectionX()
    q2d_jet_eff = eff(q2d_jet_eff_denom, q2d_jet_eff_num, ROOT.EColor.kRed, 2)
    q2d_jet_eff.Draw('pe')

    print "qcd rew jet eff"
    q2d_rew_jet_eff_denom = q2d_mu_reweight.ProjectionX()
    q2d_rew_jet_eff_num = q2d_num_mu_reweight.ProjectionX()
    q2d_rew_jet_eff = eff(q2d_rew_jet_eff_denom, q2d_rew_jet_eff_num, ROOT.EColor.kBlue, 2)
    q2d_rew_jet_eff.Draw('pe')

    canvas.SaveAs("fake_study_plots/reweighted_jet_eff.png")
