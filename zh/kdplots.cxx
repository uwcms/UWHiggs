
/* AUTHOR: Stephane 
 Usage: kdplots(channel), where channel is the desired zh channel IN ALL CAPS. For example...
 kdplots("MMMT") will return normalized kinematic discriminator comparison plots for the MMMT channel.

 generates comparison plots (signal vs. WZ vs. Z+jets) for two kinematic discriminators, with each shape normalized. */

void kdplots(char* channel) {
  TFile *signal = TFile::Open(TString::Format("results/2012-10-24-8TeV-v1-Higgs/ZHAnalyze%s/VH_H2Tau_M-120.root",channel));
  TFile *wz = TFile::Open(TString::Format("results/2012-10-24-8TeV-v1-Higgs/ZHAnalyze%s/WZJetsTo3LNu_pythia.root",channel));
  TFile *zjets = TFile::Open(TString::Format("results/2012-10-24-8TeV-v1-Higgs/ZHAnalyze%s/Zjets_M50.root",channel));

  TH1F *kd1_signal = signal->Get("os/All_Passed/kinematicDiscriminant1");
  TH1F *kd2_signal = signal->Get("os/All_Passed/kinematicDiscriminant2");
  TH1F *kd1_zjets = zjets->Get("os/All_Passed/kinematicDiscriminant1");
  TH1F *kd2_zjets = zjets->Get("os/All_Passed/kinematicDiscriminant2");
  TH1F *kd1_wz = wz->Get("os/All_Passed/kinematicDiscriminant1");
  TH1F *kd2_wz = wz->Get("os/All_Passed/kinematicDiscriminant2");
  
  // Normalize histograms
  kd1_signal = kd1_signal->Scale(1/kd1_signal->Integral());
  kd2_signal = kd2_signal->Scale(1/kd2_signal->Integral());
  kd1_zjets = kd1_zjets->Scale(1/kd1_zjets->Integral());
  kd2_zjets = kd2_zjets->Scale(1/kd2_zjets->Integral());
  kd1_wz = kd1_wz->Scale(1/kd1_wz->Integral());
  kd2_wz = kd2_wz->Scale(1/kd2_wz->Integral());
 
  // Plot normalized histograms
  TCanvas *c1 = new TCanvas("c1","kd1");
  kd1_zjets->SetMarkerColor(2);
  kd1_wz->SetMarkerColor(4);
  kd1_signal->Draw();
  kd1_zjets->Draw("same");
  kd1_wz->Draw("same");
  TLegend *l1 = new TLegend(0.5,0.95,0.65,0.8,"");
  l1->AddEntry(kd1_zjets,"zjets","p");
  l1->AddEntry(kd1_wz,"wz","p");
  l1->AddEntry(kd1_signal, "VH_H2Tau(MH=120)", "p");
  l1->SetBorderSize(0);
  l1->SetFillColor(0);
  l1->Draw();
  c1->SaveAs(TString::Format("results/2012-10-24-8TeV-v1-Higgs/kd_plots/kd1_%s.png",channel));

  TCanvas *c2 = new TCanvas("c2","kd2");
  kd2_zjets->SetMarkerColor(2);
  kd2_wz->SetMarkerColor(4);
  kd2_signal->Draw();
  kd2_zjets->Draw("same");
  kd2_wz->Draw("same");
  TLegend *l2 = new TLegend(0.5,0.95,0.65,0.8,"");
  l2->AddEntry(kd2_zjets,"zjets","p");
  l2->AddEntry(kd2_wz,"wz","p");
  l2->AddEntry(kd2_signal, "VH_H2Tau(MH=120)", "p");
  l2->SetBorderSize(0);
  l2->SetFillColor(0);
  l2->Draw();
  c2->SaveAs(TString::Format("results/2012-10-24-8TeV-v1-Higgs/kd_plots/kd2_%s.png",channel));

}
