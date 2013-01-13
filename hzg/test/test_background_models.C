void test_background_models() {

  gROOT->ProcessLine(".X CMSStyle.C");
  gSystem->Load("libRooFitCore.so");

  RooWorkspace test("test");

  char fName[] = "/home/lgray/HZG_analysis/10JAN2013_vanilla/muon/ZGToLLG-8TeV/ZGToLLG_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM.root";
  char fName2[] = "/home/lgray/HZG_analysis/10JAN2013_vanilla/muon/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-8TeV/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM.root";

  TFile* f = TFile::Open(fName,"READ");  
  TFile* f2 = TFile::Open(fName2,"READ");
  TTree* tree = (TTree*)f->Get("selected_zg");
  TTree* tree2 = (TTree*)f2->Get("selected_zg");   

  test.factory("procWeight[0]");
  test.factory("puWeight[0]");
  test.factory("weight[0]");
  test.factory("Mzg[90,200]");
  test.factory("Mz[60,120]");
  test.factory("dMzg[0,25]");
  test.factory("dMz[0,25]");
  test.factory("r94cat[cat1=1,cat2=2,cat3=3,cat4=4]");
  test.defineSet("observables","Mzg,Mz,dMzg,dMz,r94cat,procWeight,puWeight");
  test.defineSet("observables_weight","Mzg,Mz,dMzg,dMz,r94cat,weight");

  RooArgSet* obs = test.set("observables");
  RooArgSet* wobs = test.set("observables_weight");

  char channel[] = "muon_cat4";
  char subscriptzg[] = "#mu#mu#gamma";
  char subscriptz[] = "#mu#mu";
  

  RooDataSet dllg( "zgdata","M_{Z#gamma} with Errors",
		   tree, *obs);
  RooConstVar llgtot("llgtot","llgtot",6.58816100000000000e+06);
  RooFormulaVar weightllg("weight","weight","@0*@1/@2",RooArgList(*test.var("puWeight"),
								  *test.var("procWeight"),
								  llgtot));
  dllg.addColumn(weightllg);    
  
  RooDataSet dzjets( "zjdata","M_{Z Jets} with Errors",
		     tree2, *obs);
  RooConstVar zjetstot("zjetstot","zjetstot",3.04255770000000000e+07);
  RooFormulaVar weightzjets("weight","weight","@0*@1/@2",RooArgList(*test.var("puWeight"),
								    *test.var("procWeight"),
								    zjetstot));
  dzjets.addColumn(weightzjets);

  dllg.append(dzjets);

  RooDataSet d( "data","M_{Z Jets} with Errors", &dllg, *wobs," Mz + Mzg > 185 && r94cat == 4 ","weight" );

  test.import(d);
  
  test.factory("RooGaussModel::MzgResoShape(Mzg,bias[120,90,150],sigma[1,0.01,10])");
  test.factory("RooDecay::MzgBkgShape(Mzg,tau[5,0,50],MzgResoShape,RooDecay::SingleSided)");

  /*
  test.factory("RooLandau::MzgTruthShape(Mzg,lanMPV[130,100,180],lanWidth[10,0.01,50])");
  test.factory("RooGaussian::MzgResoShape(Mzg,0,sigma[0.01,10])");
  test.factory("FCONV::MzgBkgShape(Mzg,MzgTruthShape,MzgResoShape)");
  */

  test.pdf("MzgBkgShape")->fitTo(*(test.data("data")),
				 RooFit::ConditionalObservables(*(test.var("dMzg"))),
				 RooFit::SumW2Error(kTRUE));
  
  TCanvas canv("test","test",600,600);

  tlx = TLatex();
  tlx.SetNDC();

  RooPlot* frame = test.var("Mzg")->frame(90,200,150);
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle(Form("M_{%s} (GeV)",subscriptzg));
  frame->GetYaxis()->SetTitle("Entries");
  test.data("data")->plotOn(frame);
  test.pdf("MzgBkgShape")->plotOn(frame,
				  RooFit::ProjWData(*(test.var("dMzg")),
						    *(test.data("data"))));
  
  frame->GetYaxis()->SetRangeUser(0,70);
  frame->Draw();
  tlx.DrawLatex(0.63,0.88,Form("#chi^{2} = %.3f",frame->chiSquare()));
  canv.Print(Form("bkgshapetest_test_zg_%s.pdf",channel));  
  canv.Clear();
  canv.SetLogy();
  frame->GetYaxis()->SetRangeUser(1e-3,70);
  frame->Draw();
  tlx.DrawLatex(0.63,0.88,Form("#chi^{2} = %.3f",frame->chiSquare()));
  canv.Print(Form("bkgshapetest_test_zg_log_%s.pdf",channel));
  canv.Clear();

  delete frame;  

  f->Close();
  f2->Close();
}

