void test_background_models() {

  gSystem->Load("libRooFitCore.so");

  RooWorkspace test("test");

  test.factory("Mzg[90,200]");
  test.factory("Mz[60,120]");
  test.factory("dMzg[0,25]");
  test.factory("dMz[0,25]");
  test.factory("r94cat[cat1=1,cat2=2,cat3=3,cat4=4]");
  test.defineSet("observables","Mzg,Mz,dMzg,dMz,r94cat");

  char channel[] = "electron_cat1";
  char subscriptzg[] = "ee#gamma";
  char subscriptz[] = "ee";
  char fName[] = "/home/lgray/HZG_analysis/06JAN2013_vanilla/electron/ZGToLLG-8TeV/ZGToLLG_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM.root";

  TFile* f = TFile::Open(fName,"READ");
  TTree* tree = (TTree*)f->Get("selected_zg");
  RooArgSet* obs = test.set("observables");

  std::cout << f << ' ' << tree << std::endl;

  RooDataSet d( "zgdata","M_{Z#gamma} with Errors",
		tree, *obs, "dMz/Mz < 1.0 && r94cat == 1 && Mzg + Mz > 185" );

  test.import(d);
  
  test.factory("RooGaussModel::MzgResoShape(Mzg,bias[120,90,150],sigma[1,0.01,10])");
  test.factory("RooDecay::MzgBkgShape(Mzg,tau[5,0,50],MzgResoShape,RooDecay::SingleSided)");

  /*
  test.factory("RooLandau::MzgTruthShape(Mzg,lanMPV[130,100,180],lanWidth[10,0.01,50])");
  test.factory("RooGaussian::MzgResoShape(Mzg,0,sigma[0.01,10])");
  test.factory("FCONV::MzgBkgShape(Mzg,MzgTruthShape,MzgResoShape)");
  */

  test.pdf("MzgBkgShape")->fitTo(*(test.data("zgdata")),
				 RooFit::ConditionalObservables(*(test.var("dMzg"))));
  
  TCanvas canv("test","test",600,600);

  RooPlot* frame = test.var("Mzg")->frame(90,200,150);
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle(Form("M_{%s} (GeV)",subscriptzg));
  frame->GetYaxis()->SetTitle("Entries");
  test.data("zgdata")->plotOn(frame);
  test.pdf("MzgBkgShape")->plotOn(frame,
				  RooFit::ProjWData(*(test.var("dMzg")),
						    *(test.data("zgdata"))));

  frame->Draw();
  canv.Print(Form("bkgshapetest_test_zg_%s.pdf",channel));
  canv.SetLogy();
  canv.Print(Form("bkgshapetest_test_zg_log_%s.pdf",channel));
  canv.Clear();

  delete frame;  
}

