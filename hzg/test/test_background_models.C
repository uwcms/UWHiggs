void test_background_models() {

  gROOT->ProcessLine(".X CMSStyle.C");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("libUWHiggshzg.so");

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

  int category = 4;

  char channel[] = Form("muon_cat%i",category);
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

  RooDataSet d( "data","M_{Z Jets} with Errors", &dllg, *wobs,
		Form(" Mz + Mzg > 185 && r94cat == %i ",category) ,"weight" );

  RooDataSet dold( "data_old","M_{Z Jets} with Errors", &dllg, *wobs,
		Form(" r94cat == %i ",category) ,"weight" );

  test.import(d);
  test.import(dold);
  
  test.factory("RooGaussModel::MzgResoShape(Mzg,bias[120,90,150],sigma[1,0.01,10])");
  test.factory("RooDecay::MzgBkgShape(Mzg,tau[5,0,50],MzgResoShape,RooDecay::SingleSided)");

  /*
  test.factory("RooLandau::MzgTruthShape(Mzg,lanMPV[130,100,180],lanWidth[10,0.01,50])");
  test.factory("RooGaussian::MzgResoShape(Mzg,0,sigma[0.01,10])");
  test.factory("FCONV::MzgBkgShape(Mzg,MzgTruthShape,MzgResoShape)");
  */

  test.pdf("MzgBkgShape")->fitTo(*(test.data("data")),
				 RooFit::SumW2Error(kTRUE));

  
  if(category == 1) {
    test.factory("RooBernstein::MzgBkgShapeOldPoly(Mzg,{c0_old[-1e-6,1.01],c1_old[-1e-6,1.01],c2_old[-1e-6,1.01],c3_old[-1e-6,1.01]})");   
    test.factory("RooStepBernstein::MzgBkgShapePolyTruth(Mzg,stepVal[0.1,0,1],{1.0,c1[0.5,-1e-6,1],c2[0.5,-1e-6,1],c3[0.5,-1e-6,1]})");    
  } else {
    test.factory("RooBernstein::MzgBkgShapeOldPoly(Mzg,{c0_old[-1e-6,1.01],c1_old[-1e-6,1.01],c2_old[-1e-6,1.01],c3_old[-1e-6,1.01],c4_old[-1e-6,1.01]})");    
    test.factory("RooStepBernstein::MzgBkgShapePolyTruth(Mzg,stepVal[0.1,0,1],{1.0,c1[0.5,-1e-6,1],c2[0.5,-1e-6,1],c3[0.5,-1e-6,1],c4[0.5,-1e-6,1]})");    
  }
  test.factory("RooGaussian::MzgResoShapePoly(Mzg,biasPoly[0],sigmaPoly[5,0.01,20])");
  test.factory("FCONV::MzgBkgShapePoly(Mzg,MzgBkgShapePolyTruth,MzgResoShapePoly)");
  test.var("Mzg")->setBins(15000,"cache");
  test.var("Mzg")->setRange("oldpolyfit",115,180);
  //test.var("Mzg")->setRange("NormalizationRangeForoldpolyfit",115,180);

  RooFitResult* expFitRes = test.pdf("MzgBkgShapePoly")->fitTo(*(test.data("data")),
				     RooFit::Minimizer("Minuit","simplex"),
				     RooFit::SumW2Error(kTRUE),
				     RooFit::Save(kTRUE));
  RooFitResult* polyFitRes = test.pdf("MzgBkgShapePoly")->fitTo(*(test.data("data")),
								RooFit::SumW2Error(kTRUE),
								RooFit::Save(kTRUE));
  
  /*
  test.pdf("MzgBkgShapeOldPoly")->fitTo(*(test.data("data")),
					RooFit::Minimizer("Minuit","simplex"),
					RooFit::Range("oldpolyfit"),
					RooFit::SumW2Error(kTRUE));
  test.pdf("MzgBkgShapeOldPoly")->fitTo(*(test.data("data")),
					RooFit::Range("oldpolyfit"),
					RooFit::SumW2Error(kTRUE));
  */
  

  RooRealIntegral* testint = test.pdf("MzgBkgShapePolyTruth")->createIntegral(RooArgSet(*test.var("Mzg")),
									      RooArgSet(*test.var("Mzg")));

  
  
  std::cout << "Normalized integral of step bernstein:: "<< testint->getVal() << std::endl;
  delete testint;

  TCanvas canv("test","test",600,600);

  //setup caching for making the plot
  test.var("Mzg")->setBins(100000,"cache");

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
  Double_t expchi2 = frame->chiSquare();
  /*
  test.pdf("MzgBkgShapePoly")->plotOn(frame,
				      RooFit::ProjWData(*(test.var("dMzg")),
							*(test.data("data"))),
				      RooFit::FillColor(kGreen),
				      RooFit::VisualizeError(*polyFitRes,2.0,kTRUE));
  test.pdf("MzgBkgShapePoly")->plotOn(frame,
				      RooFit::ProjWData(*(test.var("dMzg")),
							*(test.data("data"))),
				      RooFit::FillColor(kYellow),
				      RooFit::VisualizeError(*polyFitRes,1.0,kTRUE));
  */
  test.pdf("MzgBkgShapePoly")->plotOn(frame,
				      RooFit::ProjWData(*(test.var("dMzg")),
							*(test.data("data"))),
				      RooFit::LineColor(kRed));
  
  /*
  Double_t old_norm= test.data("data")->sumEntries("Mzg > 115 && Mzg < 180");
  Double_t tot_norm= test.data("data")->sumEntries();
  test.pdf("MzgBkgShapeOldPoly")->plotOn(frame,	
					 RooFit::Range("oldpolyfit"),
					 RooFit::LineColor(kGreen-3));
  */

  
  frame->GetYaxis()->SetRangeUser(0,70);
  frame->Draw();
  tlx.DrawLatex(0.60,0.88,Form("#chi^{2}_{Pol} = %.3f",frame->chiSquare()));
  tlx.DrawLatex(0.59,0.77,Form("#chi^{2}_{Exp} = %.3f",expchi2));
  tlx.DrawLatex(0.59,0.65,Form("cat = %i",category));
  
  canv.Print(Form("bkgshapetest_test_zg_%s.pdf",channel));  
  canv.Clear();
  canv.SetLogy();
  frame->GetYaxis()->SetRangeUser(1e-3,70);
  frame->Draw();
  tlx.DrawLatex(0.60,0.88,Form("#chi^{2}_{Pol} = %.3f",frame->chiSquare()));
  tlx.DrawLatex(0.59,0.77,Form("#chi^{2}_{Exp} = %.3f",expchi2));
  canv.Print(Form("bkgshapetest_test_zg_log_%s.pdf",channel));
  canv.Clear();

  delete frame;  

  f->Close();
  f2->Close();
}

