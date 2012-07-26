using namespace RooFit;

void build2dPdf(char* sigTemplateFile="../datafiles/Dsignal_3signals.root",
		 char* bkgTemplateFile="../datafiles/Dbackground_3signals.root"){

  gSystem->Load("../PDFs/RooMzzModelBkg_cc.so");

  RooRealVar nasaLD("nasaLD","nasaLD",0,1);
  RooRealVar zzmass("zzmass","zzmass",100,180);

  RooRealVar* a0;
  RooRealVar* a1;
  RooRealVar* a2;
  RooRealVar* a3;
  RooRealVar* a4;
  RooRealVar* a5;
  RooRealVar* a6;
  RooRealVar* a7;
  RooRealVar* a8;
  RooRealVar* a9;

  RooMzzModelBkg* mzzBkgProj;
  TFile* bkgTempFile;
  TH2F*  bkgTemplate;
  RooDataHist* bkgTempDataHist;
  RooHistPdf* bkgTemplatePdf;

  RooRealVar* mean;
  RooRealVar* sigma;
  RooRealVar* mean2;
  RooRealVar* sigma2;
  RooRealVar* frac;

  RooGaussian* gauss1;
  RooGaussian* gauss2;
  RooAddPdf* mzzSigProj;
  TFile* sigTempFile;
  TH2F*  sigTemplate;
  RooDataHist* sigTempDataHist;
  RooHistPdf* sigTemplatePdf;

  RooProdPdf* sigPDF;
  RooProdPdf* bkgPDF;

  frac=new RooRealVar("frac","frac",      6.28560e-01  ,0.0,1.0);
  mean=new RooRealVar("mean","mean",      1.24853e+02  ,110.,180.);
  mean2=new RooRealVar("mean2","mean2",   1.22672e+02  ,110.,180.);
  sigma=new RooRealVar("sigma","sigma",   1.40047e+00  ,0.,100.);
  sigma2=new RooRealVar("sigma2","sigma2",4.02136e+00  ,0.,100.);

  a0 = new RooRealVar("a0","a0",1.11879e+02 ,100,150);
  a1 = new RooRealVar("a1","a1",1.37619e+01 ,5,150);
  a2 = new RooRealVar("a2","a2",6.15976e+01 ,50,1000);
  a3 = new RooRealVar("a3","a3",6.65936e-02 ,0.05,0.5);
  a4 = new RooRealVar("a4","a4",1.94269e+02 ,150,200);
  a5 = new RooRealVar("a5","a5",1.50000e+01 ,5,15);
  a6 = new RooRealVar("a6","a6",2.50002e+01 ,25,45);
  a7 = new RooRealVar("a7","a7",3.27125e-01 ,-20,20);
  a8 = new RooRealVar("a8","a8",9.99441e+01 ,25,100);
  a9 = new RooRealVar("a9","a9",2.75410e-01 ,-10,10);

  gauss1 = new RooGaussian("gauss1","gauss2",zzmass,*mean,*sigma);
  gauss2 = new RooGaussian("gauss2","gauss2",zzmass,*mean2,*sigma2);

  mzzBkgProj = new RooMzzModelBkg("mzzBkgProj","mzzBkgProj",zzmass,*a0,*a1,*a2,*a3,*a4,*a5,*a6,*a7,*a8,*a9);
  mzzSigProj = new RooAddPdf("mzzSigProj","mzzSigProj",*gauss1,*gauss2,*frac);

  sigTempFile = new TFile(sigTemplateFile);
  sigTemplate = (TH2F*) sigTempFile->Get("h_mzzD");
  sigTempDataHist = new RooDataHist("sigTempDataHist","sigTempDataHist",RooArgSet(zzmass,nasaLD),sigTemplate);
  sigTemplatePdf = new RooHistPdf("sigTemplatePdf","sigTemplatePdf",RooArgSet(zzmass,nasaLD),*sigTempDataHist);

  bkgTempFile = new TFile(bkgTemplateFile);
  bkgTemplate = (TH2F*) bkgTempFile->Get("h_mzzD");
  bkgTempDataHist = new RooDataHist("bkgTempDataHist","bkgTempDataHist",RooArgSet(zzmass,nasaLD),bkgTemplate);
  bkgTemplatePdf = new RooHistPdf("bkgTemplatePdf","bkgTemplatePdf",RooArgSet(zzmass,nasaLD),*bkgTempDataHist);

  bkgPDF = new RooProdPdf("bkgPDF","bkgPDF",*bkgTemplatePdf,*mzzBkgProj);
  sigPDF = new RooProdPdf("sigPDF","sigPDF",*sigTemplatePdf,*mzzSigProj);

  cout << "sig data" << endl;

  TFile* dataFileSig = new TFile("/scratch0/hep/whitbeck/4lHelicity/datafiles/7TeV/testBuildModel/SMHiggs_125_JHU_v3_wResolution_withDiscriminants.root");
  TTree* dataTreeSig = (TTree*) dataFileSig->Get("angles");
  RooDataSet* dataSig = new RooDataSet("dataSig","dataSig",dataTreeSig,RooArgSet(zzmass,nasaLD));

  cout << "bkg data" << endl;

  TFile* dataFileBkg = new TFile("/scratch0/hep/whitbeck/4lHelicity/datafiles/7TeV/training/EWKZZ4l_Powheg_mZ4GeV_total_v30_wResolution_withDiscriminants.root");
  TTree* dataTreeBkg = (TTree*) dataFileBkg->Get("angles");
  RooDataSet* dataBkg = new RooDataSet("dataBkg","dataBkg",dataTreeBkg,RooArgSet(zzmass,nasaLD));

  cout << "rooplots" << endl;

  RooPlot* mzzPlot = zzmass.frame(20);
  RooPlot* dPlot = nasaLD.frame(20);
  dataBkg->plotOn(mzzPlot,Rescale((double)1./dataBkg->sumEntries()),MarkerColor(2));
  bkgPDF->plotOn(mzzPlot,Normalization((double)1./dataBkg->sumEntries()),LineColor(2));
  dataSig->plotOn(mzzPlot,Rescale((double)1./dataSig->sumEntries()),MarkerColor(4));
  sigPDF->plotOn(mzzPlot,Normalization((double)1./dataSig->sumEntries()),LineColor(4));
  mzzPlot->GetYaxis()->SetRangeUser(0,1);

  //dataBkg->plotOn(dPlot,Rescale((double)1./dataBkg->sumEntries()),MarkerColor(2));
  bkgPDF->plotOn(dPlot,Normalization((double)1./dataBkg->sumEntries()),LineColor(2));
  //dataSig->plotOn(dPlot,Rescale((double)1./dataSig->sumEntries()),MarkerColor(4));
  sigPDF->plotOn(dPlot,Normalization((double)1./dataSig->sumEntries()),LineColor(4));
  dPlot->GetYaxis()->SetRangeUser(0,.3);

  TCanvas* mzzCan = new TCanvas("mzzCan","mzz",400,400);
  mzzPlot->Draw();
  TCanvas* dCan = new TCanvas("dCan","D",400,400);
  dPlot->Draw();

  
}
