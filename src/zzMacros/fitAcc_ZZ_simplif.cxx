// fit for graviton acceptance

// C++ includes
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooBinning.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooRandom.h"

// #include "RooJetsPentaSpinTwo.h"
#include "RooPentaSpinTwo.h"

using namespace RooFit;

double mZ = 91.1876;

double fpm( double mG ) {

  double apm = mG*mG;
  double app = 4*mZ*mZ/sqrt(6.);
  double ap0 = 2*mG*mZ/sqrt(2.);
  double a00 = (mG*mG*mG*mG/(mZ*mZ*sqrt(6.)))*(mZ*mZ/(mG*mG) + 4*mZ*mZ*mZ*mZ/(mG*mG*mG*mG));
  return apm*apm/(2*apm*apm + 2*app*app + 4*ap0*ap0 + a00*a00);

}

double fpp( double mG ) {

  double apm = mG*mG;
  double app = 4*mZ*mZ/sqrt(6.);
  double ap0 = 2*mG*mZ/sqrt(2.);
  double a00 = (mG*mG*mG*mG/(mZ*mZ*sqrt(6.)))*(mZ*mZ/(mG*mG) + 4*mZ*mZ*mZ*mZ/(mG*mG*mG*mG));
  return app*app/(2*apm*apm + 2*app*app + 4*ap0*ap0 + a00*a00); 

}

double fp0( double mG ) {

  double apm = mG*mG;
  double app = 4*mZ*mZ/sqrt(6.);
  double ap0 = 2*mG*mZ/sqrt(2.);
  double a00 = (mG*mG*mG*mG/(mZ*mZ*sqrt(6.)))*(mZ*mZ/(mG*mG) + 4*mZ*mZ*mZ*mZ/(mG*mG*mG*mG));
  return ap0*ap0/(2*apm*apm + 2*app*app + 4*ap0*ap0 + a00*a00); 

}

double fz1( double mG ) {

  ifstream nhanFile("./fqqTable.txt");
  double theMass;
  double qqFrac;
  double rightqqFrac;

  while (nhanFile >> theMass >> qqFrac) {
    if (fabs(theMass - mG) < 10.) rightqqFrac = qqFrac;
  }
  
  return rightqqFrac;
}

int main(int argc, char* argv[]) {

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");	

  char *filename;
  float mass = 0.; 
  float minAcc = 0.;
  float maxcos = 1.0;
  float maxcos2 = 1.0;
  float maxcosS = 1.0;

  for(Int_t i=1;i<argc;i++){
    char *pchar = argv[i];

    switch(pchar[0]){

    case '-':{

      switch(pchar[1]){

      case 'f':
        filename = argv[i+1];
        cout << "File name for fitted data is " << filename << endl;
        break;

      case 'm':
	mass = atof(argv[i+1]);
	cout << "Graviton mass is " << mass << " GeV" << endl;
	break;

      case 'a':
	minAcc = atof(argv[i+1]);
	cout << "Minimum acceptance considered is " << minAcc << endl;
	break;	
	
      case 'c':
	
	switch(pchar[2]){
	  
	  case '1':
	    maxcos = atof(argv[i+1]);
	    cout << "Maximum cosine of helicity angle 1 is " << maxcos << endl;
	    break;
       
          case '2':
	    maxcos2 = atof(argv[i+1]);
	    cout << "Maximum cosine of helicity angle 2 is " << maxcos2 << endl;
	    break;

	  case 's':
	    maxcosS = atof(argv[i+1]);
	    cout << "Maximum cosine of theta* is " << maxcosS << endl;
	    break;  
	}

      } 

    }
    }
  }

  RooWorkspace *ws = new RooWorkspace("ws");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *data = (RooDataSet*)fIn.Get("data");
  ws->import(*data);
  
  RooRealVar fppVal("fppVal","fppVal",fpp(mass));
  RooRealVar fpmVal("fpmVal","fpmVal",fpm(mass));
  RooRealVar fp0Val("fp0Val","fp0Val",fp0(mass));
  RooRealVar fz1Val("fz1Val","fz1Val",fz1(mass));
  RooFormulaVar fz2Val("fz2Val","1-@0",RooArgList(fz1Val));
  RooRealVar R1Val("R1Val","R1Val",0.15);
  RooRealVar R2Val("R2Val","R2Val",0.);

  RooRealVar zero("zero","zero",0.);
  RooRealVar one("one","one",1.);

  ws->var("cosTheta1")->setMin(-maxcos);
  ws->var("cosTheta1")->setMax(maxcos);
  // ws->var("cosTheta2")->setMin(0.05);
  ws->var("cosTheta2")->setMax(maxcos2);
  ws->var("cosThetaStar")->setMin(-maxcosS);
  ws->var("cosThetaStar")->setMax(maxcosS);

  /* RooJetsPentaSpinTwo fitf("fitf", "A complicated function",
			   *ws->var("cosTheta1"),
			   *ws->var("cosTheta2"),
			   *ws->var("Phi"),
			   *ws->var("cosThetaStar"),
			   *ws->var("Phi1"),
			   fppVal, fppVal, fpmVal, fp0Val, fp0Val,
			   zero, zero, zero, zero, zero,
			   fz1Val, fz2Val, R1Val, R2Val,
			   zero, zero, zero, zero,  // para2, para4 para6, para8
			   one, zero, zero, zero,  // acca0, acca1, acca2, acca4
			   zero, zero, zero, zero,  // a2, a4, b2, b4,  
			   zero, zero, one);  // cutOff, g, N
  */
  
  RooPentaSpinTwo fitf("fitf", "A complicated function",
		       *ws->var("cosTheta1"),
		       *ws->var("cosTheta2"),
		       *ws->var("Phi"),
		       *ws->var("cosThetaStar"),
		       *ws->var("Phi1"),
		       fppVal, fppVal, fpmVal, fp0Val, fp0Val,
		       zero, zero, zero, zero, zero,
		       fz1Val, fz2Val, R1Val, R2Val,
		       zero, zero, zero, zero,  // para2, para4 para6, para8
		       one, zero, zero, zero,  // acca0, acca1, acca2, acca4
		       zero, zero, zero, zero,  // a2, a4, cutOff, g 
		       zero, zero, one);  // b2, b4, N

  cout << "***" << endl << "GENERATING" << endl << "***" << endl;
  RooRandom::randomGenerator()->SetSeed(4.);
  RooDataSet *uncorr = fitf.generate(RooArgSet(*ws->var("cosTheta1"),
					       *ws->var("cosTheta2"),
					       *ws->var("Phi"),
					       *ws->var("cosThetaStar"),
					       *ws->var("Phi1")), 
				     NumEvents(100000),Name("uncorr"));

  TCanvas c1("c1","c1",10,10,1200,800);
  c1.Divide(3,2);

  cout << "***" << endl << "TEST PLOTTING" << endl << "***" << endl;
  c1.cd(1);
  RooPlot *t1frame = ws->var("cosTheta1")->frame();
  uncorr->plotOn(t1frame,DataError(RooAbsData::SumW2));
  t1frame->Draw();
  
  c1.cd(2);
  RooPlot *t2frame = ws->var("cosTheta2")->frame();
  uncorr->plotOn(t2frame,DataError(RooAbsData::SumW2));
  t2frame->Draw(); 

  c1.cd(3);
  RooPlot *tStarframe = ws->var("cosThetaStar")->frame();
  uncorr->plotOn(tStarframe,DataError(RooAbsData::SumW2));
  tStarframe->Draw();

  c1.cd(4);
  RooPlot *pframe = ws->var("Phi")->frame();
  uncorr->plotOn(pframe,DataError(RooAbsData::SumW2));
  pframe->Draw();

  c1.cd(5);
  RooPlot *p1frame = ws->var("Phi1")->frame();
  uncorr->plotOn(p1frame,DataError(RooAbsData::SumW2));
  p1frame->Draw(); 

  char fileout[200];
  sprintf(fileout,"genResult%d.ps",int(mass));
  c1.SaveAs(fileout);
  
  TH1D* acc_cosTheta1 = new TH1D("acc_cosTheta1","acc_cosTheta1",50,ws->var("cosTheta1")->getMin(),ws->var("cosTheta1")->getMax());  
  TH1D* acc_cosTheta2 = new TH1D("acc_cosTheta2","acc_cosTheta2",50,ws->var("cosTheta2")->getMin(),ws->var("cosTheta2")->getMax());    
  TH1D* acc_cosThetaStar = new TH1D("acc_cosThetaStar","acc_cosThetaStar",50,ws->var("cosThetaStar")->getMin(),ws->var("cosThetaStar")->getMax());  
  TH1D* acc_Phi = new TH1D("acc_Phi","acc_Phi",50,ws->var("Phi")->getMin(),ws->var("Phi")->getMax());  
  TH1D* acc_Phi1 = new TH1D("acc_Phi1","acc_Phi1",50,ws->var("Phi1")->getMin(),ws->var("Phi1")->getMax()); 

  acc_cosTheta1->Sumw2();
  acc_cosTheta2->Sumw2();
  acc_cosThetaStar->Sumw2();
  acc_Phi->Sumw2();
  acc_Phi1->Sumw2();

  data->fillHistogram(acc_cosTheta1,RooArgList(*ws->var("cosTheta1")));
  data->fillHistogram(acc_cosTheta2,RooArgList(*ws->var("cosTheta2")));
  data->fillHistogram(acc_cosThetaStar,RooArgList(*ws->var("cosThetaStar")));
  data->fillHistogram(acc_Phi,RooArgList(*ws->var("Phi")));
  data->fillHistogram(acc_Phi1,RooArgList(*ws->var("Phi1")));

  TH1D* den_cosTheta1 = new TH1D("den_cosTheta1","den_cosTheta1",50,ws->var("cosTheta1")->getMin(),ws->var("cosTheta1")->getMax());  
  TH1D* den_cosTheta2 = new TH1D("den_cosTheta2","den_cosTheta2",50,ws->var("cosTheta2")->getMin(),ws->var("cosTheta2")->getMax());    
  TH1D* den_cosThetaStar = new TH1D("den_cosThetaStar","den_cosThetaStar",50,ws->var("cosThetaStar")->getMin(),ws->var("cosThetaStar")->getMax());  
  TH1D* den_Phi = new TH1D("den_Phi","den_Phi",50,ws->var("Phi")->getMin(),ws->var("Phi")->getMax());  
  TH1D* den_Phi1 = new TH1D("den_Phi1","den_Phi1",50,ws->var("Phi1")->getMin(),ws->var("Phi1")->getMax()); 

  den_cosTheta1->Sumw2();
  den_cosTheta2->Sumw2();
  den_cosThetaStar->Sumw2();
  den_Phi->Sumw2();
  den_Phi1->Sumw2();

  uncorr->fillHistogram(den_cosTheta1,RooArgList(*ws->var("cosTheta1")));
  uncorr->fillHistogram(den_cosTheta2,RooArgList(*ws->var("cosTheta2")));
  uncorr->fillHistogram(den_cosThetaStar,RooArgList(*ws->var("cosThetaStar")));
  uncorr->fillHistogram(den_Phi,RooArgList(*ws->var("Phi")));
  uncorr->fillHistogram(den_Phi1,RooArgList(*ws->var("Phi1")));

  acc_cosTheta1->Divide(den_cosTheta1);
  acc_cosTheta2->Divide(den_cosTheta2);
  acc_cosThetaStar->Divide(den_cosThetaStar);
  acc_Phi->Divide(den_Phi);
  acc_Phi1->Divide(den_Phi1);

  TCanvas c2("c2","c2",10,10,1200,800);
  c2.Divide(3,2);

  cout << "***" << endl << "ACCEPTANCE FITTING" << endl << "***" << endl;

  c2.cd(1);
  acc_cosTheta1->SetMinimum(0.);
  TF1* fcosTheta1 = new TF1("fcosTheta1","pol4",ws->var("cosTheta1")->getMin(),ws->var("cosTheta1")->getMax());
  fcosTheta1->SetParameters(0.5,0.,-0.4,0.,-0.15); 
  fcosTheta1->SetLineWidth(2);
  fcosTheta1->SetLineColor(kBlue); 
  acc_cosTheta1->Fit(fcosTheta1); 
  // acc_cosTheta1->Draw();
  
  c2.cd(2);
  acc_cosTheta2->SetMinimum(0.);
  // TF1* fcosTheta2 = new TF1("fcosTheta2","pol4",ws->var("cosTheta2")->getMin(),ws->var("cosTheta2")->getMax());
  TF1* fcosTheta2 = new TF1("fcosTheta2","([0]+[2]*x+[4]*x*x)/(1.+exp((x-[1])*[3]))",ws->var("cosTheta2")->getMin(),ws->var("cosTheta2")->getMax());
  fcosTheta2->SetParameters(0.37,0.79,0.4,32.,-0.5);
  fcosTheta2->SetLineWidth(2);
  fcosTheta2->SetLineColor(kBlue); 
  acc_cosTheta2->Fit(fcosTheta2); 
  // acc_cosTheta2->Draw();

  c2.cd(3);
  acc_cosThetaStar->SetMinimum(0.);
  TF1* fcosThetaStar;
  if (mass > 900.) {
    fcosThetaStar = new TF1("fcosThetaStar","pol6",ws->var("cosThetaStar")->getMin(),ws->var("cosThetaStar")->getMax());
    // fcosThetaStar->FixParameter(2,12.);
  }
  else if (mass > 450.) fcosThetaStar = new TF1("fcosThetaStar","pol4",ws->var("cosThetaStar")->getMin(),ws->var("cosThetaStar")->getMax());
  else fcosThetaStar = new TF1("fcosThetaStar","pol2",ws->var("cosThetaStar")->getMin(),ws->var("cosThetaStar")->getMax());
  // fcosThetaStar->FixParameter(1,0.);
  fcosThetaStar->SetLineWidth(2);
  fcosThetaStar->SetLineColor(kBlue);
  acc_cosThetaStar->Fit(fcosThetaStar); 
  // acc_cosThetaStar->Draw();

  c2.cd(4);
  acc_Phi->SetMinimum(0.);
  TF1* fPhi = new TF1("fPhi","pol0",ws->var("Phi")->getMin(),ws->var("Phi")->getMax());
  // fPhi->FixParameter(1,0.);
  fPhi->SetLineWidth(2);
  fPhi->SetLineColor(kBlue);
  acc_Phi->Fit(fPhi);
  // acc_Phi->Draw();

  c2.cd(5);
  acc_Phi1->SetMinimum(0.);
  TF1* fPhi1 = new TF1("fPhi1","[0]*(1.-[1]*cos(2*x))",ws->var("Phi1")->getMin(),ws->var("Phi1")->getMax());
  // fPhi1->FixParameter(1,0.);
  fPhi1->SetLineWidth(2);
  fPhi1->SetLineColor(kBlue);
  acc_Phi1->Fit(fPhi1);
  // acc_Phi1->Draw();
  
  sprintf(fileout,"accResult%d.ps",int(mass));
  c2.SaveAs(fileout);

  // AVOID NEGATIVE NUMBERS
  // Find at which value the acceptance is = x% of the maximum using Newton's method

  float cosMax = 1.0;
  float cosMin = 0.0;
  float cosTheta1Lim = -1.0;
  float cosTheta2Lim = -1.0;
  float cosThetaStarLim = -1.0;
  for (Int_t iN = 0; iN <= 5000; iN++) {
    float theInterval = (cosMax-cosMin)/5000.;
    float thePoint = cosMin + iN*theInterval;  
    if (fcosTheta1->Eval(thePoint) > minAcc*fcosTheta1->GetParameter(0) && fcosTheta1->Eval(thePoint+theInterval) < minAcc*fcosTheta1->GetParameter(0)) cosTheta1Lim = thePoint;
    if (fcosTheta2->Eval(thePoint) > minAcc*fcosTheta2->GetParameter(0) && fcosTheta2->Eval(thePoint+theInterval) < minAcc*fcosTheta2->GetParameter(0)) cosTheta2Lim = thePoint;
    if (fcosThetaStar->Eval(thePoint) > minAcc*fcosThetaStar->GetParameter(0) && fcosThetaStar->Eval(thePoint+theInterval) < minAcc*fcosThetaStar->GetParameter(0)) cosThetaStarLim = thePoint;
  }

  cout << "***" << endl << "PARAMETERS FOR COMPLETE FIT" << endl << "***" << endl;

  float para2red = fcosThetaStar->GetParameter(2)/fcosThetaStar->GetParameter(0);
  float para2rederr = fcosThetaStar->GetParError(2)/fcosThetaStar->GetParameter(0);
  float para4red = 0.;
  float para6red = 0.;
  float para8red = 0.;
  float para4rederr, para6rederr, para8rederr;

  if (mass > 450.) {
    para4red = fcosThetaStar->GetParameter(4)/fcosThetaStar->GetParameter(0);
    para4rederr = fcosThetaStar->GetParError(4)/fcosThetaStar->GetParameter(0);
    if (mass > 900.) {
      para6red = fcosThetaStar->GetParameter(6)/fcosThetaStar->GetParameter(0);
      para6rederr = fcosThetaStar->GetParError(6)/fcosThetaStar->GetParameter(0);
      // para8red = fcosThetaStar->GetParameter(8)/fcosThetaStar->GetParameter(0);
      // para8rederr = fcosThetaStar->GetParError(8)/fcosThetaStar->GetParameter(0);
     }
  }

  float cutOffred = fcosTheta2->GetParameter(1) /* /fcosTheta2->GetParameter(0) */;
  float cutOffrederr = fcosTheta2->GetParError(1)/* /fcosTheta2->GetParameter(0)*/;
  float a2red = fcosTheta2->GetParameter(2)/fcosTheta2->GetParameter(0);
  float a2rederr = fcosTheta2->GetParError(2)/fcosTheta2->GetParameter(0);
  float gred = fcosTheta2->GetParameter(3) /* /fcosTheta2->GetParameter(0) */;
  float grederr = fcosTheta2->GetParError(3) /* /fcosTheta2->GetParameter(0) */;
  float a4red = fcosTheta2->GetParameter(4)/fcosTheta2->GetParameter(0);
  float a4rederr = fcosTheta2->GetParError(4)/fcosTheta2->GetParameter(0);

  float b2red = fcosTheta1->GetParameter(2)/fcosTheta1->GetParameter(0);
  float b2rederr = fcosTheta1->GetParError(2)/fcosTheta1->GetParameter(0);
  float b4red = fcosTheta1->GetParameter(4)/fcosTheta1->GetParameter(0);
  float b4rederr = fcosTheta1->GetParError(4)/fcosTheta1->GetParameter(0);

  cout << "para2 = " << para2red << endl;
  if (mass > 450.) {
    cout << "para4 = " << para4red << endl;
    if (mass > 900.) {
      cout << "para6 = " << para6red << endl;
      // cout << "para8 = " << para8red << endl;
    }
  }
  if (cosThetaStarLim >= 0.) cout << "cosThetaStarLim = " << cosThetaStarLim << endl << endl; 
  else cout << "No cosThetaStarLim" << endl << endl;

  cout << "acca0 = " << fPhi1->GetParameter(0) << endl;
  cout << "acca2 = " << -(fPhi1->GetParameter(1)) << endl << endl;
  
  cout << "cutOff = " << cutOffred << endl;
  // cout << "cutOff = " << fcosTheta2->GetParameter(1) << endl;
  cout << "a2 = " << a2red << endl;
  cout << "g = " << gred << endl;
  // cout << "g = " << fcosTheta2->GetParameter(3) << endl;
  cout << "a4 = " << a4red << endl;
  if (cosTheta2Lim >= 0.) cout << "cosTheta2Lim = " << cosTheta2Lim << endl << endl; 
  else cout << "No cosTheta2Lim" << endl << endl; 

  cout << "b2 = " << b2red << endl;
  cout << "b4 = " << b4red << endl;
  if (cosTheta1Lim >= 0.) cout << "cosTheta1Lim = " << cosTheta1Lim << endl << endl; 
  else cout << "No cosTheta1Lim" << endl << endl;

  // DUMP FOR ROOFIT
  sprintf(fileout,"acc_pars%d.txt",int(mass));
  ofstream of(fileout);
  of << "[Acceptance] \n";
  of << "para2 = " << para2red << " +/- " << para2rederr << " L(-5 - 5) \n";
  if (mass > 900.) {
    of << "para4 = " << para4red << " +/- " << para4rederr << " L(-20 - 20) \n";
    of << "para6 = " << para6red << " +/- " << para6rederr << " L(-20 - 10) \n";
    of << "para8 = " << para8red << " C \n";
  } else if (mass > 450.) {
    of << "para4 = " << para4red << " +/- " << para4rederr << " L(-20 - 20) \n";
    of << "para6 = " << para6red << " C \n";
    of << "para8 = " << para8red << " C \n";
  } else {
    of << "para4 = " << para4red << " C \n";
    of << "para6 = " << para6red << " C \n";
    of << "para8 = " << para8red << " C \n";
  }
  of << "acca2 = " << -(fPhi1->GetParameter(1)) << " +/- " << fPhi1->GetParError(1) << " L(-2 - 2) \n";
  of << "a2 = " << a2red << " +/- " << a2rederr << " L(" <<  a2red-(a2rederr/2.) << " - " << a2red+(a2rederr/2.) << ") \n";
  // of << "cutOff = " << cutOffred << " +/- " << cutOffrederr << " L(" <<  cutOffred-cutOffrederr << " - " << cutOffred+cutOffrederr << ") \n";
  of << "cutOff = " << cutOffred << " C \n";
  of << "g = " << gred << " +/- " << grederr << " L(" <<  gred-grederr << " - " << gred+grederr << ") \n";
  of << "a4 = " << a4red << " +/- " << a4rederr << " L(" <<  a4red-(a4rederr/2.) << " - " << a4red+(a4rederr/2.) << ") \n";
  // of << "b2 = " << b2red << " +/- " << b2rederr << " L(-5 - 5) \n";
  // of << "b4 = " << b4red << " +/- " << b4rederr << " L(-5 - 5) \n";
  of << "b2 = " << b2red << " C \n";
  of << "b4 = " << b4red << " C \n";
  of.close();

  return 1;
}
