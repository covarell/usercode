// fit for graviton acceptance

// C++ includes
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
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
  float maxcos = 1.0;
  float maxcos2 = 1.0;
  float maxcosS = 1.0;
  float minAcc = 0.;
  bool constrCostheta1 = false;
  bool constrCostheta2 = false;
  bool constrCosthetaStar = false;
  bool onlyDraw = false;

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

      case 'o':
	onlyDraw = true;
	cout << "Only draw, no fit" << endl;
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

	break;

      case 'u':
	
	switch(pchar[2]){
	  
	  case '1':
	    constrCostheta1 = true;
	    cout << "cos(theta1) acceptance is constrained!" << endl;
	    break;
       
          case '2':
	    constrCostheta2 = true;
	    cout << "cos(theta2) acceptance is constrained!" << endl;
	    break;

	  case 's':
	    constrCosthetaStar = true;
	    cout << "cos(theta*) acceptance is constrained!" << endl;
	    break;  
	}
	
	break;
      } 
    }
    }
  }
  

  RooWorkspace *ws = new RooWorkspace("ws");

  TFile fIn(filename);
  fIn.cd();
  RooDataSet *origdata = (RooDataSet*)fIn.Get("data");
  origdata->SetName("origdata");

  char theCut[200];
  sprintf(theCut,"cosTheta2 < %f && abs(cosTheta1) < %f && abs(cosThetaStar) < %f",maxcos2,maxcos,maxcosS);
  RooDataSet *data = (RooDataSet*)origdata->reduce(theCut);
  // ws->import(*data);
  
  // DO A COMPLETELY NEW RooDataSet (FOR CHECKING)
  RooRealVar* cosThetaStarNew = new RooRealVar("cosThetaStarNew","cos(#theta^{*})",-maxcosS,maxcosS);
  RooRealVar* cosTheta1New = new RooRealVar("cosTheta1New","cos(#theta_{1})",-maxcos,maxcos);	 
  RooRealVar* cosTheta2New = new RooRealVar("cosTheta2New","cos(#theta_{2})",0.001,maxcos2);	 
  RooRealVar* Phi1New = new RooRealVar("Phi1New","#Phi_{1}",-3.14,3.14); 
  RooRealVar* PhiNew = new RooRealVar("PhiNew","#Phi",-3.14,3.14);    

  RooArgList varlist(*cosThetaStarNew, *cosTheta1New, *cosTheta2New, *Phi1New, *PhiNew);
  RooDataSet *dataNew = new RooDataSet("dataNew","A sample",varlist);
  const RooArgSet* aRow;

  for (int ev=0; ev < data->numEntries(); ++ev) {

    aRow = data->get(ev);
    RooRealVar* thiscs = (RooRealVar*)aRow->find("cosThetaStar");
    cosThetaStarNew->setVal( thiscs->getVal() );
    RooRealVar* thisc1 = (RooRealVar*)aRow->find("cosTheta1");
    cosTheta1New->setVal( thisc1->getVal() );
    RooRealVar* thisc2 = (RooRealVar*)aRow->find("cosTheta2");
    cosTheta2New->setVal( thisc2->getVal() );
    RooRealVar* thisp1 = (RooRealVar*)aRow->find("Phi1");
    Phi1New->setVal( thisp1->getVal() );
    RooRealVar* thisp = (RooRealVar*)aRow->find("Phi");
    PhiNew->setVal( thisp->getVal() );

    dataNew->add(RooArgSet(*cosThetaStarNew, *cosTheta1New, *cosTheta2New, *Phi1New, *PhiNew));
 
  }
  ws->import(*dataNew);

  RooRealVar fppVal("fppVal","fppVal",fpp(mass));
  RooRealVar fpmVal("fpmVal","fpmVal",fpm(mass));
  RooRealVar fp0Val("fp0Val","fp0Val",fp0(mass));
  RooRealVar fz1Val("fz1Val","fz1Val",fz1(mass));
  RooFormulaVar fz2Val("fz2Val","1-@0",RooArgList(fz1Val));
  RooRealVar R1Val("R1Val","R1Val",0.15);
  RooRealVar R2Val("R2Val","R2Val",0.);

  RooRealVar zero("zero","zero",0.);
  RooRealVar one("one","one",1.);

  RooRealVar para4("para4","para4",-20.,20.);   ws->import(para4);
  RooRealVar para6("para6","para6",-20.,20.);   ws->import(para6);
  RooRealVar para8("para8","para8",-20.,20.);   ws->import(para8);
  RooRealVar acca0("acca0","acca0",0.1,100.); 
  RooRealVar acca1("acca1","acca1",-2.,2.);  
  RooRealVar acca2("acca2","acca2",-2.,2.); 
  RooRealVar acca4("acca4","acca4",-2.,2.);  
  RooRealVar a2("a2","a2",-10.,10.);   ws->import(a2);
  RooRealVar a4("a4","a4",-10.,10.);   ws->import(a4);
  RooRealVar g("g","g",5.,100.);   ws->import(g);   
  RooRealVar b2("b2","b2",-1.,1.);   ws->import(b2);

  RooRealVar minAccVar("minAccVar","minAccVar",minAcc);
  ws->import(minAccVar);

  if (constrCostheta1) {
    RooRealVar maxCosVar1("maxCosVar1","maxCosVar1",maxcos);
    ws->import(maxCosVar1);
    RooFormulaVar b4("b4","(-@0*pow(@1,2)-1.0+@2)/pow(@1,4)",RooArgList(*ws->var("b2"),*ws->var("maxCosVar1"),*ws->var("minAccVar")));
    ws->import(b4);
  } else {
    RooRealVar b4("b4","b4",-10.,10.);
    // b4.setVal(-0.8006);   // b4.setConstant(kTRUE);
    ws->import(b4);
  }
  
  if (constrCostheta2) {
    RooRealVar maxCosVar2("maxCosVar2","maxCosVar2",maxcos2);
    ws->import(maxCosVar2);
    RooFormulaVar cutOff("cutOff","(-@0*pow(@1,2)-@2*pow(@1,3)-@3*pow(@1,4)-1.0+@4)/@1",RooArgList(*ws->var("a2"),*ws->var("maxCosVar2"),*ws->var("g"),*ws->var("a4"),*ws->var("minAccVar")));
    ws->import(cutOff);
  } else {
    RooRealVar cutOff("cutOff","cutOff",0.5,1.);
    ws->import(cutOff);
  }

  if (constrCosthetaStar) {
    RooRealVar maxCosVarS("maxCosVarS","maxCosVarS",maxcosS);
    ws->import(maxCosVarS);
    RooFormulaVar para2("para2","(-1.0+@0-@1*pow(@2,4)-@3*pow(@2,6)-@4*pow(@2,8))/pow(@2,2)",RooArgList(*ws->var("minAccVar"),*ws->var("para4"),*ws->var("maxCosVarS"),*ws->var("para6"),*ws->var("para8")));
    ws->import(para2);
  } else {
    RooRealVar para2("para2","para2",-10.,10.);
    // para2.setVal(-0.6051);  // para2.setConstant(kTRUE);
    ws->import(para2);
  }

  RooArgSet thePars(RooArgList(acca2,*ws->var("para4"), *ws->var("para6"),
			       *ws->var("para8"),
			       *ws->var("a2"), *ws->var("a4"), 
			       *ws->var("g"), *ws->var("b2")));
  if (!constrCostheta1) thePars.add(*ws->var("b4"));
  if (!constrCostheta2) thePars.add(*ws->var("cutOff"));
  if (!constrCosthetaStar) thePars.add(*ws->var("para2"));

  char filein[200];
  if (onlyDraw) sprintf(filein,"acc_parsEstim%d.txt",int(mass));
  else sprintf(filein,"acc_pars%d.txt",int(mass));
  thePars.readFromFile(filein,0,"Acceptance");

  // para4.setVal(0.0);  para4.setConstant(kTRUE);
  acca0.setVal(1.0);  acca0.setConstant(kTRUE);
  acca1.setVal(0.0);  acca1.setConstant(kTRUE);
  // acca2.setVal(-0.2865);  // acca2.setConstant(kTRUE);
  acca4.setVal(0.0);  acca4.setConstant(kTRUE); 
  // ws->var("a2")->setVal(3.22617);  ws->var("a2")->setConstant(kTRUE);
  // ws->var("g")->setVal(-9.85547);  ws->var("g")->setConstant(kTRUE);
  // ws->var("a4")->setVal(6.00711);  // ws->var("a4")->setConstant(kTRUE);
  // ws->var("b2")->setVal(-0.2637);  // ws->var("b2")->setConstant(kTRUE);

  /* RooJetsPentaSpinTwo fitf("fitf", "A complicated function",
			   *ws->var("cosTheta1New"),   
			   *ws->var("cosTheta2New"),   
			   *ws->var("PhiNew"),	    
			   *ws->var("cosThetaStarNew"),
			   *ws->var("Phi1New"),        
			   fppVal, fppVal, fpmVal, fp0Val, fp0Val,
			   zero, zero, zero, zero, zero,
			   fz1Val, fz2Val, R1Val, R2Val,
			   *ws->function("para2"), para4, 
			   zero, zero, // para6, para8
			   acca0, acca1, acca2, acca4,
			   *ws->var("a2"), *ws->var("a4"), 
			   *ws->var("b2"), *ws->function("b4"),
			   *ws->function("cutOff"), *ws->var("g"), one);  // N
			   */

  /* RooPentaSpinTwo fitf("fitf", "A complicated function",
		       *ws->var("cosTheta1New"),   
		       *ws->var("cosTheta2New"),   
		       *ws->var("PhiNew"),	    
		       *ws->var("cosThetaStarNew"),
		       *ws->var("Phi1New"),        
		       fppVal, fppVal, fpmVal, fp0Val, fp0Val,
		       zero, zero, zero, zero, zero,
		       fz1Val, fz2Val, R1Val, R2Val,
		       *ws->function("para2"), *ws->var("para4"),
		       *ws->var("para6"), *ws->var("para8"),
		       acca0, acca1, acca2, acca4,
		       *ws->var("a2"), *ws->var("a4"), 
		       *ws->function("cutOff"), *ws->var("g"), 
		       *ws->var("b2"), *ws->function("b4"), one);  // N
		       */

  RooPentaSpinTwo fitf("fitf", "A complicated function",
		       *ws->var("cosTheta1New"),   
		       *ws->var("cosTheta2New"),   
		       *ws->var("PhiNew"),	    
		       *ws->var("cosThetaStarNew"),
		       *ws->var("Phi1New"),        
		       fppVal, fppVal, fpmVal, fp0Val, fp0Val,
		       zero, zero, zero, zero, zero,
		       fz1Val, fz2Val, R1Val, R2Val,
		       *ws->var("para2"), *ws->var("para4"),
		       *ws->var("para6"), *ws->var("para8"),
		       acca0, acca1, acca2, acca4,
		       *ws->var("a2"), *ws->var("a4"), 
		       *ws->var("cutOff"), *ws->var("g"), 
		       *ws->var("b2"), *ws->var("b4"), one);  // N

  // MAXIMUM LIKELIHOOD FIT
  if (!onlyDraw) {
    cout << "***" << endl << "FITTING" << endl << "***" << endl;
    fitf.fitTo(*dataNew,Minos(0),Hesse(0),SumW2Error(kTRUE),NumCPU(4),Save(1));
  }

  char fileout[200];
  TCanvas c1("c1","c1",10,10,1200,800);
  c1.Divide(3,2);

  // DUMP FOR ROOT
  if (!onlyDraw) {
    
    cout << "***" << endl << "WRITING" << endl << "***" << endl;
    sprintf(fileout,"acc_parsFinal%d.txt",int(mass));
    ofstream of(fileout);
    /* of << "para2 " << ws->var("para2")->getVal() << " " << ws->var("para2")->getError() << endl;
    of << "para4 " << ws->var("para4")->getVal() << " " << ws->var("para4")->getError() << endl;
    of << "para6 " << ws->var("para6")->getVal() << " " << ws->var("para6")->getError() << endl;
    of << "para8 " << ws->var("para8")->getVal() << " " << ws->var("para8")->getError()  << endl;
    of << "acca2 " << acca2.getVal() << " " << acca2.getError() << endl;
    of << "a2 " << ws->var("a2")->getVal() << " " << ws->var("a2")->getError() << endl;
    of << "cutOff " << ws->var("cutOff")->getVal() << " 0.01" << endl;
    of << "g " << ws->var("g")->getVal() << " " << ws->var("g")->getError() << endl;
    of << "a4 " << ws->var("a4")->getVal() << " " << ws->var("a4")->getError() << endl;
    of << "b2 " << ws->var("b2")->getVal() << " " << 0.14*fabs(ws->var("b2")->getVal()) << endl;
    of << "b4 " << ws->var("b4")->getVal() << " " << 0.06*fabs(ws->var("b4")->getVal()) << endl; */
    of << "para2 " << ws->var("para2")->getVal() << " " << 0.01*fabs(ws->var("para2")->getVal()) << endl;
    of << "para4 " << ws->var("para4")->getVal() << " " << 0.01*fabs(ws->var("para4")->getVal()) << endl;
    of << "para6 " << ws->var("para6")->getVal() << " " << 0.01*fabs(ws->var("para6")->getVal()) << endl;
    // of << "para8 " << ws->var("para8")->getVal() << " " << 0.01*fabs(ws->var("para8")->getVal()) << endl;
    of << "acca2 " << acca2.getVal() << " " << acca2.getError() << endl;
    of << "a2 " << ws->var("a2")->getVal() << " " << 0.01*fabs(ws->var("a2")->getVal()) << endl;
    of << "cutOff " << ws->var("cutOff")->getVal() << " 0.01" << endl;
    of << "g " << ws->var("g")->getVal() << " " << ws->var("g")->getError() << endl;
    of << "a4 " << ws->var("a4")->getVal() << " " << 0.01*fabs(ws->var("a4")->getVal()) << endl;
    of << "b2 " << ws->var("b2")->getVal() << " " << 0.14*fabs(ws->var("b2")->getVal()) << endl;
    of << "b4 " << ws->var("b4")->getVal() << " " << 0.06*fabs(ws->var("b4")->getVal()) << endl;
    of.close();
  }

  cout << "***" << endl << "PLOTTING" << endl << "***" << endl;
  c1.cd(1);
  RooPlot *t1frame = ws->var("cosTheta1New")->frame();
  dataNew->plotOn(t1frame,DataError(RooAbsData::SumW2));
  fitf.plotOn(t1frame,Normalization(dataNew->sumEntries(),RooAbsReal::NumEvent),ProjWData(*dataNew),NumCPU(4));
  t1frame->Draw();
  
  c1.cd(2);
  RooPlot *t2frame = ws->var("cosTheta2New")->frame();
  dataNew->plotOn(t2frame,DataError(RooAbsData::SumW2));
  fitf.plotOn(t2frame,Normalization(dataNew->sumEntries(),RooAbsReal::NumEvent),ProjWData(*dataNew),NumCPU(4));
  t2frame->Draw(); 

  c1.cd(3);
  RooPlot *tStarframe = ws->var("cosThetaStarNew")->frame();
  dataNew->plotOn(tStarframe,DataError(RooAbsData::SumW2));
  fitf.plotOn(tStarframe,Normalization(dataNew->sumEntries(),RooAbsReal::NumEvent),ProjWData(*dataNew),NumCPU(4));
  tStarframe->Draw();

  c1.cd(4);
  RooPlot *pframe = ws->var("PhiNew")->frame();
  dataNew->plotOn(pframe,DataError(RooAbsData::SumW2));
  fitf.plotOn(pframe,Normalization(dataNew->sumEntries(),RooAbsReal::NumEvent),ProjWData(*dataNew),NumCPU(4));
  pframe->Draw();

  c1.cd(5);
  RooPlot *p1frame = ws->var("Phi1New")->frame();
  dataNew->plotOn(p1frame,DataError(RooAbsData::SumW2));
  fitf.plotOn(p1frame,Normalization(dataNew->sumEntries(),RooAbsReal::NumEvent),ProjWData(*dataNew));
  p1frame->Draw(); 
  
  if (onlyDraw) sprintf(fileout,"fitCheck%d.ps",int(mass));
  else sprintf(fileout,"fitResult%d.ps",int(mass));
  c1.SaveAs(fileout);

  return 1;
}
