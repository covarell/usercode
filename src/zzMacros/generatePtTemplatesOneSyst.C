/* 
 * Create 2D (mass, LD) templates. Script imported from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JHU/MELA/scripts/generateTemplates.C?revision=1.1.2.1&view=markup&pathrev=post_unblinding
  * usage: 
 * -set input paths variables in Config.h
 * -run with:
 * > root -q -b 
 * > .L RooModifTsallis.cc+ 
 * > .x generatePtTemplatesOneSyst.C+
 * 2D templates are written to "destDir"
 *
 */

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include <sstream>
#include <vector>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
//#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif

#include "RooRealVar.h"
#include "RooModifTsallis.h"
//----------> SET INPUT VARIABLES in Config.h
#include "Config.h"
//<----------

//---
int useSqrts=0;              //0=use 7+8TeV; 1=use 7TeV only, 2 use 8TeV only
TString melaName = "ZZLD"; // name of MELA branch to be used. Possibilities are ZZLD,ZZLD_analBkg,ZZLD_postICHEP,ZZLD_PtY,pseudoMelaLD, spinTwoMinimalMelaLD 

bool extendToHighMass = true; // Include signal samples above 600 GeV

float highMzz=(extendToHighMass?1000:800);
float mBinSize=2.;

const TString destDir = "../CreateDatacards/templates2D/";
static const int nsamp = 8;
const TString dataFileNames[nsamp] = {"gg","vbf","wh","zh","tth","zz","zx","ggzz"};
TString systSources[nsamp][5];

//=======================================================================

TH2F* fillTemplate(TString channel="4mu", int sampleIndex=0,int LHCsqrts= 7,bool isLowMass=true, TString systName = "Default", bool down = false){  
 
  static const int ptmbins = 20;
  double ptmbinLimits[ptmbins+1] = {0.,.1,.2,.3,.4,.5,.6,.7,
                                    .8,.9,1.0,1.1,1.2,1.3,1.4,1.5,
                                    1.75,2.0,2.25,2.5,4.};

  TH2F* bkgHist;
  if(!isLowMass)
    // Fine binning (constant)
    // bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-100.)/mBinSize+0.5),100,highMzz,500,0,4);
    // Coarse binning (constant)
    bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-100.)/mBinSize+0.5),100,highMzz,50,0,3);
    // Coarse binning (variable)
    // bkgHist = new TH2F("bkgHisto","bkgHisto",int((highMzz-100.)/mBinSize+0.5),100,highMzz,ptmbins,ptmbinLimits);
  else
    // Fine binning (constant)
    // bkgHist = new TH2F("bkgHisto","bkgHisto",int(35/mBinSize+0.5),106,141,500,0,4);
    // Coarse binning (constant)
    bkgHist = new TH2F("bkgHisto","bkgHisto",int(35/mBinSize+0.5),106,141,50,0,3);
    // Coarse binning (variable)
    // bkgHist = new TH2F("bkgHisto","bkgHisto",int(35/mBinSize+0.5),106,141,ptmbins,ptmbinLimits);

  bkgHist->Sumw2();

  if (sampleIndex < 0 || sampleIndex > nsamp) {
    cout << "Samples are numbered from 0 to " << nsamp-1 << endl;
    return bkgHist;
  }

  // fill histogram
  RooRealVar* ptoverm = new RooRealVar("ptoverm","p_{T}/M^{4l}",0.,4.,"GeV/c");

  RooRealVar mup("mup","emme", 1.,-1000000.,1000000.);
  RooRealVar nup("nup","enne", 0.93, -1000000.,1000000.);
  RooRealVar n2up("n2up","enne2", 0.75, -1000000.,1000000.);
  RooRealVar bbup("bbup","bibi",0.02, -1000000.,1000000.);
  RooRealVar Tup("Tup","tti",0.2,-1000000.,1000000.);
  RooRealVar bb2up("bb2up","bibi2",0.02, -1000000.,1000000.);
  RooRealVar fexpup("fexpup","f_exp",0.02, -1000000.,1000000.);
 
  RooModifTsallis* rtup = new RooModifTsallis("rtup","rtup",*ptoverm,
					    mup,nup,n2up,bbup,bb2up,Tup,fexpup);
  RooRealVar m("m","emme", 1.,-1000000.,1000000.);
  RooRealVar n("n","enne", 0.93, -1000000.,1000000.);
  RooRealVar n2("n2","enne2", 0.75, -1000000.,1000000.);
  RooRealVar bb("bb","bibi",0.02, -1000000.,1000000.);
  RooRealVar T("T","tti",0.2,-1000000.,1000000.);
  RooRealVar bb2("bb2","bibi2",0.02, -1000000.,1000000.);
  RooRealVar fexp("fexp","f_exp",0.02, -1000000.,1000000.);
 
  RooModifTsallis* rt = new RooModifTsallis("rt","rt",*ptoverm,
					    m,n,n2,bb,bb2,T,fexp);

  // default
  char fileName[200];
  int nXbins=bkgHist->GetNbinsX();
  int nYbins=bkgHist->GetNbinsY();
    
  sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s125_%dTeV_Default.txt",dataFileNames[sampleIndex].Data(),LHCsqrts);
  if (sampleIndex == 4) sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s125_%dTeV_Default.txt",dataFileNames[6].Data(),LHCsqrts);
  ifstream theFile(fileName);
  (RooArgSet(mup,nup,n2up,bbup,bb2up,fexpup,Tup)).readFromStream(theFile,false);   
  if (systName == "Default") {
    
    for (Int_t i=1; i<=nYbins; i++) {      
      ptoverm->setVal(bkgHist->GetYaxis()->GetBinCenter(i));
      // if (sampleIndex == 0 && LHCsqrts == 8) cout << i << " " << ptoverm->getVal() << " " << rt->getVal() << endl; 
      for (Int_t j=1; j<=nXbins; j++) {
	// float binVolume = bkgHist->GetYaxis()->GetBinWidth(i)*bkgHist->GetXaxis()->GetBinWidth(j);
	bkgHist->SetBinContent(j,i,rtup->getVal());
	// if (i == 20 && j == 20) cout << "Default " << bkgHist->GetBinContent(j,i) << endl;
      }
    } 

  } else {

    float mtotVar = mup.getError()*mup.getError();
    float ntotVar = nup.getError()*nup.getError();
    float n2totVar = n2up.getError()*n2up.getError();
    float bbtotVar = bbup.getError()*bbup.getError();
    float bb2totVar = bb2up.getError()*bb2up.getError();
    float fexptotVar = fexpup.getError()*fexpup.getError();
    float TtotVar = T.getError()*T.getError();

    for (int ss = 0; ss < 5; ss++) {
      if (systSources[sampleIndex][ss] != "") {
    
	if (systSources[sampleIndex][ss] == "Mela") {
    
	  if (down) {
	    sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s125_%dTeV_Mela00-03.txt",dataFileNames[sampleIndex].Data(),LHCsqrts);
	  } else {
	    sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s125_%dTeV_Mela06-10.txt",dataFileNames[sampleIndex].Data(),LHCsqrts);
	  }

	} else {
    
	  sprintf(fileName,"tsallisPars/paramsPTOverMCJLST_%s125_%dTeV_%s.txt",dataFileNames[sampleIndex].Data(),LHCsqrts,systSources[sampleIndex][ss].Data());
	      
	}

	ifstream theFile2(fileName);
	(RooArgSet(m,n,n2,bb,bb2,fexp,T)).readFromStream(theFile2,false); 

	mtotVar += (mup.getVal() - m.getVal())*(mup.getVal() - m.getVal());
	ntotVar += (nup.getVal() - n.getVal())*(nup.getVal() - n.getVal());
	n2totVar += (n2up.getVal() - n2.getVal())*(n2up.getVal() - n2.getVal());
	bbtotVar += (bbup.getVal() - bb.getVal())*(bbup.getVal() - bb.getVal());
	bb2totVar += (bb2up.getVal() - bb2.getVal())*(bb2up.getVal() - bb2.getVal());
	fexptotVar += (fexpup.getVal() - fexp.getVal())*(fexpup.getVal() - fexp.getVal());
	TtotVar += (Tup.getVal() - T.getVal())*(Tup.getVal() - T.getVal());
      	
      }
    }

    mtotVar = sqrt(mtotVar);
    ntotVar = sqrt(ntotVar);
    n2totVar = sqrt(n2totVar);
    bbtotVar = sqrt(bbtotVar);
    bb2totVar = sqrt(bb2totVar);
    fexptotVar = sqrt(fexptotVar);
    TtotVar = sqrt(TtotVar);

    float sign = 1.;
    if (down) sign = -1.;
    m.setVal(mup.getVal() + sign*mtotVar);
    n.setVal(nup.getVal() + sign*ntotVar);
    n2.setVal(n2up.getVal() + sign*n2totVar);
    bb.setVal(bbup.getVal() + sign*bbtotVar);
    bb2.setVal(bb2up.getVal() + sign*bb2totVar);
    fexp.setVal(fexpup.getVal() + sign*fexptotVar);
    T.setVal(Tup.getVal() + sign*TtotVar);

    for (Int_t i=1; i<=nYbins; i++) {      
      ptoverm->setVal(bkgHist->GetYaxis()->GetBinCenter(i));
      // if (sampleIndex == 0 && LHCsqrts == 8) cout << i << " " << ptoverm->getVal() << " " << rt->getVal() << endl; 
      for (Int_t j=1; j<=nXbins; j++) {
	// float binVolume = bkgHist->GetYaxis()->GetBinWidth(i)*bkgHist->GetXaxis()->GetBinWidth(j);
	bkgHist->SetBinContent(j,i,rt->getVal());
	// if (i == 20 && j == 20) cout << "OneSyst " << bkgHist->GetBinContent(j,i) << endl;
      }
    } 
  }

  // normalize slices

  double norm;
  TH1F* tempProj;
  
  for(int i=1; i<=nXbins; i++){
    
    tempProj = (TH1F*)bkgHist->ProjectionY("tempProj",i,i);
    norm = tempProj->Integral();
    // if (sampleIndex == 0 && LHCsqrts == 8) cout << i << " " << norm << endl;

    if (norm>0) { // Avoid introducing NaNs in the histogram
      for(int j=1; j<=nYbins; j++){
	bkgHist->SetBinContent(i,j,bkgHist->GetBinContent(i,j)/norm);
      }
    }

  }
  
  return bkgHist;
  
}

//=======================================================================

void makeTemplate(TString channel="4mu",bool isLowMass=true){

  char fileName[200];

  for (Int_t lhc=7; lhc<9; lhc++) {
    for (Int_t k=0; k<nsamp; k++) {
     
      TString lhcs = "7";
      if (lhc==8) lhcs="8";

      TFile* ftemp = new TFile(destDir + "PtoverM_mZZ_" + dataFileNames[k] + "_" + channel + "_"+ lhcs + "TeV.root","RECREATE");

      cout << "*** Now filling: " << dataFileNames[k] << ", Default template";
      TH2F* h_mzzD = fillTemplate(channel,k,lhc,isLowMass);
     
      ftemp->cd();
  
      h_mzzD->Write("h_Ptmzz_mzz");
      cout << " " << h_mzzD->GetBinContent(100,10) << endl;
      
      cout << "*** Now filling: " << dataFileNames[k] << ", OneSyst templates";
      
      TH2F* h_mzzDup = fillTemplate(channel,k,lhc,isLowMass,"OneSyst",false);
      sprintf(fileName,"h_Ptmzz_mzz_OneSyst_up");
      h_mzzDup->Write(fileName);
      cout << " " << h_mzzDup->GetBinContent(100,10);
      
      TH2F* h_mzzDdown = fillTemplate(channel,k,lhc,isLowMass,"OneSyst",true);
      sprintf(fileName,"h_Ptmzz_mzz_OneSyst_down");
      h_mzzDdown->Write(fileName);
      cout << " " << h_mzzDdown->GetBinContent(100,10) << endl;

      ftemp->Close();
    }
  }

}

//=======================================================================

void generatePtTemplatesOneSyst(){

  bool isLowMass = false;

  systSources[0][0] = "Resummation";
  systSources[0][1] = "TopMass";
  systSources[0][2] = "Mela";
  
  systSources[1][0] = "PDF-VBF";
  systSources[1][1] = "scale-VBF";
  systSources[1][2] = "Mela";
  
  systSources[2][0] = "NLOLO_WH";
  
  systSources[3][0] = "NLOLO_ZH";

  systSources[5][0] = "SingleZ";
  systSources[5][1] = "PDF-ZZ";
  systSources[5][2] = "scale-ZZ";
  systSources[5][3] = "Mela";

  makeTemplate("4mu",isLowMass);
  makeTemplate("4e",isLowMass);
  makeTemplate("2e2mu",isLowMass);

}



