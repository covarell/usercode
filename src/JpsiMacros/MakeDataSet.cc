#include <TStyle.h>
#include <TCanvas.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooCategory.h"

#include <iostream>
#include <unistd.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string.h>

#include <math.h>
#include <map>
#include <vector>

#include "MakeDataSet.hh"

using namespace std;
using namespace RooFit;

MakeDataSet::MakeDataSet(TTree *tree) 
  : JPsiTreeBase(tree) {

  onlyTheBest = true;
  efficiencyStore = true;

  /// SELECTION CUTS ///
 
  MIN_nhits_trk = 12;
  MAX_normchi2_trk = 5.0;
  MAX_normchi2_glb = 20.0;
  MIN_nhits_pixel = 2;
  MAX_d0_trk = 5.0;
  MAX_dz_trk = 20.0;
  MIN_vtxprob_jpsi = 0.001;

  // Tag'n'probe stuff
  if (efficiencyStore) {

    TFile *feffTrk = TFile::Open("TnPfiles/fit_result_TkFromSta.root");
    TFile *feffMuGlb = TFile::Open("TnPfiles/fit_result_MuFromTkJPsiGlb.root");
    TFile *feffMuTrk = TFile::Open("TnPfiles/fit_result_MuFromTkJPsiTkM.root");
    TFile *feffMuHLT = TFile::Open("TnPfiles/fit_result_HltFromJPsiGlb.root");
    
    heffTrk = (TH2F*)feffTrk->Get("fit_eff_pt_eta");
    heffTrk->SetName("heffTrk");
    heffMuGlb = (TH2F*)feffMuGlb->Get("fit_eff_pt_eta");
    heffMuGlb->SetName("heffmuGlb");
    heffMuTrk = (TH2F*)feffMuTrk->Get("fit_eff_pt_eta");
    heffMuTrk->SetName("heffmuTrk");
    heffMuHLT = (TH2F*)feffMuHLT->Get("fit_eff_pt_eta");
    heffMuHLT->SetName("heffmuHLT");

  }

  gammaFactor_GGnonprompt         = new TH1F("gammaFactorGG_nonprompt",  "MC correction for non-prompt decays", 100, 0.,3.);
  gammaFactor_GTnonprompt         = new TH1F("gammaFactorGT_nonprompt",  "MC correction for non-prompt decays", 100, 0.,3.);
  hMcPR_GGMass                  = new TH1F("hMcPR_GGMass",  "Inv. mass prompt - MC matched (global+global muons)", 100, 2.6,3.6);
  hMcPR_GTMass                  = new TH1F("hMcPR_GTMass",  "Inv. mass prompt - MC matched (global+tracker muons)", 100, 2.6,3.6);
  hMcNP_GGMass                  = new TH1F("hMcNP_GGMass",  "Inv. mass non-prompt - MC matched (global+global muons)", 100, 2.6,3.6);
  hMcNP_GTMass                  = new TH1F("hMcNP_GTMass",  "Inv. mass non-prompt - MC matched (global+tracker muons)", 100, 2.6,3.6);
  hMcBK_GGMass                  = new TH1F("hMcBK_GGMass",  "Inv. mass background (global+global muons)", 100, 2.6,3.6);
  hMcBK_GTMass                  = new TH1F("hMcBK_GTMass",  "Inv. mass background (global+tracker muons)", 100, 2.6,3.6);
  hMcPR_GGLife                  = new TH1F("hMcPR_GGLife",  "c #tau prompt - MC matched (global+global muons)", 90, -1.0, 3.5);
  hMcPR_GTLife                  = new TH1F("hMcPR_GTLife",  "c #tau prompt - MC matched (global+tracker muons)", 90, -1.0, 3.5);
  hMcNP_GGLife                  = new TH1F("hMcNP_GGLife",  "c #tau non-prompt - MC matched (global+global muons)", 90, -1.0, 3.5);
  hMcNP_GTLife                  = new TH1F("hMcNP_GTLife",  "c #tau non-prompt - MC matched (global+tracker muons)", 90, -1.0, 3.5);
  hMcBK_GGLife                  = new TH1F("hMcBK_GGLife",  "c #tau background (global+global muons)", 90, -1.0, 3.5);
  hMcBK_GTLife                  = new TH1F("hMcBK_GTLife",  "c #tau background (global+tracker muons)", 90, -1.0, 3.5);
  
}

MakeDataSet::~MakeDataSet(){ } 

void MakeDataSet::Loop() {

  // mass-lifetime limits
  const float JpsiMassMin = 2.6;
  const float JpsiMassMax = 3.6;
  const float JpsiCtMin = -1.0;
  const float JpsiCtMax = 3.5;
  float JpsiPtMin = 0.0;           // SET BY 
  float JpsiPtMax = 0.0;           // DEFINITION
  float JpsiEtaMin = 0.0;          // OF BIN
  float JpsiEtaMax = 0.0;          // LIMITS NOW
  
  // B-meson/baryon masses
  const float JpsiNomMass = 3.09688;
  const float B0NomMass = 5.27953;
  const float BpNomMass = 5.27915;
  const float BsNomMass = 5.3663;
  const float LbNomMass = 5.6202;

  ///////////////////////

  if (fChain == 0) return;  
  int nentries = (int)fChain->GetEntries(); 

  // loop over events
  cout << "Number of entries = " << nentries << endl;

  // counters
  int totalEvents = 0;
  int passedCandidates = 0;

  // PT category
  RooCategory JpsiPtType("JpsiPtType","Category of pT");

  //read from file the pt bins
  vector<double> ptmin;
  vector<double> ptmax;

  ifstream fpt;
  fpt.open("ptranges.txt");
  if (!fpt) { cout << "Error opening file ptranges.txt" << endl; assert(0);  }

  Double_t min = 0., max = 0.;  
  //Read in the cache file and store back to array
  while(!fpt.eof()){
    fpt >> min >> max;
    ptmin.push_back(min);
    ptmax.push_back(max);
  }

  JpsiPtMin = ptmin[0];  cout << "Pt min = " << JpsiPtMin << endl;
  JpsiPtMax = ptmax[ptmax.size()-2];  cout << "Pt max = " << JpsiPtMax << endl;

  for(int i=0;i<ptmin.size()-1;i++){
    char catname[100];
    sprintf(catname,"P%d",i+1);
    JpsiPtType.defineType(catname,i+1); 
    cout << "Pt bin " << i+1 << ": Min = " << ptmin[i] << " Max = " << ptmax[i] << endl;   
  }

  // eta category
  RooCategory JpsiEtaType("JpsiEtaType","Category of eta");

  //read from file the pt bins
  vector<double> etamin;
  vector<double> etamax;

  ifstream feta;
  feta.open("etaranges.txt");
  if (!feta){ cout << "Error opening file etaranges.txt" << endl; assert(0); }

  //Read in the cache file and store back to array
  while(!feta.eof()){
    feta >> min >> max;
    etamin.push_back(min);
    etamax.push_back(max);    
  }

  JpsiEtaMin = etamin[0];   cout << "Eta min = " << JpsiEtaMin << endl;
  JpsiEtaMax = etamax[etamax.size()-2];   cout << "Eta max = " << JpsiEtaMax << endl;

  for(int i=0;i<etamin.size()-1;i++){
    char catname[100];
    sprintf(catname,"E%d",i+1);
    JpsiEtaType.defineType(catname,i+1);
    cout << "Eta bin " << i+1 << ": Min = " << etamin[i] << " Max = " << etamax[i] << endl;
  }

  // RooFit stuff
  RooRealVar* JpsiMass = new RooRealVar("JpsiMass","J/psi mass",JpsiMassMin,JpsiMassMax,"GeV/c^{2}");
  RooRealVar* JpsiPt = new RooRealVar("JpsiPt","J/psi pt",JpsiPtMin,JpsiPtMax,"GeV/c");
  RooRealVar* JpsiEta = new RooRealVar("JpsiEta","J/psi eta",JpsiEtaMin,JpsiEtaMax);
  RooRealVar* Jpsict = new RooRealVar("Jpsict","J/psi ctau",JpsiCtMin,JpsiCtMax,"mm");
  RooRealVar* JpsictTrue = new RooRealVar("JpsictTrue","J/psi ctau true",JpsiCtMin,JpsiCtMax,"mm");
  RooRealVar* TNPeff = new RooRealVar("TNPeff","Tag and probe efficiency",0.,1.);
  RooRealVar* TNPefferr = new RooRealVar("TNPefferr","Tag and probe efficiency uncertainty",0.,1.);

  RooRealVar* MCweight = new RooRealVar("MCweight","Monte Carlo Weight",0.,100.);

  //categories for reconstruction
  RooCategory JpsiType("JpsiType","Category of muons");
  JpsiType.defineType("GG",0);
  JpsiType.defineType("GT",1);
  JpsiType.defineType("GC",3);

  //categories for MC production
  RooCategory MCType("MCType","Category of MC");
  MCType.defineType("PR",0);
  MCType.defineType("NP",1);
  MCType.defineType("BK",2);

  int MCcat = -999;

  float weight = 0.;
  float tnpeff = 0.;
  float tnpefferr = 0.;

  RooArgList varlist(*JpsiMass,*Jpsict,*JpsiPt,*JpsiEta,*MCweight,*TNPeff,*TNPefferr,JpsiType,MCType);
  varlist.add(*JpsictTrue);
  varlist.add(JpsiPtType);
  varlist.add(JpsiEtaType);

  RooDataSet* data = new RooDataSet("data","A sample",varlist);

  for (int jentry=0; jentry< nentries; jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);

    if (jentry%100000 == 0) cout << ">>> Processing event # " << jentry << endl;
    
    totalEvents++;

    // TRIGGER CUTS 
    if (!HLTBits_accept[0]) continue;    // SingleMu3

    //estrablish which kind of MC this is
    MCcat = -999;
    TString filestring(fChain->GetCurrentFile()->GetName());

    if (filestring.Contains("promptJpsiMuMu")) {
       MCcat = 0;
       // weight = 0.1089;
       weight = 0.09176;
    } else if (filestring.Contains("BJpsiMuMu")) {
       MCcat = 1;
       // weight = 0.1226;
       weight = 0.0309;
    } else if (filestring.Contains("ppMuX")) {
       MCcat = 2;
       // weight = 12.76;
       weight = 8.180;
    } else if (filestring.Contains("ppMuMu")) {
       MCcat = 2;
       // weight = 1.108;
       weight = 0.892;
    }
     
    //exclude real Jpsi in background MC
    if(MCcat == 2){
      bool aRealJpsiEvent = false;
      for (int imugen=0; imugen<Mc_mu_size; imugen++) {
	if (Mc_mumoth_id[imugen] == 443) aRealJpsiEvent = true;
      }
      if (aRealJpsiEvent) continue;
    }
    // Calculate correction factor/true lifetime for non-prompt
    float corrFact = 1.0;
    float trueLife = -999.;
    bool hasABMother = false;

    for (int iQQgen=0; iQQgen<Mc_QQ_size; iQQgen++) {

      TVector3 *theMcQQvtx = (TVector3*)Mc_QQ_3vec->At(iQQgen);
      TVector3 *theMcBvtx = (TVector3*)Mc_QQmoth_3vec->At(iQQgen);

      TVector3 theDiff = &theMcQQvtx - &theMcBvtx;
      
      TLorentzVector *theMcQQmom = (TLorentzVector*)Mc_QQ_4mom->At(iQQgen);
      TLorentzVector *theMcBmom = (TLorentzVector*)Mc_QQmoth_4mom->At(iQQgen);

      if (abs(Mc_QQmoth_id[iQQgen]) == 511) {
	corrFact = (BpNomMass*theMcQQmom->Perp())/(JpsiNomMass*theMcBmom->Perp());
	hasABMother = true;
      }
      else if (abs(Mc_QQmoth_id[iQQgen]) == 521) {
	corrFact = (B0NomMass*theMcQQmom->Perp())/(JpsiNomMass*theMcBmom->Perp());
      	hasABMother = true;
      }					
      else if (abs(Mc_QQmoth_id[iQQgen]) == 531) {
	corrFact = (BsNomMass*theMcQQmom->Perp())/(JpsiNomMass*theMcBmom->Perp());
	hasABMother = true;
      }
      else if (abs(Mc_QQmoth_id[iQQgen]) == 5122) {
	corrFact = (LbNomMass*theMcQQmom->Perp())/(JpsiNomMass*theMcBmom->Perp());
	hasABMother = true;
      }
      trueLife = theDiff.Perp()*10./corrFact;
    }

    int theMCMatchedGlbMu1 = -1;
    int theMCMatchedGlbMu2 = -1;
    int theMCMatchedTrkMu = -1;
    int theMCMatchedCalMu = -1;
    unsigned int nGlobMuMatched = 0;

    //
    // Find MC matched muons (all QQ categories)
    //
    for (int imugen=0; imugen<Mc_mu_size; imugen++) {
      
      if (Mc_mumoth_id[imugen] == 443) {
        TLorentzVector *theMcMumom = (TLorentzVector*)Mc_mu_4mom->At(imugen);
        for (int imugl=0; imugl<Reco_mu_glb_size; imugl++) {
         
          TLorentzVector *theGlMumom = (TLorentzVector*)Reco_mu_glb_4mom->At(imugl);
          if (deltaR(theMcMumom,theGlMumom) < 0.03) {
            if (nGlobMuMatched == 0) {
	      theMCMatchedGlbMu1 = imugl;
              nGlobMuMatched++;
	    } else {
	      theMCMatchedGlbMu2 = imugl;
	      break;
	    }
	  }
	}
        for (int imutr=0; imutr<Reco_mu_trk_size; imutr++) {
          
          TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(imutr);
          if (deltaR(theMcMumom,theTrMumom) < 0.03) {
	    theMCMatchedTrkMu = imutr;      break;
	  }
	}
         for (int imuca=0; imuca<Reco_mu_cal_size; imuca++) {
	 
          TLorentzVector *theCaMumom = (TLorentzVector*)Reco_mu_cal_4mom->At(imuca);
          if (deltaR(theMcMumom,theCaMumom) < 0.03) {
	    theMCMatchedCalMu = imuca;       break;
	  }
	}
      }      
    }

    // Find the best candidate (if needed)
    int myBest = 0;
    if (onlyTheBest) myBest = theBestQQ();

    for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

      if (Reco_QQ_sign[iqq] != 0) continue;

      // if(Reco_QQ_type[iqq] == 2 || Reco_QQ_type[iqq] > 3) continue;
      if(Reco_QQ_type[iqq] > 1) continue;

      if (onlyTheBest && iqq != myBest) continue;

      TLorentzVector *theQQ4mom = (TLorentzVector*)Reco_QQ_4mom->At(iqq);
      float theMass = theQQ4mom->M();
      float theCtau = Reco_QQ_ctau[iqq]*10.;
      int theLoPtMu = Reco_QQ_mulpt[iqq];  
      int theHiPtMu = Reco_QQ_muhpt[iqq];
      
        // Now, AFTER setting the weight, change to consider MC truth!
      if (filestring.Contains("promptJpsiMuMu") || filestring.Contains("BJpsiMuMu")) {
	bool isMatchedGlbGlb = (theHiPtMu == theMCMatchedGlbMu1 && theLoPtMu == theMCMatchedGlbMu2) || (theHiPtMu == theMCMatchedGlbMu2 && theLoPtMu == theMCMatchedGlbMu1);
	bool isMatchedGlbTrk = (theLoPtMu == theMCMatchedTrkMu && (theHiPtMu == theMCMatchedGlbMu1 || theHiPtMu == theMCMatchedGlbMu2) );
	bool isMatchedGlbCal = (theLoPtMu == theMCMatchedCalMu && (theHiPtMu == theMCMatchedGlbMu1 || theHiPtMu == theMCMatchedGlbMu2) ); 
	if (!isMatchedGlbGlb && !isMatchedGlbTrk && !isMatchedGlbCal) MCcat = 2;
      }
      
      if (theMass > JpsiMassMin && theMass < JpsiMassMax && theCtau > JpsiCtMin && theCtau < JpsiCtMax && theQQ4mom->Perp() > JpsiPtMin && theQQ4mom->Perp() < JpsiPtMax && fabs(theQQ4mom->Eta()) > JpsiEtaMin && fabs(theQQ4mom->Eta()) < JpsiEtaMax) {
	
	passedCandidates++;
	
        // Fill histos
	switch(Reco_QQ_type[iqq]){
	  
	case 0:{
	  
	  switch(MCcat){
	  
	  case 0 :{
	    hMcPR_GGMass->Fill(theMass,weight);
	    hMcPR_GGLife->Fill(theCtau,weight);
	    break;
	  }	    
	  case 1 :{
	    hMcNP_GGMass->Fill(theMass,weight);
	    hMcNP_GGLife->Fill(theCtau,weight);
	    if (hasABMother) gammaFactor_GGnonprompt->Fill(corrFact,weight);
	    break;
	  } 
	  case 2 :{
	    hMcBK_GGMass->Fill(theMass,weight);
	    hMcBK_GGLife->Fill(theCtau,weight);
	    break;
	  }
	  }
	  break;
	}
	
	case 1 :{
	  
	  switch(MCcat){
	    
	  case 0 :{
	    hMcPR_GTMass->Fill(theMass,weight);
	    hMcPR_GTLife->Fill(theCtau,weight);
	    break;
	  }	    
	  case 1 :{
	    hMcNP_GTMass->Fill(theMass,weight);
	    hMcNP_GTLife->Fill(theCtau,weight);
	    if (hasABMother) gammaFactor_GTnonprompt->Fill(corrFact,weight);
	    break;
	  } 
	  case 2 :{
	    hMcBK_GTMass->Fill(theMass,weight);
	    hMcBK_GTLife->Fill(theCtau,weight);
	    break;
	  }
	  }
	  break;
	}
	}
	
        // Fill RooDataSet
	JpsiPt->setVal(theQQ4mom->Perp()); 
	JpsiEta->setVal(theQQ4mom->Eta()); 
	JpsiMass->setVal(theMass);
	Jpsict->setVal(theCtau);
	JpsictTrue->setVal(trueLife);
	JpsiType.setIndex(Reco_QQ_type[iqq],kTRUE);
	MCType.setIndex(MCcat,kTRUE); 
	MCweight->setVal(weight);

	// if efficiencyStore, store efficiencies from TagNProbe
	if (efficiencyStore) {
	  TLorentzVector *theHpt4mom = (TLorentzVector*)Reco_mu_glb_4mom->At(theHiPtMu);
	  TLorentzVector *theLpt4mom;
	  if (Reco_QQ_type[iqq] == 0) 
	    theLpt4mom = (TLorentzVector*)Reco_mu_glb_4mom->At(theLoPtMu);
	  else { 
	    theLpt4mom = (TLorentzVector*)Reco_mu_trk_4mom->At(theLoPtMu);
	  }    

	  float effTrk1 = findEff(heffTrk, theHpt4mom->Pt(), theHpt4mom->Eta(), true);
	  float effTrk2 = findEff(heffTrk, theLpt4mom->Pt(), theLpt4mom->Eta(), true);
	  float effMu1 = findEff(heffMuGlb, theHpt4mom->Pt(), theHpt4mom->Eta(), true);
	  float effMu2;
	  float effHLT1 = findEff(heffMuHLT, theHpt4mom->Pt(), theHpt4mom->Eta(), true);
	  
	  float efferrTrk1 = findEffErr(heffTrk, theHpt4mom->Pt(), theHpt4mom->Eta(), true);
	  float efferrTrk2 = findEffErr(heffTrk, theLpt4mom->Pt(), theLpt4mom->Eta(), true); 
	  float efferrMu1 = findEffErr(heffMuGlb, theHpt4mom->Pt(), theHpt4mom->Eta(), true);
	  float efferrMu2;
	  float efferrHLT1 = findEffErr(heffMuHLT, theHpt4mom->Pt(), theHpt4mom->Eta(), true);
	  
	  
	  if (Reco_QQ_type[iqq] == 0) {
	    effMu2 = findEff(heffMuGlb, theLpt4mom->Pt(), theLpt4mom->Eta(), true);
	    efferrMu2 = findEffErr(heffMuGlb, theLpt4mom->Pt(), theLpt4mom->Eta(), true);
	    float effHLT2 = findEff(heffMuHLT, theLpt4mom->Pt(), theLpt4mom->Eta(), true);
	    float efferrHLT2 = findEffErr(heffMuHLT, theLpt4mom->Pt(), theLpt4mom->Eta(), true);
	    // GLOBAL - GLOBAL : eff_Jpsi ~ eff_(Track from Standalone)^2 * 
	    // eff_(GlobalMu from Track)^2 * 
	    // (2eff_(HLT from GlobalMu) - eff_(HLT from GlobalMu)^2)
	    tnpeff = effTrk1 * effTrk2 * effMu1 * effMu2 * (effHLT1 + effHLT2 - effHLT1*effHLT2);
	    float uglyNumber = (pow((1-effHLT1)*efferrHLT2,2) + pow((1-effHLT2)*efferrHLT1,2))/pow(effHLT1 + effHLT2 - effHLT1*effHLT2,2);
	    tnpefferr = tnpeff * sqrt(pow(efferrTrk1/effTrk1,2) + pow(efferrTrk2/effTrk2,2) + pow(efferrMu1/effMu1,2) + pow(efferrMu2/effMu2,2) + uglyNumber);
	  } else {
	    effMu2 = findEff(heffMuTrk, theLpt4mom->Pt(), theLpt4mom->Eta(), true);
	    efferrMu2 = findEffErr(heffMuTrk, theLpt4mom->Pt(), theLpt4mom->Eta(), true);
	    // GLOBAL - TRACKER : eff_Jpsi ~ eff_(Track from Standalone)^2 * 
	    // eff_(GlobalMu from Track) * eff_(TrackerMu from Track) *
	    // eff_(HLT from GlobalMu)^2
	    tnpeff = effTrk1 * effTrk2 * effMu1 * effMu2 * effHLT1;
	    tnpefferr = tnpeff * sqrt(pow(efferrTrk1/effTrk1,2) + pow(efferrTrk2/effTrk2,2) + pow(efferrMu1/effMu1,2) + pow(efferrMu2/effMu2,2) + pow(efferrHLT1/effHLT1,2));
	  }
	}
	
	TNPeff->setVal(tnpeff);
	TNPefferr->setVal(tnpefferr);
	
	JpsiPtType.setIndex(get_Jpsi_var_type(theQQ4mom->Perp(),ptmin,ptmax),kTRUE);
	JpsiEtaType.setIndex(get_Jpsi_var_type(fabs(theQQ4mom->Eta()),etamin,etamax),kTRUE);
	
	// Fill RooDataSet
	RooArgSet varlist_tmp(*JpsiMass,*Jpsict,*JpsiPt,*JpsiEta,*MCweight,*TNPeff,*TNPefferr,JpsiType,MCType);
        varlist_tmp.add(*JpsictTrue);
	varlist_tmp.add(JpsiPtType);
	varlist_tmp.add(JpsiEtaType);
	  
	data->add(varlist_tmp);
	
      }
    }
  }
  
  data->setWeightVar(*MCweight);

  TCanvas c1("c1","c1",10,10,600,600);
  // TCanvas c1("c1","c1",10,10,1000,500);
  c1.Divide(2,2);
  // c1.Divide(2,1);

  c1.cd(1);
  RooPlot* frameMass = JpsiMass->frame();
  data->plotOn(frameMass,Binning(100),RooFit::Cut("JpsiType==JpsiType::GG"),DataError(RooAbsData::SumW2));
  data->plotOn(frameMass,Binning(100),RooFit::Cut("JpsiType==JpsiType::GG && (MCType==MCType::NP || MCType==MCType::BK)"),LineColor(2),MarkerColor(2),DataError(RooAbsData::SumW2));
  data->plotOn(frameMass,Binning(100),RooFit::Cut("JpsiType==JpsiType::GG && MCType==MCType::BK"),LineColor(4),MarkerColor(4),DataError(RooAbsData::SumW2));
  data->plotOn(frameMass,Binning(100),RooFit::Cut("JpsiType==JpsiType::GG"),DataError(RooAbsData::SumW2));
  frameMass->Draw(); 

  c1.cd(3);
  RooPlot* frameMass2 = JpsiMass->frame();
  data->plotOn(frameMass2,Binning(100),RooFit::Cut("JpsiType==JpsiType::GT"),DataError(RooAbsData::SumW2));
  data->plotOn(frameMass2,Binning(100),RooFit::Cut("JpsiType==JpsiType::GT && (MCType==MCType::NP || MCType==MCType::BK)"),LineColor(2),MarkerColor(2),DataError(RooAbsData::SumW2));
  data->plotOn(frameMass2,Binning(100),RooFit::Cut("JpsiType==JpsiType::GT && MCType==MCType::BK"),LineColor(4),MarkerColor(4),DataError(RooAbsData::SumW2));
  data->plotOn(frameMass2,Binning(100),RooFit::Cut("JpsiType==JpsiType::GT"),DataError(RooAbsData::SumW2));
  frameMass2->Draw(); 

  /* c1.cd(5);
  RooPlot* frameMass3 = JpsiMass->frame(); 
  data->plotOn(frameMass3,Binning(100),RooFit::Cut("JpsiType==JpsiType::GC"));
  frameMass3->Draw(); */

  // c1.cd(1);
  c1.cd(2);
  gPad->SetLogy(1);
  RooPlot* framect = Jpsict->frame();
  data->plotOn(framect,Binning(90),RooFit::Cut("JpsiType==JpsiType::GG"),DataError(RooAbsData::SumW2));
  data->plotOn(framect,Binning(90),RooFit::Cut("JpsiType==JpsiType::GG && (MCType==MCType::NP || MCType==MCType::BK)"),LineColor(2),MarkerColor(2),DataError(RooAbsData::SumW2));
  data->plotOn(framect,Binning(90),RooFit::Cut("JpsiType==JpsiType::GG && MCType==MCType::BK"),LineColor(4),MarkerColor(4),DataError(RooAbsData::SumW2));
  data->plotOn(framect,Binning(90),RooFit::Cut("JpsiType==JpsiType::GG"),DataError(RooAbsData::SumW2));
  framect->Draw();

  // c1.cd(2);
  c1.cd(4);
  gPad->SetLogy(1);
  RooPlot* framect2 = Jpsict->frame();
  data->plotOn(framect2,Binning(90),RooFit::Cut("JpsiType==JpsiType::GT"),DataError(RooAbsData::SumW2));
  data->plotOn(framect2,Binning(90),RooFit::Cut("JpsiType==JpsiType::GT && (MCType==MCType::NP || MCType==MCType::BK)"),LineColor(2),MarkerColor(2),DataError(RooAbsData::SumW2));
  data->plotOn(framect2,Binning(90),RooFit::Cut("JpsiType==JpsiType::GT && MCType==MCType::BK"),LineColor(4),MarkerColor(4),DataError(RooAbsData::SumW2));
  data->plotOn(framect2,Binning(90),RooFit::Cut("JpsiType==JpsiType::GT"),DataError(RooAbsData::SumW2));
  framect2->Draw();

  /* c1.cd(6);
  gPad->SetLogy(1);
  RooPlot* framect3 = Jpsict->frame();
  data->plotOn(framect3,Binning(90),RooFit::Cut("JpsiType==JpsiType::GC"));
  framect3->Draw(); */

  c1.SaveAs("bestCands.gif");
  // c1.SaveAs("lifeTimes.gif");

  TCanvas c2("c2","c2",10,10,600,600);
  c2.Divide(2,2);

  c2.cd(1);
  RooPlot* frameEta = JpsiEta->frame();
  data->plotOn(frameEta,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG"),DataError(RooAbsData::SumW2));
  data->plotOn(frameEta,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG && (MCType==MCType::NP || MCType==MCType::BK)"),LineColor(2),MarkerColor(2),DataError(RooAbsData::SumW2));
  data->plotOn(frameEta,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG && MCType==MCType::BK"),LineColor(4),MarkerColor(4),DataError(RooAbsData::SumW2));
  data->plotOn(frameEta,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG"),DataError(RooAbsData::SumW2));
  frameEta->Draw();

  c2.cd(3);
  RooPlot* frameEta2 = JpsiEta->frame();
  data->plotOn(frameEta2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT"),DataError(RooAbsData::SumW2));
  data->plotOn(frameEta2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT && (MCType==MCType::NP || MCType==MCType::BK)"),LineColor(2),MarkerColor(2),DataError(RooAbsData::SumW2));
  data->plotOn(frameEta2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT && MCType==MCType::BK"),LineColor(4),MarkerColor(4),DataError(RooAbsData::SumW2));
  data->plotOn(frameEta2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT"),DataError(RooAbsData::SumW2));
  frameEta2->Draw();

  /* c2.cd(5);
  RooPlot* frameMass3 = JpsiMass->frame(); 
  data->plotOn(frameMass3,Binning(50),RooFit::Cut("JpsiType==JpsiType::GC"));
  frameMass3->Draw(); */

  c2.cd(2);
  gPad->SetLogy(1);
  RooPlot* framePt = JpsiPt->frame();
  data->plotOn(framePt,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG"),DataError(RooAbsData::SumW2));
  data->plotOn(framePt,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG && (MCType==MCType::NP || MCType==MCType::BK)"),LineColor(2),MarkerColor(2),DataError(RooAbsData::SumW2));
  data->plotOn(framePt,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG && MCType==MCType::BK"),LineColor(4),MarkerColor(4),DataError(RooAbsData::SumW2));
  data->plotOn(framePt,Binning(50),RooFit::Cut("JpsiType==JpsiType::GG"),DataError(RooAbsData::SumW2));
  framePt->Draw();

  c2.cd(4);
  gPad->SetLogy(1);
  RooPlot* framePt2 = JpsiPt->frame();
  data->plotOn(framePt2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT"),DataError(RooAbsData::SumW2));
  data->plotOn(framePt2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT && (MCType==MCType::NP || MCType==MCType::BK)"),LineColor(2),MarkerColor(2),DataError(RooAbsData::SumW2));
  data->plotOn(framePt2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT && MCType==MCType::BK"),LineColor(4),MarkerColor(4),DataError(RooAbsData::SumW2));
  data->plotOn(framePt2,Binning(50),RooFit::Cut("JpsiType==JpsiType::GT"),DataError(RooAbsData::SumW2));
  framePt2->Draw();

  /* c2.cd(6);
  gPad->SetLogy(1);
  RooPlot* framect3 = Jpsict->frame();
  data->plotOn(framect3,Binning(50),RooFit::Cut("JpsiType==JpsiType::GC"));
  framect3->Draw(); */

  c2.SaveAs("bestCands2.gif");
 
  TFile fOut("DataSet.root", "RECREATE");
  fOut.cd();
  data->Write();
  fOut.Close();

  TFile fOut2("FinalHistos.root", "RECREATE");
  fOut2.cd();
  gammaFactor_GGnonprompt  ->Write();
  gammaFactor_GTnonprompt  ->Write();
  hMcPR_GGMass       ->Write();        
  hMcPR_GTMass       ->Write();        
  hMcNP_GGMass       ->Write();       
  hMcNP_GTMass       ->Write();        
  hMcBK_GGMass       ->Write();        
  hMcBK_GTMass       ->Write();        
  hMcPR_GGLife       ->Write();        
  hMcPR_GTLife       ->Write();        
  hMcNP_GGLife       ->Write();        
  hMcNP_GTLife       ->Write();        
  hMcBK_GGLife       ->Write();        
  hMcBK_GTLife       ->Write();       
  fOut2.Close();
   
  // summary
  cout << "total number of events = " << totalEvents << endl;
  cout << "number of passed candidates = " << passedCandidates << endl;
  RooDataSet *reddataPR = (RooDataSet*)data->reduce("JpsiType==JpsiType::GG && MCType==MCType::PR");
  cout << "Weighted global+global (PROMPT) = " << reddataPR->numEntries(kTRUE) << endl;
  RooDataSet *reddataNP = (RooDataSet*)data->reduce("JpsiType==JpsiType::GG && MCType==MCType::NP");
  cout << "Weighted global+global (NON-PROMPT) = " << reddataNP->numEntries(kTRUE) << endl;
  RooDataSet *reddataBK = (RooDataSet*)data->reduce("JpsiType==JpsiType::GG && MCType==MCType::BK");
  cout << "Weighted global+global (BACKGROUND) = " << reddataBK->numEntries(kTRUE) << endl;
  RooDataSet *reddata2PR = (RooDataSet*)data->reduce("JpsiType==JpsiType::GT && MCType==MCType::PR");
  cout << "Weighted global+tracker (PROMPT) = " << reddata2PR->numEntries(kTRUE) << endl;
  RooDataSet *reddata2NP = (RooDataSet*)data->reduce("JpsiType==JpsiType::GT && MCType==MCType::NP");
  cout << "Weighted global+tracker (NON-PROMPT) = " << reddata2NP->numEntries(kTRUE) << endl;
  RooDataSet *reddata2BK = (RooDataSet*)data->reduce("JpsiType==JpsiType::GT && MCType==MCType::BK");
  cout << "Weighted global+tracker (BACKGROUND) = " << reddata2BK->numEntries(kTRUE) << endl;

} // end of program
		

	       
int MakeDataSet::theBestQQ() {
    
  int theBest = -1;
  float thehighestPt = -1.;
 
  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 0 ) {

      int thehptMu = Reco_QQ_muhpt[iqq];   if (thehptMu >= Reco_mu_glb_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];   if (thelptMu >= Reco_mu_glb_size) continue;
      if (Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi && 
	  Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	  Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb && 
	  (((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thehptMu] == 1)) && 
	  fabs(Reco_mu_glb_d0[thehptMu]) < MAX_d0_trk && 
	  fabs(Reco_mu_glb_dz[thehptMu]) < MAX_dz_trk && 
	  Reco_mu_glb_nhitstrack[thelptMu] > MIN_nhits_trk && 
	  Reco_mu_glb_normChi2[thelptMu] < MAX_normchi2_glb &&
	  (((Reco_mu_glb_nhitsPixB[thelptMu] + Reco_mu_glb_nhitsPixE[thelptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thelptMu] + Reco_mu_glb_nhitsPixE[thelptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thelptMu] == 1)) &&
	  fabs(Reco_mu_glb_d0[thelptMu]) < MAX_d0_trk && 
	  fabs(Reco_mu_glb_dz[thelptMu]) < MAX_dz_trk) {
	return iqq;
      }
    }
  }

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 1 ) {
      
      int thehptMu = Reco_QQ_muhpt[iqq];  if (thehptMu >= Reco_mu_glb_size) continue;
      int thelptMu = Reco_QQ_mulpt[iqq];  if (thelptMu >= Reco_mu_trk_size) continue;

      if ( Reco_QQ_probChi2[iqq] > MIN_vtxprob_jpsi &&
	   Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb &&
	   (((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel) || ((Reco_mu_glb_nhitsPixB[thehptMu] + Reco_mu_glb_nhitsPixE[thehptMu]) > MIN_nhits_pixel-1 && Reco_mu_glb_nhitsPix1Hit[thehptMu] == 1)) &&
	   fabs(Reco_mu_glb_d0[thehptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_glb_dz[thehptMu]) < MAX_dz_trk && 
	   Reco_mu_trk_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   ((Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,5))/(int)pow(2,5) > 0 || (Reco_mu_trk_PIDmask[thelptMu] & (int)pow(2,8))/(int)pow(2,8) > 0) &&
	   (Reco_mu_trk_nhitsPixB[thelptMu] + Reco_mu_trk_nhitsPixE[thelptMu]) > MIN_nhits_pixel &&
	   Reco_mu_trk_normChi2[thelptMu] < MAX_normchi2_trk &&
	   fabs(Reco_mu_trk_d0[thelptMu]) < MAX_d0_trk && 
	   fabs(Reco_mu_trk_dz[thelptMu]) < MAX_dz_trk) {
	
        TLorentzVector *theTrMumom = (TLorentzVector*)Reco_mu_trk_4mom->At(thelptMu);
        if (theTrMumom->Perp() > thehighestPt) {
	  thehighestPt = theTrMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  if (theBest >= 0) return theBest;

  for (int iqq=0; iqq<Reco_QQ_size; iqq++) {

    if (Reco_QQ_sign[iqq] == 0 && Reco_QQ_type[iqq] == 3 ) {

      int thehptMu = Reco_QQ_muhpt[iqq];
      int thelptMu = Reco_QQ_mulpt[iqq];
      if (thelptMu >= Reco_mu_cal_size) {
	// cout << "Non deve succedere! cmIndex = " << thelptMu+1 << " cmSize = " << Reco_mu_cal_size << endl;
	continue;
      }

      if ( Reco_mu_glb_nhitstrack[thehptMu] > MIN_nhits_trk && 
	   Reco_mu_glb_normChi2[thehptMu] < MAX_normchi2_glb &&
	   Reco_mu_cal_nhitstrack[thelptMu] > MIN_nhits_trk && 
	   Reco_mu_cal_normChi2[thelptMu] < 3.0 && 
	   Reco_mu_cal_caloComp[thelptMu] > 0.89) {
	
        TLorentzVector *theCaMumom = (TLorentzVector*)Reco_mu_cal_4mom->At(thelptMu);
        if (theCaMumom->Perp() > thehighestPt) {
	  thehighestPt = theCaMumom->Perp();
          theBest = iqq;
	}
      }
    }    
  }
  
  return theBest;
} 			       

int MakeDataSet::get_Jpsi_var_type(const double jpsivar, vector<double> vectmin, vector<double> vectmax) {

  for(int i=0;i<vectmin.size()-1;i++) {
    if(jpsivar > vectmin[i] && jpsivar < vectmax[i]) return i+1;
  }

  return -999;
}


double MakeDataSet::PhiInRange(const double& phi) const {
      double phiout = phi;

      if( phiout > 2*M_PI || phiout < -2*M_PI) {
            phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;

      return phiout;
}

double MakeDataSet::deltaR(const TLorentzVector* t, const TLorentzVector* u) const {
      return sqrt(pow(t->Eta()-u->Eta(),2) +pow(PhiInRange(t->Phi()-u->Phi()),2));
}

float MakeDataSet::findEff(TH2F* effhist, float pt, float eta, bool approx) const {

  if (effhist->GetBinContent( effhist->FindBin(pt,eta) ) >  0.0001) {
    return effhist->GetBinContent( effhist->FindBin(pt,eta) );
  } else {
    if (approx) {
 
      // try close-by bins
      int binp;       int bineta;        int dummy;
      effhist->GetBinXYZ(effhist->FindBin(pt,eta),binp,bineta,dummy);
      int newbinp = binp - 1;
      if (effhist->GetBinContent( effhist->GetBin(newbinp,bineta,dummy) ) > 0.0001)
	return effhist->GetBinContent( effhist->GetBin(newbinp,bineta,dummy) );
      newbinp = binp + 1;
      if (effhist->GetBinContent( effhist->GetBin(newbinp,bineta,dummy) ) > 0.0001)
	return effhist->GetBinContent( effhist->GetBin(newbinp,bineta,dummy) );
    } 
  }
  return 1.;
}

float MakeDataSet::findEffErr(TH2F* effhist, float pt, float eta, bool approx) const {

  if (effhist->GetBinContent( effhist->FindBin(pt,eta) ) >  0.0001) {
    return effhist->GetBinError( effhist->FindBin(pt,eta) );
  } else {
    if (approx) {

      // try close-by bins
      int binp;       int bineta;        int dummy;
      effhist->GetBinXYZ(effhist->FindBin(pt,eta),binp,bineta,dummy);
      int newbinp = binp - 1;
      if (effhist->GetBinContent( effhist->GetBin(newbinp,bineta,dummy) ) > 0.0001)
	return effhist->GetBinError( effhist->GetBin(newbinp,bineta,dummy) );
      newbinp = binp + 1;
      if (effhist->GetBinContent( effhist->GetBin(newbinp,bineta,dummy) ) > 0.0001)
	return effhist->GetBinError( effhist->GetBin(newbinp,bineta,dummy) );
    } 
  }
  return 0.;
}
