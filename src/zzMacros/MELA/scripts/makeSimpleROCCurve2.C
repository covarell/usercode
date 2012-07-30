// C++ includes
#include <iostream>
#include <string>
#include <stdio.h>
#include <sstream>

// ROOT includes
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>

void makeSimpleROCCurve2(string sigfile, string bkgfile) {

  gROOT->ProcessLine(".L mytdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");

  string sigfile1 = sigfile + "_withDiscriminants.root";
  string bkgfile1 = bkgfile + "_withDiscriminants.root";

  TFile fsig(sigfile1.c_str());
  TFile fbkg(bkgfile1.c_str());

  // f.cd(); 
  TTree* tsig = (TTree*)fsig.Get("angles");
  tsig->SetName("tsig");
  TTree* tbkg = (TTree*)fbkg.Get("angles");
  tbkg->SetName("tbkg");

  float effsig[101];
  float effbkg[101];
  char cutstr[100];
  for (int i = 0; i <= 100; i++) {
     float cutval = 0.01*i;
     sprintf(cutstr,"MC_weight*(melaLD > %f)",cutval);
     tsig->Draw("melaLD >> ldh",cutstr);
     sprintf(cutstr,"MC_weight*(melaLD > 0.0)");
     tsig->Draw("melaLD >> ldh2",cutstr);
     effsig[i] = (float)ldh->Integral()/(float)ldh2->Integral();
     sprintf(cutstr,"MC_weight*(melaLD > %f)",cutval);
     tbkg->Draw("melaLD >> ldh",cutstr);
     sprintf(cutstr,"MC_weight*(melaLD > 0.0)");
     tbkg->Draw("melaLD >> ldh2",cutstr);
     effbkg[i] = (float)ldh->Integral()/(float)ldh2->Integral();
  }	
 
  TGraph groc(101,effsig,effbkg);
  groc.SetMarkerStyle(20);
  groc.SetMarkerColor(kRed);
  groc.GetXaxis()->SetTitle("#epsilon_{sig}");
  groc.GetYaxis()->SetTitle("#epsilon_{bkg}");

  /* string sigfile2 = sigfile + "_withDiscriminantsAndPt.root";
  string bkgfile2 = bkgfile + "_withDiscriminantsAndPt.root";

  TFile fsig2(sigfile2.c_str());
  TFile fbkg2(bkgfile2.c_str());  

  // f.cd(); 
  TTree* tsig2 = (TTree*)fsig2.Get("angles");
  tsig->SetName("tsig");
  TTree* tbkg2 = (TTree*)fbkg2.Get("angles");
  tbkg->SetName("tbkg");  */

  float effsig2[101];
  float effbkg2[101];
  for (int i = 0; i <= 100; i++) {
    float cutval = 0.01*i;
    sprintf(cutstr,"MC_weight*(melaLDWithPt > %f)",cutval);
    tsig->Draw("melaLDWithPt >> ldh",cutstr);
    sprintf(cutstr,"MC_weight*(melaLDWithPt > 0.0)");
    tsig->Draw("melaLDWithPt >> ldh2",cutstr);
    effsig2[i] = (float)ldh->Integral()/(float)ldh2->Integral();
    sprintf(cutstr,"MC_weight*(melaLDWithPt > %f)",cutval);
    tbkg->Draw("melaLDWithPt >> ldh",cutstr);
    sprintf(cutstr,"MC_weight*(melaLDWithPt > 0.0)");
    tbkg->Draw("melaLDWithPt >> ldh2",cutstr);
    effbkg2[i] = (float)ldh->Integral()/(float)ldh2->Integral();
  }	
 
  TGraph groc2(101,effsig2,effbkg2);
  groc2.SetMarkerStyle(21);
  groc2.SetMarkerColor(kBlue);
  groc2.GetXaxis()->SetTitle("#epsilon_{sig}");
  groc2.GetYaxis()->SetTitle("#epsilon_{bkg}");

  TFile fchris("../datafiles/Roc_curves_120_130.root");
  TGraph* chris1 = (TGraph*)fchris.Get("MELAroc");
  chris1->SetName("chris1");
  chris1->SetMarkerStyle(5);
  chris1->SetMarkerColor(kBlack);
  TGraph* chris2 = (TGraph*)fchris.Get("MELApTroc");
  // TGraph* chris2 = (TGraph*)fchris.Get("pTroc");
  chris2->SetName("chris2");
  chris2->SetMarkerStyle(5);
  chris2->SetMarkerColor(kBlack);

  TCanvas d;
  gPad->SetGridx();
  gPad->SetGridy();
  groc.Draw("APL");
  chris1->Draw("PLSAME");
  TLegend * leg = new TLegend(0.18,0.56,0.49,0.835,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->AddEntry(&groc,"RC: w/o pT","p");
  leg->AddEntry(chris1,"CM: w/o pT","p");
  leg->Draw("same");

  d.SaveAs("comparMela.gif");
    
  groc2.Draw("APL");
  chris2->Draw("PLSAME");
  TLegend * leg2 = new TLegend(0.18,0.56,0.49,0.835,NULL,"brNDC");
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetShadowColor(0);
  // leg2->AddEntry(&groc2,"RC: w/ pT","p");
  // leg2->AddEntry(chris2,"CM: w/ pT","p");
  leg2->AddEntry(&groc2,"RC: pT only","p");
  leg2->AddEntry(chris2,"CM: pT only","p");
  leg2->Draw("same");
  
  d.SaveAs("comparMelaPt.gif");

  sprintf(cutstr,"MC_weight*(melaLDWithPt > 0.0)");
  tsig->Draw("melaLDWithPt",cutstr);

  TH1F* chris3 = (TH1F*)fchris.Get("p_{T sig}");
  chris3->SetName("chris3");
  chris3->SetLineStyle(kDashed);
  // chris3->Draw("SAME");

  d.SaveAs("comparPsig.gif");

  sprintf(cutstr,"MC_weight*(melaLDWithPt > 0.0)");
  tbkg->Draw("melaLDWithPt",cutstr);  

  TH1F* chris4 = (TH1F*)fchris.Get("p_{T bkg}");
  chris4->SetName("chris4");
  chris4->SetLineStyle(kDashed);
  // chris4->Draw("SAME");

  d.SaveAs("comparPbkg.gif");

  return;
}

