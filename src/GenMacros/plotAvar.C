TH1F* theDiff(TH1F* ref, TH1F* con) {
  TH1F* normref = (TH1F*)ref->Clone();
  TH1F* normcon = (TH1F*)con->Clone();
  normref->Scale(1./normref->Integral());
  normcon->Scale(1./normcon->Integral());
  TH1F* result = (TH1F*)normcon->Clone();
  TH1F* temp   = (TH1F*)normcon->Clone();
  temp->Add(normcon,normref,1,-1);
  result->Divide(temp, normref);
  return result;
}

void plotAvar(TString whichVar = "pTZ", bool log = false) {
  
  TFile py6("pythia6.root");
  TH1F* theVar6 = (TH1F*)py6.Get(whichVar.Data());
  
  gROOT->ProcessLine(".L ~/myscripts/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  gStyle->SetOptStat(0);

  theVar6->SetLineWidth(3); 
  theVar6->SetLineColor(kBlack);
  theVar6->SetMinimum(1.);
  // theVar6->Sumw2();
 
  TFile py8n("pythia8_noEV.root");
  TH1F* theVar8n = (TH1F*)py8n.Get(whichVar.Data());
  
  theVar8n->SetLineWidth(3); 
  theVar8n->SetLineColor(kRed+1);
  theVar8n->SetMinimum(1.);
  // theVar8n->Sumw2();

  TFile py8("pythia8_EV.root");
  TH1F* theVar8 = (TH1F*)py8.Get(whichVar.Data());
  
  theVar8->SetLineWidth(3); 
  theVar8->SetLineColor(kGreen+2);
  theVar8->SetMinimum(1.);
  // theVar8->Sumw2();

  TCanvas c1("c1","c1",10,10,600,800);
  c1.cd();
  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
  pad2->Draw();

  pad1->cd(); 

  TLegend aa(0.55,0.65,0.95,0.95,"Z jets POWHEG +");
  aa.AddEntry(theVar6, "PYTHIA 6", "l");
  aa.AddEntry(theVar8n, "PYTHIA 8 no e.v.", "l");
  aa.AddEntry(theVar8, "PYTHIA 8 e.v.", "l");

  int imax = 2;
  float theMax = -1.;
  if (theVar6->GetMaximum() > theMax) {
    imax = 0;   theMax = theVar6->GetMaximum();
  }
  if (theVar8n->GetMaximum() > theMax) {
    imax = 1;   theMax = theVar8n->GetMaximum();
  }
  if (theVar8->GetMaximum() > theMax) imax = 2;

  if (imax = 0) {
    theVar6->Draw();
    theVar8->Draw("same");
    theVar8n->Draw("same");
  } else if (imax = 1);
    theVar8n->Draw();
    theVar6->Draw("same");
    theVar8->Draw("same");
  else {
    theVar8->Draw();
    theVar8n->Draw("same");
    theVar6->Draw("same");
  }
  if (log) gPad->SetLogy();
  aa.Draw("same");   

  pad2->cd();
  TH1F* flatline = theDiff(theVar8,theVar8);
  flatline->SetMinimum(-0.1);
  flatline->SetMaximum(0.1);
  flatline->Draw();
  TH1F* comp8 = theDiff(theVar8,theVar8n);
  comp8->Draw("same");
  TH1F* comp6 = theDiff(theVar8,theVar6);
  comp6->Draw("same");

  char namefile[100];
  sprintf(namefile,"~/www/powheg/emissionVetoTest/%s.gif",whichVar.Data());
  c1.SaveAs(namefile);
  
  

}
