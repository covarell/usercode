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

void plotAvar2(TString whichVar = "pTZ", bool log = false) {
  
  TFile py6("njetsMGOld.root");
  TH1F* theVar6 = (TH1F*)py6.Get(whichVar.Data());
  
  gROOT->ProcessLine(".L ~/myscripts/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
  gStyle->SetOptStat(0);

  theVar6->SetLineWidth(3); 
  theVar6->SetLineColor(kBlack);  
  theVar6->SetMarkerColor(kBlack);
  theVar6->SetMinimum(1.);
  theVar6->Sumw2();
  theVar6->Scale(1./theVar6->Integral());
  
  /* TFile py8n("njetsMGNew.root");
  TH1F* theVar8n = (TH1F*)py8n.Get(whichVar.Data());
  
  theVar8n->SetLineWidth(3); 
  theVar8n->SetLineColor(kRed+1); 
  theVar8n->SetMarkerColor(kRed+1);
  theVar8n->SetMinimum(1.);
  theVar8n->Sumw2();
  theVar8n->Scale(1./theVar8n->Integral());  */

  /* TFile py8p("njetsMGJosh_cteq.root");
  TH1F* theVar8p = (TH1F*)py8p.Get(whichVar.Data());
  
  theVar8p->SetLineWidth(3); 
  theVar8p->SetLineColor(kGreen+2);
  theVar8p->SetMarkerColor(kGreen+2);
  theVar8p->SetMinimum(1.);
  theVar8p->Sumw2();
  theVar8p->Scale(1./theVar8p->Integral());  */

  TFile py8("njetsMGJosh_nnpdf2norwgt.root");
  TH1F* theVar8 = (TH1F*)py8.Get(whichVar.Data());
  
  theVar8->SetLineWidth(3); 
  theVar8->SetLineColor(kMagenta+1);  
  theVar8->SetMarkerColor(kMagenta+1);
  theVar8->SetMinimum(1.);
  theVar8->Sumw2();
  theVar8->Scale(1./theVar8->Integral());

  TFile py56("njetsMGJosh_nnpdf3norwgt.root");
  TH1F* theVar56 = (TH1F*)py56.Get(whichVar.Data());
  
  theVar56->SetLineWidth(3); 
  theVar56->SetLineColor(kBlue+1); 
  theVar56->SetMarkerColor(kBlue+1);
  theVar56->SetMinimum(1.);
  theVar56->Sumw2();
  theVar56->Scale(1./theVar56->Integral());
 

  TFile pynorw("njetsMGJosh_cteqnorwgt.root");
  TH1F* theVarnorw = (TH1F*)pynorw.Get(whichVar.Data());
  
  theVarnorw->SetLineWidth(3); 
  theVarnorw->SetLineColor(kRed-2);  
  theVarnorw->SetMarkerColor(kRed-2);
  theVarnorw->SetMinimum(1.);
  theVarnorw->Sumw2();
  theVarnorw->Scale(1./theVarnorw->Integral());

  TCanvas c1("c1","c1",10,10,600,800);
  c1.cd();
  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.35,0.95,0.97);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.02,0.95,0.35);
  pad2->Draw();

  pad1->cd(); 

  TLegend aa(0.55,0.65,0.95,0.95,"Z+jets (HADRONIZATION LEVEL)");
  aa.AddEntry(theVar6, "MG5+py6 (CTEQ6L)", "l");
  // aa.AddEntry(theVar8n, "MG5_aMC (NNPDF2.3 from GenVal)", "l");
  // aa.AddEntry(theVar8p, "MG5_aMC (CTEQ6L)", "l");
  aa.AddEntry(theVar8, "MG5_aMC+py8 (NNPDF2.3 no reweighting)", "l");
  aa.AddEntry(theVar56, "MG5_aMC+py8 (NNPDF3.0 no reweighting)", "l");
  aa.AddEntry(theVarnorw, "MG5_aMC+py8 (CTEQ6L no reweighting)", "l");

  int imax = 5;
  float theMax = -1.;
  if (theVar6->GetMaximum() > theMax) {
    imax = 0;   theMax = theVar6->GetMaximum();
  }
  /* if (theVar8n->GetMaximum() > theMax) {
    imax = 1;   theMax = theVar8n->GetMaximum();
  }
  if (theVar8p->GetMaximum() > theMax) { 
    imax = 2;   theMax = theVar8p->GetMaximum();
    }*/
  if (theVar8->GetMaximum() > theMax) { 
    imax = 3;   theMax = theVar8->GetMaximum();
  } 
  if (theVar56->GetMaximum() > theMax) { 
    imax = 4;   theMax = theVar56->GetMaximum();
  } 
  if (theVarnorw->GetMaximum() > theMax) imax = 5;

  if (imax == 0) {
    theVar6->Draw("ehist");                                                     
    theVar8->Draw("ehistsame");
    // theVar8n->Draw("ehistsame");
    // theVar8p->Draw("ehistsame");
    theVar56->Draw("ehistsame");
    theVarnorw->Draw("ehistsame");                                                    
    /* } else if (imax == 1) {
    theVar8n->Draw("ehist");
    theVar6->Draw("ehistsame");
    theVar8->Draw("ehistsame");
    theVar8p->Draw("ehistsame");
    theVar56->Draw("ehistsame");
    theVarnorw->Draw("ehistsame");
  } else if (imax == 2) {
    theVar8p->Draw("ehist");
    theVar6->Draw("ehistsame");
    theVar8->Draw("ehistsame");
    theVar8n->Draw("ehistsame");
    theVar56->Draw("ehistsame");
    theVarnorw->Draw("ehistsame");  */
  } else if (imax == 3) {
    theVar8->Draw("ehist");
    theVar6->Draw("ehistsame");
    // theVar8p->Draw("ehistsame");
    // theVar8n->Draw("ehistsame");
    theVar56->Draw("ehistsame");
    theVarnorw->Draw("ehistsame");
  } else if (imax == 4) {
    theVar56->Draw("ehist");
    theVar6->Draw("ehistsame");
    // theVar8p->Draw("ehistsame");
    // theVar8n->Draw("ehistsame");
    theVar8->Draw("ehistsame");
    theVarnorw->Draw("ehistsame");
  } else {
    theVarnorw->Draw("ehist");
    theVar8->Draw("ehistsame");
    // theVar8n->Draw("ehistsame");
    theVar6->Draw("ehistsame");
    // theVar8p->Draw("ehistsame");
    theVar56->Draw("ehistsame");
  }
  if (log) gPad->SetLogy();
  aa.Draw("ehistsame");   

  pad2->cd();
  TH1F* flatline = theDiff(theVar6,theVar6);
  flatline->SetMinimum(-0.2);
  flatline->SetMaximum(0.2);
  flatline->Draw("ehist");
  // TH1F* comp8 = theDiff(theVar6,theVar8n);
  // comp8->Draw("ehistsame");
  TH1F* comp6 = theDiff(theVar6,theVar56);
  comp6->Draw("ehistsame");
  // TH1F* comp8p = theDiff(theVar6,theVar8p);
  // comp8p->Draw("ehistsame");
  TH1F* comp56 = theDiff(theVar6,theVar8);
  comp56->Draw("ehistsame");
  TH1F* compnorw = theDiff(theVar6,theVarnorw);
  compnorw->Draw("ehistsame");

  char namefile[100];
  if (log) 
    sprintf(namefile,"~/www/multiplTest/%s_errors_log.gif",whichVar.Data());
  else sprintf(namefile,"~/www/multiplTest/%s_errors.gif",whichVar.Data());
  c1.SaveAs(namefile);
  
  

}
