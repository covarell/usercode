Double_t myFunc(Double_t x) {
    if (x >= 0. && x < 28.) return 0.0277*x+0.001605*x*x;
    else if (x >= 28. && x < 54.) return 2.+0.1175*(x-28.)+0.01436*(x-28.)*(x-28.);
    else if (x >= 54. && x < 74.) return 15.+0.864*(x-54.)+0.0207*(x-54.)*(x-54.);
    else if (x >= 74. && x < 91.) return 40.+1.693*(x-74.)+0.00203*(x-74.)*(x-74.); 
    else if (x >= 91. && x < 104.) return 70.+1.763*(x-91.)+0.1709*(x-91.)*(x-91.); 
    return 0;
}
// --> myFunc(36) = 3.86
 
void doCarlosPlots() {
   
   TCanvas *c1 = new TCanvas("c1","Luminosity",200,10,700,500);

   // c1->SetFillColor(42);
   c1->SetGrid();
   // c1->GetFrame()->SetFillColor(21);
   // c1->GetFrame()->SetBorderSize(12);

   const Int_t n = 6;
   Float_t x[n]  = {0., 28., 54., 74., 91., 100.};
   // Float_t y[n]  = {8.0E29, 3.4E30, 2.5E31, 4.9E31, 5.1E31, 1.4E32};
   Float_t y[n]  = {8.0E29, 3.4E30, 2.5E31, 4.9E31, 5.1E31, 6.278E31};

   Float_t ex[n] = {.0001,.0001,.0001,.0001,.0001,.0001};
   Float_t ey[n] = {0.,0.,0.,0.,0.,0.};

   Float_t y2[n]  = {0., 2., 15., 40.,70.,99.7};

   Float_t xd[2]  = {0., 100.};
   Float_t yd[2]  = {-10000.,1500000.};
   Float_t yd2[2]  = {0.001,0.3};

   Float_t exd[n] = {.0001,.0001};
   Float_t eyd[n] = {0.,0.}; 

   TGraphErrors *gr3 = new TGraphErrors(n,x,y,ex,ey);    
   gr3->SetMarkerColor(2);
   gr3->SetLineColor(2);
   gr3->SetLineWidth(3);
   gr3->SetMarkerStyle(20);
   gr3->SetMinimum(-2.E30);
   gr3->GetXaxis()->SetTitle("Time (days)");
   gr3->GetYaxis()->SetTitle("L_{inst} (cm^{-2} s^{-1})");
   gr3->Draw("ALP");

   TLine *hline = new TLine(-10.,1E31,36.,1E31);
   hline->SetLineColor(4);
   hline->SetLineWidth(3);
   hline->Draw("SAME");
   TLine *vline = new TLine(36.,-2E30,36.,1E31);
   vline->SetLineColor(4);
   vline->SetLineWidth(3);
   vline->Draw("SAME");

   c1->Update();
   c1->SaveAs("LumiInst.gif");

   TGraphErrors *gr2 = new TGraphErrors(n,x,y2,ex,ey);
   gr2->SetMarkerColor(2);     
   gr2->SetMinimum(-5.0);
   gr2->SetMaximum(110.0);
   gr2->SetMarkerStyle(20);
   gr2->GetXaxis()->SetTitle("Time (days)");
   gr2->GetYaxis()->SetTitle("L_{integ} (pb^{-1})");
   gr2->Draw("AP");

   TF1 *ftotal = new TF1("ftotal","myFunc(x)",0.,101.); 
   ftotal->SetLineWidth(3);
   ftotal->SetLineColor(2);
   ftotal->Draw("SAME");

   TLine *hline2 = new TLine(-10.,3.85,36.,3.85);
   hline2->SetLineColor(4);
   hline2->SetLineWidth(3);
   hline2->Draw("SAME");
   TLine *vline2 = new TLine(36.,-5.0,36.,3.85);
   vline2->SetLineColor(4);
   vline2->SetLineWidth(3);
   vline2->Draw("SAME");

   c1->Update();
   c1->SaveAs("LumiInteg.gif");
   c1->SetLogy();
   c1->Update();
   c1->SaveAs("LumiIntegLog.gif");

   TCanvas *c2 = new TCanvas("c2","Events",300,100,700,500);  
   c2->cd();
   c2->SetGrid();
 
   TGraphErrors *gr4 = new TGraphErrors(2,xd,yd,exd,eyd);   // DUMMY
   gr4->SetMarkerColor(kWhite);     
   gr4->GetXaxis()->SetTitle("Time (days)");
   gr4->GetYaxis()->SetTitle("N (prompt J/#psi reco)");
  
   TLine *vline3 = new TLine(36.,-20000.0,36.,900000.);  // MENU CHANGE
   vline3->SetLineColor(kMagenta);
   vline3->SetLineStyle(kDashed);
   vline3->SetLineWidth(2);
 
   /// FUNCTIONS: HLTMu3
   TF1 *fhighHLTMu3 = new TF1("fhighHLTMu3","23887*myFunc(x)",0.,36.);
   fhighHLTMu3->SetLineWidth(4);
   fhighHLTMu3->SetLineColor(2);
   TF1 *fmediumHLTMu3 = new TF1("fmediumHLTMu3","30993*myFunc(x)",0.,36.);
   fmediumHLTMu3->SetLineWidth(2);
   fmediumHLTMu3->SetLineColor(4);
   TF1 *flowHLTMu3 = new TF1("flowHLTMu3","66119*myFunc(x)",0.,36.);
   flowHLTMu3->SetLineWidth(2);
   flowHLTMu3->SetLineColor(3);
   TF1 *ftotalHLTMu3 = new TF1("ftotalHLTMu3","121000*myFunc(x)",0.,36.);
   ftotalHLTMu3->SetLineWidth(4);
   ftotalHLTMu3->SetLineColor(1);
   ///
   /// FUNCTIONS: HLTMu5
   TF1 *fhighHLTMu5 = new TF1("fhighHLTMu5","491*myFunc(x)+90309",36.,100.);
   fhighHLTMu5->SetLineWidth(4);
   fhighHLTMu5->SetLineColor(2);
   TF1 *fmediumHLTMu5 = new TF1("fmediumHLTMu5","475*myFunc(x)+117799",36.,100.);
   fmediumHLTMu5->SetLineWidth(2);
   fmediumHLTMu5->SetLineColor(4);
   TF1 *flowHLTMu5 = new TF1("flowHLTMu5","1112*myFunc(x)+250927",36.,100.);
   flowHLTMu5->SetLineWidth(2);
   flowHLTMu5->SetLineColor(3);
   TF1 *ftotalHLTMu5 = new TF1("ftotalHLTMu5","2078*myFunc(x)+459038",36.,100.);
   ftotalHLTMu5->SetLineWidth(4);
   ftotalHLTMu5->SetLineColor(1);
   ///
   /// FUNCTIONS: HLTMu9
   TF1 *fhighHLTMu9 = new TF1("fhighHLTMu9","2430*myFunc(x)+82824",36.,100.);
   fhighHLTMu9->SetLineWidth(4);
   fhighHLTMu9->SetLineColor(2);
   TF1 *fmediumHLTMu9 = new TF1("fmediumHLTMu9","1401*myFunc(x)+114225",36.,100.);
   fmediumHLTMu9->SetLineWidth(2);
   fmediumHLTMu9->SetLineColor(4);
   TF1 *flowHLTMu9 = new TF1("flowHLTMu9","2833*myFunc(x)+244283",36.,100.);
   flowHLTMu9->SetLineWidth(2);
   flowHLTMu9->SetLineColor(3);
   TF1 *ftotalHLTMu9 = new TF1("ftotalHLTMu9","6664*myFunc(x)+441337",36.,100.);
   ftotalHLTMu9->SetLineWidth(4);
   ftotalHLTMu9->SetLineColor(1);
   ///
   /// FUNCTIONS: HLT2Mu3
   TF1 *fhighHLT2Mu3 = new TF1("fhighHLT2Mu3","10271*myFunc(x)+52258",36.,100.);
   fhighHLT2Mu3->SetLineWidth(4);
   fhighHLT2Mu3->SetLineColor(2);
   TF1 *fmediumHLT2Mu3 = new TF1("fmediumHLT2Mu3","1157*myFunc(x)+115167",36.,100.);
   fmediumHLT2Mu3->SetLineWidth(2);
   fmediumHLT2Mu3->SetLineColor(4);
   TF1 *flowHLT2Mu3 = new TF1("flowHLT2Mu3","52*myFunc(x)+255122",36.,100.);
   flowHLT2Mu3->SetLineWidth(2);
   flowHLT2Mu3->SetLineColor(3);
   TF1 *ftotalHLT2Mu3 = new TF1("ftotalHLT2Mu3","11480*myFunc(x)+422747",36.,100.);
   ftotalHLT2Mu3->SetLineWidth(4);
   ftotalHLT2Mu3->SetLineColor(1);
   ///
   
   /// FIRST PLOT: HLTMu3 + HLTMu5
   gr4->Draw("AP");
   fhighHLTMu3->Draw("SAME"); 
   fmediumHLTMu3->Draw("SAME");
   flowHLTMu3->Draw("SAME");
   ftotalHLTMu3->Draw("SAME");
   fhighHLTMu5->Draw("SAME"); 
   fmediumHLTMu5->Draw("SAME");
   flowHLTMu5->Draw("SAME");
   ftotalHLTMu5->Draw("SAME");
   vline3->Draw("SAME"); 

   // LEGENDS
   leg = new TLegend(0.20,0.65,0.60,0.9);
   leg->AddEntry(fhighHLTMu3,"2 global muons","l");
   leg->AddEntry(fmediumHLTMu3,"1 global + 1 tracker muon","l");
   leg->AddEntry(flowHLTMu3,"1 global + 1 calo muon","l");
   leg->AddEntry(ftotalHLTMu3,"all muons","l");
   leg->AddEntry(vline3,"Trigger menu switch","l");
   leg->Draw("SAME");
   hltmu3 = new TPaveLabel(-5.,400000.,20.,550000.,"HLT_Mu3 (1x1)");
   hltmu3->SetTextColor(kMagenta);
   hltmu3->Draw("SAME");
   hltmu5 = new TPaveLabel(65.,1100000.,95.,1250000.,"HLT_Mu5 (25x1)");
   hltmu5->SetTextColor(kMagenta);
   hltmu5->Draw("SAME");

   c2->Update();
   c2->SaveAs("Njpsi_Mu3Mu5.gif");

   /// SECOND PLOT: HLTMu3 + HLTMu9
   gr4->Draw("AP");
   fhighHLTMu3->Draw("SAME"); 
   fmediumHLTMu3->Draw("SAME");
   flowHLTMu3->Draw("SAME");
   ftotalHLTMu3->Draw("SAME");
   fhighHLTMu9->Draw("SAME"); 
   fmediumHLTMu9->Draw("SAME");
   flowHLTMu9->Draw("SAME");
   ftotalHLTMu9->Draw("SAME");
   vline3->Draw("SAME"); 

   // LEGENDS
   leg->Draw("SAME");
   hltmu3->Draw("SAME");
   hltmu9 = new TPaveLabel(65.,1100000.,90.,1250000.,"HLT_Mu9 (1x1)");
   hltmu9->SetTextColor(kMagenta);
   hltmu9->Draw("SAME");

   c2->Update();
   c2->SaveAs("Njpsi_Mu3Mu9.gif");

   /// THIRD PLOT: HLTMu3 + HLT2Mu3
   gr4->Draw("AP");
   fhighHLTMu3->Draw("SAME"); 
   fmediumHLTMu3->Draw("SAME");
   flowHLTMu3->Draw("SAME");
   ftotalHLTMu3->Draw("SAME");
   fhighHLT2Mu3->Draw("SAME"); 
   fmediumHLT2Mu3->Draw("SAME");
   flowHLT2Mu3->Draw("SAME");
   ftotalHLT2Mu3->Draw("SAME");
   vline3->Draw("SAME"); 

   // LEGENDS
   leg->Draw("SAME");
   hltmu3->Draw("SAME");
   hltmu9 = new TPaveLabel(60.,1200000.,87.,1350000.,"HLT_2Mu3 (1x1)");
   hltmu9->SetTextColor(kMagenta);
   hltmu9->Draw("SAME");

   c2->Update();
   c2->SaveAs("Njpsi_Mu32Mu3.gif");

   // S/B
   TCanvas *c3 = new TCanvas("c3","Events",500,300,700,500);  
   c3->cd();
   c3->SetGrid();
   c3->SetLogy();
 
   TGraphErrors *grsb4 = new TGraphErrors(2,xd,yd2,exd,eyd);   // DUMMY
   grsb4->SetMarkerColor(kWhite);     
   grsb4->GetXaxis()->SetTitle("Time (days)");
   grsb4->GetYaxis()->SetTitle("sqrt(S+B)/S (prompt J/#psi reco)");
  
   TLine *vline6 = new TLine(36.,0.0,36.,0.3);  // MENU CHANGE
   vline6->SetLineColor(kMagenta);
   vline6->SetLineStyle(kDashed);
   vline6->SetLineWidth(2);
 
   /// FUNCTIONS: HLTMu3
   TF1 *sbhighHLTMu3 = new TF1("sbhighHLTMu3","sqrt(31649*myFunc(x))/(23887*myFunc(x))",0.,36.);
   sbhighHLTMu3->SetLineWidth(4);
   sbhighHLTMu3->SetLineColor(2);
   TF1 *sbmediumHLTMu3 = new TF1("sbmediumHLTMu3","sqrt(215502*myFunc(x))/(30993*myFunc(x))",0.,36.);
   sbmediumHLTMu3->SetLineWidth(2);
   sbmediumHLTMu3->SetLineColor(4);
   TF1 *sblowHLTMu3 = new TF1("sblowHLTMu3","sqrt(493892*myFunc(x))/(66119*myFunc(x))",0.,36.);
   sblowHLTMu3->SetLineWidth(2);
   sblowHLTMu3->SetLineColor(3);
   TF1 *sbtotalHLTMu3 = new TF1("sbtotalHLTMu3","sqrt(741043*myFunc(x))/(121000*myFunc(x))",0.,36.);
   sbtotalHLTMu3->SetLineWidth(4);
   sbtotalHLTMu3->SetLineColor(1);
   ///
   /// FUNCTIONS: HLTMu5
   TF1 *sbhighHLTMu5 = new TF1("sbhighHLTMu5","sqrt(636*myFunc(x)+119747)/(491*myFunc(x)+90309)",36.,100.);
   sbhighHLTMu5->SetLineWidth(4);
   sbhighHLTMu5->SetLineColor(2);
   TF1 *sbmediumHLTMu5 = new TF1("sbmediumHLTMu5","sqrt(3640*myFunc(x)+821259)/(475*myFunc(x)+117799)",36.,100.);
   sbmediumHLTMu5->SetLineWidth(2);
   sbmediumHLTMu5->SetLineColor(4);
   TF1 *sblowHLTMu5 = new TF1("sblowHLTMu5","sqrt(7050*myFunc(x)+1674140)/(1112*myFunc(x)+250927)",36.,100.);
   sblowHLTMu5->SetLineWidth(2);
   sblowHLTMu5->SetLineColor(3);
   TF1 *sbtotalHLTMu5 = new TF1("sbtotalHLTMu5","sqrt(11326*myFunc(x)+2615146)/(2078*myFunc(x)+459038)",36.,100.);
   sbtotalHLTMu5->SetLineWidth(4);
   sbtotalHLTMu5->SetLineColor(1);
   ///
   /// FUNCTIONS: HLTMu9
   TF1 *sbhighHLTMu9 = new TF1("sbhighHLTMu9","sqrt(3339*myFunc(x)+109315)/(2430*myFunc(x)+82824)",36.,100.);
   sbhighHLTMu9->SetLineWidth(4);
   sbhighHLTMu9->SetLineColor(2);
   TF1 *sbmediumHLTMu9 = new TF1("sbmediumHLTMu9","sqrt(17310*myFunc(x)+768495)/(1401*myFunc(x)+114225)",36.,100.);
   sbmediumHLTMu9->SetLineWidth(2);
   sbmediumHLTMu9->SetLineColor(4);
   TF1 *sblowHLTMu9 = new TF1("sblowHLTMu9","sqrt(12969*myFunc(x)+1856362)/(2833*myFunc(x)+244283)",36.,100.);
   sblowHLTMu9->SetLineWidth(2);
   sblowHLTMu9->SetLineColor(3);
   TF1 *sbtotalHLTMu9 = new TF1("sbtotalHLTMu9","sqrt(33438*myFunc(x)+2734172)/(6664*myFunc(x)+441337)",36.,100.);
   sbtotalHLTMu9->SetLineWidth(4);
   sbtotalHLTMu9->SetLineColor(1);
   ///
   /// FUNCTIONS: HLT2Mu3
   TF1 *sbhighHLT2Mu3 = new TF1("sbhighHLT2Mu3","sqrt(12498*myFunc(x)+73661)/(10271*myFunc(x)+52258)",36.,100.);
   sbhighHLT2Mu3->SetLineWidth(4);
   sbhighHLT2Mu3->SetLineColor(2);
   TF1 *sbmediumHLT2Mu3 = new TF1("sbmediumHLT2Mu3","sqrt(6702*myFunc(x)+809442)/(1157*myFunc(x)+115167)",36.,100.);
   sbmediumHLT2Mu3->SetLineWidth(2);
   sbmediumHLT2Mu3->SetLineColor(4);
   TF1 *sblowHLT2Mu3 = new TF1("sblowHLT2Mu3","sqrt(5233*myFunc(x)+1886327)/(52*myFunc(x)+255122)",36.,100.);
   sblowHLT2Mu3->SetLineWidth(2);
   sblowHLT2Mu3->SetLineColor(3);
   TF1 *sbtotalHLT2Mu3 = new TF1("sbtotalHLT2Mu3","sqrt(24433*myFunc(x)+2769430)/(11480*myFunc(x)+422747)",36.,100.);
   sbtotalHLT2Mu3->SetLineWidth(4);
   sbtotalHLT2Mu3->SetLineColor(1);
   ///
   
   /// FIRST PLOT: HLTMu3 + HLTMu5
   grsb4->Draw("AP");
   sbhighHLTMu3->Draw("SAME"); 
   sbmediumHLTMu3->Draw("SAME");
   sblowHLTMu3->Draw("SAME");
   sbtotalHLTMu3->Draw("SAME");
   sbhighHLTMu5->Draw("SAME"); 
   sbmediumHLTMu5->Draw("SAME");
   sblowHLTMu5->Draw("SAME");
   sbtotalHLTMu5->Draw("SAME");
   vline6->Draw("SAME"); 

   // LEGENDS
   leg = new TLegend(0.50,0.65,0.90,0.9);
   leg->AddEntry(sbhighHLTMu3,"2 global muons","l");
   leg->AddEntry(sbmediumHLTMu3,"1 global + 1 tracker muon","l");
   leg->AddEntry(sblowHLTMu3,"1 global + 1 calo muon","l");
   leg->AddEntry(sbtotalHLTMu3,"all muons","l");
   leg->AddEntry(vline6,"Trigger menu switch","l");
   leg->Draw("SAME");
   hltmu3 = new TPaveLabel(-5.,0.002,20.,0.004,"HLT_Mu3 (1x1)");
   hltmu3->SetTextColor(kMagenta);
   hltmu3->Draw("SAME");
   hltmu5 = new TPaveLabel(65.,0.02,95.,0.03,"HLT_Mu5 (25x1)");
   hltmu5->SetTextColor(kMagenta);
   hltmu5->Draw("SAME");

   c3->Update();
   c3->SaveAs("SBjpsi_Mu3Mu5.gif");

   /// SECOND PLOT: HLTMu3 + HLTMu9
   grsb4->Draw("AP");
   sbhighHLTMu3->Draw("SAME"); 
   sbmediumHLTMu3->Draw("SAME");
   sblowHLTMu3->Draw("SAME");
   sbtotalHLTMu3->Draw("SAME");
   sbhighHLTMu9->Draw("SAME"); 
   sbmediumHLTMu9->Draw("SAME");
   sblowHLTMu9->Draw("SAME");
   sbtotalHLTMu9->Draw("SAME");
   vline6->Draw("SAME"); 

   // LEGENDS
   leg->Draw("SAME");
   hltmu3->Draw("SAME");
   hltmu9 = new TPaveLabel(65.,0.02,90.,0.03,"HLT_Mu9 (1x1)");
   hltmu9->SetTextColor(kMagenta);
   hltmu9->Draw("SAME");

   c3->Update();
   c3->SaveAs("SBjpsi_Mu3Mu9.gif");

   /// THIRD PLOT: HLTMu3 + HLT2Mu3
   grsb4->Draw("AP");
   sbhighHLTMu3->Draw("SAME"); 
   sbmediumHLTMu3->Draw("SAME");
   sblowHLTMu3->Draw("SAME");
   sbtotalHLTMu3->Draw("SAME");
   sbhighHLT2Mu3->Draw("SAME"); 
   sbmediumHLT2Mu3->Draw("SAME");
   sblowHLT2Mu3->Draw("SAME");
   sbtotalHLT2Mu3->Draw("SAME");
   vline6->Draw("SAME"); 

   // LEGENDS
   leg->Draw("SAME");
   hltmu3->Draw("SAME");
   hltmu9 = new TPaveLabel(60.,0.02,87.,0.03,"HLT_2Mu3 (1x1)");
   hltmu9->SetTextColor(kMagenta);
   hltmu9->Draw("SAME");

   c3->Update();
   c3->SaveAs("SBjpsi_Mu32Mu3.gif"); 

   // S/B (pt < 6)
   TCanvas *c4 = new TCanvas("c4","Events",600,400,700,500);  
   c4->cd();
   c4->SetGrid();
   c4->SetLogy();

   /// FUNCTIONS: HLTMu3
   TF1 *sbpt6highHLTMu3 = new TF1("sbpt6highHLTMu3","sqrt(4301*myFunc(x))/(3301*myFunc(x))",0.,36.);
   sbpt6highHLTMu3->SetLineWidth(4);
   sbpt6highHLTMu3->SetLineColor(2);
   TF1 *sbpt6mediumHLTMu3 = new TF1("sbpt6mediumHLTMu3","sqrt(85523*myFunc(x))/(14660*myFunc(x))",0.,36.);
   sbpt6mediumHLTMu3->SetLineWidth(2);
   sbpt6mediumHLTMu3->SetLineColor(4);
   TF1 *sbpt6lowHLTMu3 = new TF1("sbpt6lowHLTMu3","sqrt(174973*myFunc(x))/(35700*myFunc(x))",0.,36.);
   sbpt6lowHLTMu3->SetLineWidth(2);
   sbpt6lowHLTMu3->SetLineColor(3);
   TF1 *sbpt6totalHLTMu3 = new TF1("sbpt6totalHLTMu3","sqrt(264797*myFunc(x))/(53661*myFunc(x))",0.,36.);
   sbpt6totalHLTMu3->SetLineWidth(4);
   sbpt6totalHLTMu3->SetLineColor(1);
   /// FUNCTIONS: HLTMu5 (NON PRESCALED)
   TF1 *sbpt6highHLTMu5 = new TF1("sbpt6highHLTMu5","sqrt(402*myFunc(x)+15050)/(311*myFunc(x)+11541)",36.,100.);
   sbpt6highHLTMu5->SetLineWidth(4);
   sbpt6highHLTMu5->SetLineColor(2);
   TF1 *sbpt6mediumHLTMu5 = new TF1("sbpt6mediumHLTMu5","sqrt(10083*myFunc(x)+291198)/(2902*myFunc(x)+45385)",36.,100.);
   sbpt6mediumHLTMu5->SetLineWidth(2);
   sbpt6mediumHLTMu5->SetLineColor(4);
   TF1 *sbpt6lowHLTMu5 = new TF1("sbpt6lowHLTMu5","sqrt(38740*myFunc(x)+525859)/(11650*myFunc(x)+92833)",36.,100.);
   sbpt6lowHLTMu5->SetLineWidth(2);
   sbpt6lowHLTMu5->SetLineColor(3);
   TF1 *sbpt6totalHLTMu5 = new TF1("sbpt6totalHLTMu5","sqrt(49225*myFunc(x)+832107)/(14863*myFunc(x)+149759)",36.,100.);
   sbpt6totalHLTMu5->SetLineWidth(4);
   sbpt6totalHLTMu5->SetLineColor(1);

   /// FIRST PLOT: HLTMu3 + HLTMu5
   grsb4->GetYaxis()->SetTitle("sqrt(S+B)/S (p_{T} < 6 GeV)");
   grsb4->Draw("AP");
   sbpt6highHLTMu3->Draw("SAME"); 
   sbpt6mediumHLTMu3->Draw("SAME");
   sbpt6lowHLTMu3->Draw("SAME");
   sbpt6totalHLTMu3->Draw("SAME");
   sbpt6highHLTMu5->Draw("SAME"); 
   sbpt6mediumHLTMu5->Draw("SAME");
   sbpt6lowHLTMu5->Draw("SAME");
   sbpt6totalHLTMu5->Draw("SAME");
   vline6->Draw("SAME"); 

   // LEGENDS
   leg->Draw("SAME");
   hltmu3->Draw("SAME");
   hltmu5 = new TPaveLabel(65.,0.02,90.,0.03,"HLT_Mu5 (1x1)");
   hltmu5->SetTextColor(kMagenta);
   hltmu5->Draw("SAME");

   c4->Update();
   c4->SaveAs("SBjpsi_Mu3Mu5_ptlt6.gif"); 

}


