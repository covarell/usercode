void doFabrizioPlots() {
   
   TCanvas *c1 = new TCanvas("c1","c1",200,10,800,400);
   c1->Divide(2,1);

   TFile f("FinalHistos.root");
   TH1F* hMcPR_GGMass = (TH1F*)f.Get("hMcPR_GGMass");
   TH1F* hMcNP_GGMass = (TH1F*)f.Get("hMcNP_GGMass");
   TH1F* hMcBK_GGMass = (TH1F*)f.Get("hMcBK_GGMass");

   TH1F* hMcBKNP_GGMass = (TH1F*)hMcBK_GGMass->Clone();
   hMcBKNP_GGMass->Add(hMcNP_GGMass);

   TH1F* hMcALL_GGMass = (TH1F*)hMcBKNP_GGMass->Clone();
   hMcALL_GGMass->Add(hMcPR_GGMass);

   hMcBK_GGMass->SetFillColor(kRed);
   // hMcBK_GGMass->SetLineColor(kRed);
   hMcBK_GGMass->SetFillStyle(1001);

   hMcBKNP_GGMass->SetFillColor(kBlue);
   // hMcBKNP_GGMass->SetLineColor(kBlue);
   hMcBKNP_GGMass->SetFillStyle(1001);

   hMcALL_GGMass->SetFillColor(kGreen);
   // hMcALL_GGMass->SetLineColor(kGreen);
   hMcALL_GGMass->SetFillStyle(1001);

   c1->cd(1);
   hMcALL_GGMass->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
   hMcALL_GGMass->GetYaxis()->SetTitle("Events/0.01 GeV/c^{2}");
   hMcALL_GGMass->GetYaxis()->SetTitleOffset(1.5);
   hMcALL_GGMass->Draw("");
   hMcBKNP_GGMass->Draw("SAME");
   hMcBK_GGMass->Draw("SAME");

   leg = new TLegend(0.70,0.45,0.95,0.65);
   leg->AddEntry(hMcALL_GGMass,"prompt J/#psi","f");
   leg->AddEntry(hMcBKNP_GGMass,"b #rightarrow J/#psi","f");
   leg->AddEntry(hMcBK_GGMass,"QCD background","f");
   leg->Draw("same");

   TH1F* hMcPR_GTMass = (TH1F*)f.Get("hMcPR_GTMass");
   TH1F* hMcNP_GTMass = (TH1F*)f.Get("hMcNP_GTMass");
   TH1F* hMcBK_GTMass = (TH1F*)f.Get("hMcBK_GTMass");

   TH1F* hMcBKNP_GTMass = (TH1F*)hMcBK_GTMass->Clone();
   hMcBKNP_GTMass->Add(hMcNP_GTMass);

   TH1F* hMcALL_GTMass = (TH1F*)hMcBKNP_GTMass->Clone();
   hMcALL_GTMass->Add(hMcPR_GTMass);

   hMcBK_GTMass->SetFillColor(kRed);
   // hMcBK_GTMass->SetLineColor(kRed);
   hMcBK_GTMass->SetFillStyle(1001);

   hMcBKNP_GTMass->SetFillColor(kBlue);
   // hMcBKNP_GTMass->SetLineColor(kBlue);
   hMcBKNP_GTMass->SetFillStyle(1001);

   hMcALL_GTMass->SetFillColor(kGreen);
   // hMcALL_GTMass->SetLineColor(kGreen);
   hMcALL_GTMass->SetFillStyle(1001);

   c1->cd(2);
   hMcALL_GTMass->GetXaxis()->SetTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
   hMcALL_GTMass->GetYaxis()->SetTitle("Events/0.01 GeV/c^{2}");
   hMcALL_GTMass->GetYaxis()->SetTitleOffset(1.5);
   hMcALL_GTMass->Draw("");
   hMcBKNP_GTMass->Draw("SAME");
   hMcBK_GTMass->Draw("SAME");

   leg = new TLegend(0.70,0.45,0.95,0.65);
   leg->AddEntry(hMcALL_GTMass,"prompt J/#psi","f");
   leg->AddEntry(hMcBKNP_GTMass,"b #rightarrow J/#psi","f");
   leg->AddEntry(hMcBK_GTMass,"QCD background","f");
   leg->Draw("same");

   c1->SaveAs("jpsimassMC.pdf");

   TCanvas *c2 = new TCanvas("c2","c2",200,10,800,400);
   c2->Divide(2,1);

   TH1F* hMcPR_GGLife = (TH1F*)f.Get("hMcPR_GGLife");
   TH1F* hMcNP_GGLife = (TH1F*)f.Get("hMcNP_GGLife");
   // TH1F* hMcBK_GGLife = (TH1F*)f.Get("hMcBK_GGLife");

   // TH1F* hMcBKNP_GGLife = (TH1F*)hMcBK_GGLife->Clone();
   // hMcBKNP_GGLife->Add(hMcNP_GGLife);

   // TH1F* hMcALL_GGLife = (TH1F*)hMcBKNP_GGLife->Clone();
   // hMcALL_GGLife->Add(hMcPR_GGLife);

   // hMcBK_GGLife->SetFillColor(kRed);
   // hMcBK_GGLife->SetLineColor(kRed);
   // hMcBK_GGLife->SetFillStyle(1001);

   // hMcBKNP_GGLife->SetFillColor(kBlue);
   hMcNP_GGLife->SetLineColor(kBlue);
   hMcNP_GGLife->SetLineStyle(kDashed);
   hMcNP_GGLife->SetLineWidth(2);
   // hMcBKNP_GGLife->SetFillStyle(1001);

   hMcPR_GGLife->SetLineColor(kRed);
   hMcPR_GGLife->SetLineStyle(kDotted);
   hMcPR_GGLife->SetLineWidth(2);
   // hMcALL_GGLife->SetLineColor(kGreen);
   // hMcALL_GGLife->SetFillStyle(1001);

   c2->cd(1);
   c2_1->SetLogy(1);
   hMcPR_GGLife->GetXaxis()->SetTitle("J/#psi c#tau [mm]");
   hMcPR_GGLife->GetYaxis()->SetTitle("Events/0.05 mm");
   hMcPR_GGLife->Draw("");
   hMcNP_GGLife->Draw("SAME");
   // hMcBK_GGLife->Draw("SAME");

   leg = new TLegend(0.70,0.55,0.95,0.70);
   leg->AddEntry(hMcPR_GGLife,"prompt J/#psi","l");
   leg->AddEntry(hMcNP_GGLife,"b #rightarrow J/#psi","l");
   // leg->AddEntry(hMcBK_GGLife,"QCD background","f");
   leg->Draw("same");

   TH1F* hMcPR_GTLife = (TH1F*)f.Get("hMcPR_GTLife");
   TH1F* hMcNP_GTLife = (TH1F*)f.Get("hMcNP_GTLife");
   // TH1F* hMcBK_GTLife = (TH1F*)f.Get("hMcBK_GTLife");

   // TH1F* hMcBKNP_GTLife = (TH1F*)hMcBK_GTLife->Clone();
   // hMcBKNP_GTLife->Add(hMcNP_GTLife);

   // TH1F* hMcALL_GTLife = (TH1F*)hMcBKNP_GTLife->Clone();
   // hMcALL_GTLife->Add(hMcPR_GTLife);

   // hMcBK_GTLife->SetFillColor(kRed);
   // hMcBK_GTLife->SetLineColor(kRed);
   // hMcBK_GTLife->SetFillStyle(1001);

   // hMcBKNP_GTLife->SetFillColor(kBlue);
   hMcNP_GTLife->SetLineColor(kBlue);
   hMcNP_GTLife->SetLineStyle(kDashed);
   hMcNP_GTLife->SetLineWidth(2);
   // hMcBKNP_GTLife->SetFillStyle(1001);

   hMcPR_GTLife->SetLineColor(kRed);
   hMcPR_GTLife->SetLineStyle(kDotted);
   hMcPR_GTLife->SetLineWidth(2);
   // hMcALL_GTLife->SetLineColor(kGreen);
   // hMcALL_GTLife->SetFillStyle(1001);

   c2->cd(2);
   c2_2->SetLogy(1);
   hMcPR_GTLife->GetXaxis()->SetTitle("J/#psi c#tau [mm]");
   hMcPR_GTLife->GetYaxis()->SetTitle("Events/0.05 mm");
   hMcPR_GTLife->Draw("");
   hMcNP_GTLife->Draw("SAME");
   // hMcBK_GTLife->Draw("SAME");

   leg = new TLegend(0.70,0.55,0.95,0.70);
   leg->AddEntry(hMcPR_GTLife,"prompt J/#psi","l");
   leg->AddEntry(hMcNP_GTLife,"b #rightarrow J/#psi","l");
   // leg->AddEntry(hMcBK_GTLife,"QCD background","f");
   leg->Draw("same");

   c2->SaveAs("jpsilifeMC.pdf");
}


