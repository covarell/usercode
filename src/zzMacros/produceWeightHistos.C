void produceWeightHistos(string toReweight = "",
			 string toBeReweighted = "",
			 string result = "") {

  TFile tbr(toBeReweighted.c_str());
 
  TH1F* h1 = (TH1F*)((TH2F*)tbr.Get("Pt_sig"))->ProjectionY("h1");
  h1->SetLineColor(2);
  h1->Sumw2();
  h1->Scale(1./h1->Integral());

  TFile tr(toReweight.c_str());
 
  TH1F* h2 = (TH1F*)((TH2F*)tr.Get("Pt_sig"))->ProjectionY("h2");
  h2->Sumw2();
  h2->Scale(1./h2->Integral());

  TH1F* wei = (TH1F*)h2->Clone();
  wei->SetName("wei");
  wei->Divide(h2,h1);

  TFile f(result.c_str(),"RECREATE");
  h1->Write();
  h2->Write();
  wei->Write();
  f.Close();

}
