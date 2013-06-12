void createAllTemplatesPtMWithSlices(TString type = "0")
{
  gROOT->ProcessLine(".L makeTemplatesPtSyst.C+");  

  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",0," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-2," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-3," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-5," + type + ",true)");  
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-4," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-6," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-7," + type + ",true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",1," + type + ",true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",3," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",4," + type + ",true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",5," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",6," + type + ",true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",7," + type + ",true)");

  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",0," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-2," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-3," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-5," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-4," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-6," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-7," + type + ",false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",1," + type + ",false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",3," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",4," + type + ",false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",5," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",6," + type + ",false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",7," + type + ",false)"); 

}
