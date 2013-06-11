void createAllTemplatesPtMWithSlices()
{
  gROOT->ProcessLine(".L makeTemplatesPtSyst.C+");  

  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",0,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-2,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-3,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-5,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-4,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-6,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-7,true,false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",1,true,false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",3,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",4,true,false)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",5,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",6,true,false)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",7,true,false)"); 

  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",0,true,true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-2,true,true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-3,true,true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-5,true,true)");  
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-4,true,true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-6,true,true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",-7,true,true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",1,true,true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",3,true,true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",4,true,true)");
  // gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",5,true,true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",6,true,true)");
  gROOT->ProcessLine("makeTemplatesPtSyst(\"2e2mu\",7,true,true)");
}
