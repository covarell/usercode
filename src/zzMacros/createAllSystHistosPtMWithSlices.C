void createAllSystHistosPtMWithSlices(TString mass = "125")
{
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",0,true,false,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",-2,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",-3,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",-5,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",-4,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",1,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",3,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",4,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",5,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",6,true,false)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",7,true,false)"); 

  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",0,true,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",-2,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",-3,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",-5,true,true)");  
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",-4,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",1,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",3,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",4,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",5,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",6,true,true)");
  gROOT->ProcessLine(".x studyPtSyst.C(" + mass + ",7,true,true)");
}
