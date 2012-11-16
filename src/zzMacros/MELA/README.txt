==================
version notes:  ||
==================

V00-00-01 - 

first stable version.  Low mass LD only.  

-------------------

V00-00-02 - 

corrected bug in calculateAngles() and added protection against 
possibility of m2>m1 (pdfs implicitly assume m1>m2 below threshold)

LD included for highmass (>180), PDFs still only generated for 
100<mZZ<180.

-------------------

V00-00-03 - 

8D template PDF for qq->ZZ background has been extended to 
80<mZZ<185

PDF can now be generated for arbitrary values of mZZ by setting 
global variables mZZbins, lowMzz, and highMzz.  Low m2 cut can 
be set with lowM2.  

-------------------

V00-00-04 - 

Adding script for generating 2d PDF.  Instead of hard coding 
template as 2d array, PDF is built as a prod PDF between 
1D mZZ shape (RooMzzBkg) and 2D RooHistPdf.  

Add_pT_y

Adding possibility to include pT and rapidity in MELA 
evaluation (aka nloMELA)

==================

Content:
--------
PDFs:
RooXZsZs_5D.cxx  RooXZsZs_5D.h 
-> 5D PDF for signal

RooRapidityBkg.cxx RooRapidityBkg.h
-> PDF for bkg Rapidity

RooRapiditySig.cxx RooRapiditySig.h
-> PDF for sig Rapidity

RooTsallis(Exp).cxx RooTsallis(Exp).h
-> PDF for pT

src:
AngularPdfFactory.cc  
-> Utility class to initialize properly the 5D signal PDF

datafiles:
my8DTemplateNotNorm.root
-> 8D template PDF for qq->ZZ background (m2>4 GeV, mZZ 80-185 GeV)

allParams*.txt
-> parameters for signal/background pT description

scripts:
MELA.C  
-> script to evaluate MELA on a given sample (root tree)

Instructions:
------------
root 
 gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L ../PDFs/RooXZsZs_5D.cxx+
.L ../src/AngularPdfFactory.cc+
.L ../PDFs/RooqqZZ_JHU.cxx+
.L ../PDFs/RooTsallis.cxx+              <-- ***RC***
.L ../PDFs/RooTsallisExp.cxx+           <-- ***RC***
.L ../PDFs/RooRapiditySig.cxx+	       <-- ***CM&CY***
.L ../PDFs/RooRapidityBkg.cxx+	       <-- ***CM&CY***
.L MELA.C+

addDtoTree("nameOfYourFile",LHCenergy)

Where nameOfYourFiles.root is the name of the file which contain
your tree. A new file will be created containing
a new tree where the value of LD has been added. LHC energy is either
7 or 8 (TeV).

Notice the tree format should be:
  sigTree->SetBranchAddress("Z1Mass",&m1); 
  sigTree->SetBranchAddress("Z2Mass",&m2);
  sigTree->SetBranchAddress("ZZMass",&mzz);
  sigTree->SetBranchAddress("helcosthetaZ1",&h1);
  sigTree->SetBranchAddress("helcosthetaZ2",&h2);
  sigTree->SetBranchAddress("costhetastar",&hs);
  sigTree->SetBranchAddress("helphi",&phi);
  sigTree->SetBranchAddress("phistarZ1",&phi1);
  sigTree->SetBranchAddress("MC_weight",&w);

  if (containsPt) sigTree->SetBranchAddress("ZZPt",&ZZPt);
  if (containsY) sigTree->SetBranchAddress("ZZRapidity", &ZZY);
  [see below]

Or you may change the macro addDtoTree to adapt to the
format of your tree

note: by default this function only write events which pass the 
following (loose) cuts:

80<mZZ<1000
mz2>4

More options are available:

addDtoTree("nameOfYourFile",LHCenergy,<minimumMzz>,<maximumMzz>)

will only consider events in the range minimumMzz<mZZ<maximumMzz, overwriting
the default. 

will write the original MELA discriminant to the file.

addDtoTree("nameOfYourFile",LHCenergy,<minimumMzz>,<maximumMzz>,true)

will write discriminant to the output file called "melaLDWithPt" that
also includes pT variable in it.

addDtoTree("nameOfYourFile",LHCenergy,<minimumMzz>,<maximumMzz>,false,true)

will write a discriminant to the output file called "melaLDWithY" 
that includes Y variable in it.

addDtoTree("nameOfYourFile",LHCenergy,<minimumMzz>,<maximumMzz>,true,true)

will write a discriminant to the output file called "melaLDWithPtY" 
that includes Y  and pT variable in it.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MELA.C contains also a macro ("calculateAngles") which computes the 5 angles starting from the 4-momenta of the Higgs, the 2 Zeds and the 4 Leptons

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generating 2D PDF - 

MELA.C can also be used to generate a 2D pdf: mZZ vs D.
The mZZ projection is hard coded in RooMELAModel*.tpl. 
This can be changed to your favorite mZZ projection.  The
D portion is configured such that the projection onto 
mZZ is exactly the function you input. Ranges to be used 
for mZZ and mZ2 can be set via the following global variables:

mZZbins
lowMzz
highMzz
lowM2

To generate signal and background PDFs:

root
 gSystem->AddIncludePath("-I/$ROOFITSYS/include/");
.L ../PDFs/RooXZsZs_5D.cxx+
.L ../src/AngularPdfFactory.cc+
.L ../PDFs/RooqqZZ_JHU.cxx+
.L MELA.C+
//store D vs mZZ template for signal
storeLDDistribution(true,"mySignalFile.root")  
//store D vs mZZ template for background
storeLDDistribution(false,"myBackgroundFile.root")  
//- - - - - CODA - - - -
genMELApdf(true)
genMELApdf(false)

/* Note the input file names are passed to 
TChain::Add(), so wild cards can be used.  The
tree name is assumed to be angles but this can
be changed in LDDistributionSignal and 
LDDistributionBackground */

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Alternate method to generate 2D PDF - 

Follow instructions above until "CODA", then skip
to here. 

Once Dsignal*.root and Dbackground*.root are 
created, one can build a fully correlated 2D PDF
using Dsignal*.root and Dbackground*.root to build
a 2D RooHistPdf then multiplying it by a 1D PDF which
represent the desired mZZ project.  An example is 
shown in buil2dPdf.C where the background mZZ 
projection is taken from a custom RooAbsPdf, 
RooMzzBkg.cc, and the signal mZZ projection is 
a sum of 2 gaussians.  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculating the angles and PDFs also works in CMSSW
Building the 2dPDF is not yet supported

cvs co -d JHU/MELA/src UserCode/JHU/MELA
rm JHU/MELA/src/scripts/build2dPdf.C
scramv1 b

the package that uses the resulting library needs to have
an entry:
<use name="JHUMELA"/>
in its BuildFile.xml

