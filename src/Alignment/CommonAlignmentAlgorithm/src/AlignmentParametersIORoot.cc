// this class's header
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParametersIORoot.h"

#include <TTree.h>

#include "Alignment/CommonAlignmentParametrization/interface/CompositeRigidBodyAlignmentParameters.h"
#include "Alignment/CommonAlignmentParametrization/interface/RigidBodyAlignmentParameters.h"

#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"



// ----------------------------------------------------------------------------
// constructor

AlignmentParametersIORoot::AlignmentParametersIORoot()
{
  treename = "AlignmentParameters";
  treetxt = "Alignment Parameters";
}

// ----------------------------------------------------------------------------

void AlignmentParametersIORoot::createBranches(void) 
{
  tree->Branch("parSize",   &theCovRang,   "CovRang/I");
  tree->Branch("Id",        &theId,        "Id/I");
  tree->Branch("Par",       &thePar,       "Par[CovRang]/D");
  tree->Branch("covarSize", &theCovarRang, "CovarRang/I");
  tree->Branch("Cov",       &theCov,       "Cov[CovarRang]/D");
  tree->Branch("ObjId",     &theObjId,     "ObjId/I");
}

// ----------------------------------------------------------------------------

void AlignmentParametersIORoot::setBranchAddresses(void) 
{
  tree->SetBranchAddress("parSize",   &theCovRang);
  tree->SetBranchAddress("covarSize", &theCovarRang);
  tree->SetBranchAddress("Id",        &theId);
  tree->SetBranchAddress("Par",       &thePar);
  tree->SetBranchAddress("Cov",       &theCov);
  tree->SetBranchAddress("ObjId",     &theObjId);
}

// ----------------------------------------------------------------------------
// find tree entry based on detID and typeID

int AlignmentParametersIORoot::findEntry(unsigned int detId,int comp)
{
  double noAliPar = tree->GetEntries();
  for (int ev = 0;ev<noAliPar;ev++) {
    tree->GetEntry(ev); 
    if ( theId==detId && theObjId==comp ) return (ev);
  }
  return -1;
}

// ----------------------------------------------------------------------------
int AlignmentParametersIORoot::writeOne(Alignable* ali)
{
  AlignmentParameters* ap =ali->alignmentParameters();
  AlgebraicVector params  = ap->parameters();
  AlgebraicSymMatrix cov  = ap->covariance();

  theCovRang   = params.num_row();
  theCovarRang = theCovRang*(theCovRang+1)/2;
  int count=0;
  for(int row=0;row<theCovRang;row++){
    thePar[row]=params[row];
    for(int col=0;col<theCovRang;col++){
      if(row-1<col) { theCov[count] = cov[row][col]; count++; }
    }
  }

  TrackerAlignableId converter;
  theId = converter.alignableId(ali);
  theObjId = converter.alignableTypeId(ali);

  tree->Fill();
  return 0;
}

// ----------------------------------------------------------------------------

AlignmentParameters* AlignmentParametersIORoot::readOne( Alignable* ali, int& ierr )
{
  
  TrackerAlignableId converter;
  int obj = converter.alignableTypeId(ali);
  unsigned int detInt = converter.alignableId(ali);

  AlignmentParameters* alipar;
  AlgebraicVector par(nParMax,0);
  AlgebraicSymMatrix cov(nParMax,0);
  std::vector<bool> sel = ali->alignmentParameters()->selector();
 
  int entry = findEntry(detInt,obj);
  if( entry != -1 ) 
	{
	  tree->GetEntry(entry);
	  int covsize = theCovRang;
	  int count=0;
	  for(int row=0;row<covsize;row++) 
		{
		  par[row]=thePar[row];
		  for(int col=0; col < covsize;col++) {
			if(row-1<col) {cov[row][col]=theCov[count];count++;}
		  }
		} 
	  if (obj == 0 ) // create parameters for sensor
		alipar = new RigidBodyAlignmentParameters(ali,par,cov,sel);
	  else // create parameters for composed object
		alipar = new CompositeRigidBodyAlignmentParameters(ali,par,cov,sel);
	  alipar->setValid(true); 
	  ierr=0;
	  return alipar;
	}

  ierr=-1;
  return(0);
}


