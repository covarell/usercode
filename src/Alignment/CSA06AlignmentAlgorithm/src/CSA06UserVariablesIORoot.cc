#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"
#include "Alignment/CSA06AlignmentAlgorithm/interface/CSA06UserVariables.h"


// this class's header
#include "Alignment/CSA06AlignmentAlgorithm/interface/CSA06UserVariablesIORoot.h"

// ----------------------------------------------------------------------------
// constructor

CSA06UserVariablesIORoot::CSA06UserVariablesIORoot()
{
  treename = "T9";
  treetxt = "CSA06 User Variables";
}

// ----------------------------------------------------------------------------

void CSA06UserVariablesIORoot::createBranches(void) 
{
  tree->Branch("Id",        &Id,        "Id/i");
  tree->Branch("ObjId",     &ObjId,     "ObjId/I");

  tree->Branch("Nhit",      &Nhit,      "Nhit/I");
  tree->Branch("Nparj",     &Nparj,     "Nparj/I");
  tree->Branch("Jtvj",      &Jtvj,      "Jtvj[Nparj]/D");
  tree->Branch("Npare",     &Npare,     "Npare/I");
  tree->Branch("Jtve",      &Jtve,      "Jtve[Npare]/D");
}

// ----------------------------------------------------------------------------

void CSA06UserVariablesIORoot::setBranchAddresses(void) 
{
  tree->SetBranchAddress("Id",        &Id);
  tree->SetBranchAddress("ObjId",     &ObjId);

  tree->SetBranchAddress("Nhit",      &Nhit);
  tree->SetBranchAddress("Nparj",     &Nparj);
  tree->SetBranchAddress("Jtvj",      &Jtvj);
  tree->SetBranchAddress("Npare",     &Npare);
  tree->SetBranchAddress("Jtve",      &Jtve);
}

// ----------------------------------------------------------------------------
// find tree entry based on detID and typeID

int CSA06UserVariablesIORoot::findEntry(unsigned int detId,int comp)
{
  if (newopen) { // we're here first time
    edm::LogInfo("Alignment") <<"[CSA06UserVariablesIORoot::findEntry] fill map ...";
    treemap.erase(treemap.begin(),treemap.end());
    for (int ev = 0;ev<tree->GetEntries();ev++) {
      tree->GetEntry(ev); 
      treemap[std::make_pair(Id,ObjId)]=ev;
    }
    newopen=false;
  }
  
  // now we have filled the map
  treemaptype::iterator imap = treemap.find(std::make_pair(detId,comp));
  int result=-1;
  if (imap != treemap.end()) result=(*imap).second;
  return result;


  //double noAliPar = tree->GetEntries();
  //for (int ev = 0;ev<noAliPar;ev++) {
  //  tree->GetEntry(ev); 
  //  if(Id==detId&&comp==ObjId) return (ev);
  //}
  //return(-1);
}

// ----------------------------------------------------------------------------

int CSA06UserVariablesIORoot::writeOne(Alignable* ali)
{
  AlignmentParameters* ap=ali->alignmentParameters();

  if ((ap->userVariables())==0) { 
    edm::LogError("Alignment") <<"UserVariables not found!"; 
    return -1; 
  }

  CSA06UserVariables* uvar = 
    dynamic_cast<CSA06UserVariables*>(ap->userVariables());

  AlgebraicSymMatrix jtvj = uvar->jtvj;
  AlgebraicVector jtve = uvar->jtve;
  int nhit=uvar->nhit;
  int np=jtve.num_row();

  TrackerAlignableId ID;
  unsigned int detInt = ID.alignableId(ali);
  int typeInt = ID.alignableTypeId(ali);

  Nhit=nhit;
  Npare=np;
  Nparj=np*(np+1)/2;
  int count=0;
  for(int row=0;row<np;row++){
    Jtve[row]=jtve[row];
    for(int col=0;col<np;col++){
      if(row-1<col){Jtvj[count]=jtvj[row][col];count++;}
    }
  }
  Id = detInt;
  ObjId = typeInt;

  tree->Fill();
  return 0;
}

// ----------------------------------------------------------------------------

AlignmentUserVariables* CSA06UserVariablesIORoot::readOne(Alignable* ali, 
  int& ierr)
{
  ierr=0;
  CSA06UserVariables* uvar;

  TrackerAlignableId ID;
  int obj = ID.alignableTypeId(ali);
  unsigned int detInt = ID.alignableId(ali);
  int entry = findEntry(detInt,obj);
  if(entry!=-1) {
    tree->GetEntry(entry);

    int np=Npare;
    AlgebraicVector jtve(np,0);
    AlgebraicSymMatrix jtvj(np,0);
    int count=0;
    for(int row=0;row<np;row++) {
      jtve[row]=Jtve[row];
      for(int col=0; col < np;col++) {
 	if(row-1<col) {jtvj[row][col]=Jtvj[count];count++;}
      }
    } 

    uvar = new CSA06UserVariables(np);
    uvar->jtvj=jtvj;
    uvar->jtve=jtve;
    uvar->nhit=Nhit;

    return uvar;
  }

  //  ierr=-1;
  return 0 ;
}

//-----------------------------------------------------------------------------

void 
CSA06UserVariablesIORoot::writeCSA06UserVariables (const Alignables& alivec, 
  const char* filename, int iter, bool validCheck, int& ierr)
{
  ierr=0;
  int iret;
  iret = open(filename,iter,true);
  if (iret!=0) { ierr=-1; return;}
  iret = write(alivec,validCheck);
  if (iret!=0) { ierr=-2; return;}
  iret = close();
  if (iret!=0) { ierr=-3; return;}
}

//-----------------------------------------------------------------------------

std::vector<AlignmentUserVariables*> 
CSA06UserVariablesIORoot::readCSA06UserVariables (const Alignables& alivec, 
  const char* filename, int iter, int& ierr)
{
  std::vector<AlignmentUserVariables*> result;
  ierr=0;
  int iret;
  iret = open(filename,iter,false);
  if (iret!=0) { ierr=-1; return result;}
  result = read(alivec,iret);
  if (iret!=0) { ierr=-2; return result;}
  iret = close();
  if (iret!=0) { ierr=-3; return result;}

  return result;
}
