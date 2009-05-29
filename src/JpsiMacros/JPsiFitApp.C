// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "fitBfrac.cc"

int main(int argc, char* argv[]) {

  char inputFileName[150];
  if ( argc < 2 ){
    std::cout << "missing argument: insert inputFile with list of root files" << std::endl; 
    return 1;
  }
  strcpy(inputFileName,argv[1]);
  
  TChain *theChain = new TChain("T1");
  char Buffer[500];
  char MyRootFile[2000];
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);

  int nfiles=1;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer))) { 
      sscanf(Buffer,"%s",MyRootFile);
      theChain->Add(MyRootFile);
      std::cout << "chaining " << MyRootFile << std::endl;
      nfiles++;
    }
  }
  inputFile->close();
  delete inputFile;

  fitBfrac jpsiStudy(theChain);
  jpsiStudy.Loop();
  
  return 0;

}
