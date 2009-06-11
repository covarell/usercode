// C++ includes
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "MakeDataSet.cc"

int main(int argc, char* argv[]) {

  char inputFileName[150];
  if ( argc < 2 ){
    std::cout << "missing argument: insert inputFile with list of root files" << std::endl; 
    return 1;
  }

 TChain *theChain=new TChain("T1");

  //char* output=argv[1];
  for(int i=1;i<argc;i++) {
    theChain->Add(argv[i]);
    //printf("Processing:%s\n",argv[i]);
  }

  MakeDataSet jpsiStudy(theChain);
  jpsiStudy.Loop();
  
  return 0;

}
