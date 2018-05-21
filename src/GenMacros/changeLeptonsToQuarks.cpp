// g++ -Wall -o changeLeptonsToQuarks changeLeptonsToQuarks.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

bool replacce(std::string& str, const std::string& from, const std::string& to) {
  // replace substring inside string
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}

std::string randFlav() {
  // choose Z->qq decays according to BRs
  std::string quarkType;
  float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
  if (r < 0.22286) quarkType = "1    1";       // d
  else if (r < 0.4457) quarkType = "3    1";   // s
  else if (r < 0.6686) quarkType = "5    1";   // b
  else if (r < 0.8342) quarkType = "2    1";   // u
  else quarkType = "4    1";                   // c
  return quarkType;
}

int main(int argc, char** argv)
{
  // set Number of lhe files
  static const int nLHE = 500;
  // set number of events per job
  static const int nEvents = 1000;

  if(argc < 2) {
    std::cout << ">>>changeLeptonsToQuarks.cpp::Usage:   " << argv[0] << "   initialFile.lhe " << std::endl;
    return -1;
  }

  srand (time(NULL));   // initialize random seed
  
  char* initialFileName = argv[1]; 
  std::cout << "initialFileName = " << initialFileName << std::endl;
  
  // open lhe files
  char outputName[20];
  std::ifstream initialFile(initialFileName, std::ios::in);
  std::ofstream *outFile[nLHE];
  for (int i= 0; i < nLHE; i++) {
    sprintf(outputName, "lheFiles/newFile%d.lhe",i);
    outFile[i] = new std::ofstream(outputName, std::ios::out);
  }
   
  int ievent = 0;
  int whichFile = 0;
  int whichFilePrevEvent = 0;
  std::string line, line2;
  std::string leptType = "11    1";
  std::string leptTypeMinus = "-11    1";
  std::string quarkType = randFlav();

  while(!initialFile.eof()) {
    
    getline(initialFile, line);
    if( !initialFile.good() ) break;   
    
    whichFile = ievent/nEvents;
    if (whichFile > whichFilePrevEvent) {
      *outFile[whichFilePrevEvent] << "</LesHouchesEvents>" << std::endl;
      outFile[whichFilePrevEvent]->close();
      std::ifstream tempinitialFile(initialFileName, std::ios::in);
      int iline2 = 0; 
      while(!tempinitialFile.eof()) {
    	getline(tempinitialFile, line2);
        *outFile[whichFile] << line2 << std::endl;
        iline2++;
	if( iline2 == 5 ) break;
      }
      tempinitialFile.close();
    }
    whichFilePrevEvent = ievent/nEvents;

    if( line.find(leptType) != std::string::npos) {
      // change
      if( line.find(leptTypeMinus) != std::string::npos) replacce(line,"    0    0 ","    0   503 ");
      else replacce(line,"    0    0 ","   503    0 ");   //color
      replacce(line, leptType, quarkType);   //flavor        
    }
      
    *outFile[whichFile] << line << std::endl;

    if (ievent == nLHE*nEvents - 1) break;
    if( line == "</event>" ) {   // new event, set new values
      ievent++;  
      quarkType = randFlav();
      if (leptType == "11    1") {
	leptType = "13    1";
	leptTypeMinus = "-13    1";
      } else {
	leptType = "11    1";
	leptTypeMinus = "-11    1";  
      }
    }

  }
  
  std::cout << "Done " << std::endl;
  return 0;
}
