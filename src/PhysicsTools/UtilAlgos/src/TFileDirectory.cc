#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "TROOT.h"
#include "TFile.h"
#include "FWCore/Utilities/interface/Exception.h"
using namespace std;

void TFileDirectory::cd() const {
  string fpath = fullPath();
  TDirectory * dir = file_->GetDirectory( fpath.c_str() );
  if ( dir == 0 ) {
    if ( ! path_.empty() ) {
      dir = file_->GetDirectory( path_.c_str() );
      if ( dir == 0 )
	throw 
	  cms::Exception( "InvalidDirectory" ) 
	    << "Can't change directory to path: " << path_;
    } else {
      dir = file_;
    }
    dir = dir->mkdir( dir_.c_str(), descr_.c_str() );  
    if ( dir == 0 )
      throw 
	cms::Exception( "InvalidDirectory" ) 
	  << "Can't create directory " << dir_ << " in path: " << path_;
  }
  bool ok = file_->cd( fpath.c_str() );
  if ( ! ok )
    throw 
      cms::Exception( "InvalidDirectory" ) 
	<< "Can't change directory to path: " << fpath;
}

std::string TFileDirectory::fullPath() const {
  return string( path_.empty() ? dir_ : path_ + "/" + dir_ );
}

TFileDirectory TFileDirectory::mkdir( const std::string & dir, const std::string & descr ) {
  cd();
  return TFileDirectory( dir, descr, file_, fullPath() );
}
