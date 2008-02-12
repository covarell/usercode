#ifndef Alignment_HIPAlignmentAlgorithm_HIPAlignmentAlgorithm_h
#define Alignment_HIPAlignmentAlgorithm_HIPAlignmentAlgorithm_h

#include "Alignment/CommonAlignment/interface/AlignableObjectId.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmBase.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentIORoot.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

#include "Alignment/CommonAlignmentProducer/interface/TrackLocalAngle.h"
 
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

class AlignableNavigator;
class TFile;
class TTree;

class HIPAlignmentAlgorithm : public AlignmentAlgorithmBase
{

 public:
  
  /// Constructor
  HIPAlignmentAlgorithm(const edm::ParameterSet& cfg);

  /// Destructor
  ~HIPAlignmentAlgorithm() {};

  /// Call at beginning of job
  void initialize( const edm::EventSetup& setup, 
                   AlignableTracker* tracker, AlignableMuon* muon, 
                   AlignmentParameterStore* store);

  /// Call at end of job
  void terminate(void);

  /// Called at start of new loop
  void startNewLoop(void);

  /// Run the algorithm on trajectories and tracks
  void run( const edm::EventSetup& setup, 
	    const ConstTrajTrackPairCollection& tracks,
	    const edm::SimTrackContainer& simcoll);

 private:

  // private member functions

  int readIterationFile(std::string filename);
  void writeIterationFile(std::string filename,int iter);
  void setAlignmentPositionError(void);
  double calcAPE(double* par, int iter,std::string param);
  void bookRoot(void);
  void fillRoot(void);
  bool calcParameters(Alignable* ali);
  void collector(void);
  // COSMICS
  std::vector<unsigned int> storeNumberingScheme(const DetId& detid, int type);

  // private data members

  AlignmentParameterStore* theAlignmentParameterStore;
  std::vector<Alignable*> theAlignables;
  AlignableNavigator* theAlignableDetAccessor;

  AlignmentIORoot    theIO;
  int ioerr;
  int theIteration;
  // COSMICS
  TrackLocalAngle* theAngleFinder;

  // steering parameters

  // verbosity flag
  bool verbose;
  // names of IO root files
  std::string outfile,outfile2,outpath,suvarfile,sparameterfile;
  std::string struefile,smisalignedfile,salignedfile,siterationfile;
  // alignment position error parameters
  double apesp[3],aperp[3];
  std::string apeparam;
  // max allowed pull (residual / uncertainty) on a hit used in alignment
  double theMaxAllowedHitPull;
  // min number of hits on alignable to calc parameters
  int theMinimumNumberOfHits;
  // max allowed rel error on parameter (else not used)
  double theMaxRelParameterError;
  // collector mode (parallel processing)
  bool isCollector;
  int theCollectorNJobs;
  std::string theCollectorPath;
  int theEventPrescale,theCurrentPrescale;

  std::vector<AlignableObjectId::AlignableObjectIdType> theLevels; // for survey residuals

  // root tree variables
  TFile* theFile;
  TTree* theTree; // event-wise tree
  TFile* theFile2;
  TTree* theTree2; // alignable-wise tree

  // variables for event-wise tree
  // COSMICS
  static const int MAXSIM = 5;
  static const int MAXREC = 5;
  static const int MAXHIT = 99*MAXREC;

  // int m_Run,m_Event;
  int allTracks;
  int m_NtracksSim;
  float m_PtSim[MAXSIM],m_EtaSim[MAXSIM],m_PhiSim[MAXSIM];
  int m_Ntracks,m_allHits,m_Nhits[MAXREC];
  float m_Pt[MAXREC],m_Eta[MAXREC],m_Phi[MAXREC],m_Chi2n[MAXREC];
  int m_isOnAli[MAXHIT];
  int m_hType[MAXHIT], m_hLayer[MAXHIT], m_hOwnerTrack[MAXHIT];
  int m_hFwBw[MAXHIT], m_hIntExt[MAXHIT], m_hStrRod[MAXHIT], m_hModule[MAXHIT];
  float m_hR[MAXHIT],m_hPhi[MAXHIT],m_hZ[MAXHIT];
  float m_hLocalX[MAXHIT],m_hLocalY[MAXHIT],m_hLocalZ[MAXHIT];
  float m_hLocalAngleMono[MAXHIT], m_hLocalAngleSter[MAXHIT], m_hChargeMono[MAXHIT], m_hChargeSter[MAXHIT], m_hBarMono[MAXHIT], m_hBarSter[MAXHIT];
  float m_Xres[MAXHIT],m_Yres[MAXHIT],m_Xerr[MAXHIT],m_Yerr[MAXHIT];
  float m_XerrHit[MAXHIT],m_XerrIP[MAXHIT];

  // variables for alignable-wise tree
  int m2_Nhit,m2_Type,m2_Layer;
  float m2_Xpos, m2_Ypos, m2_Zpos, m2_Eta, m2_Phi; 
  int m2_Id,m2_ObjId;

};

#endif
