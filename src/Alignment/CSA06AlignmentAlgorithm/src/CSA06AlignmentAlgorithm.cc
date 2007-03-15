#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"
#include "Alignment/MuonAlignment/interface/AlignableMuon.h"

#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"
#include "Alignment/CSA06AlignmentAlgorithm/interface/CSA06AlignmentAlgorithm.h"
#include "Alignment/CSA06AlignmentAlgorithm/interface/CSA06UserVariables.h"
#include "Alignment/CSA06AlignmentAlgorithm/interface/CSA06UserVariablesIORoot.h"
#include "Alignment/CommonAlignment/interface/AlignableNavigator.h"  

#include "Alignment/CSA06AlignmentAlgorithm/interface/TrajectoryMeasurementResidual.h"

#include <fstream>

using namespace std;

// Constructor ----------------------------------------------------------------

CSA06AlignmentAlgorithm::CSA06AlignmentAlgorithm(const edm::ParameterSet& cfg):
  AlignmentAlgorithmBase( cfg )
{

  // parse parameters

  verbose = cfg.getParameter<bool>("verbosity");

  outpath = cfg.getParameter<string>("outpath");
  outfile = cfg.getParameter<string>("outfile");
  outfile2 = cfg.getParameter<string>("outfile2");
  struefile = cfg.getParameter<string>("trueFile");
  smisalignedfile = cfg.getParameter<string>("misalignedFile");
  salignedfile = cfg.getParameter<string>("alignedFile");
  siterationfile = cfg.getParameter<string>("iterationFile");
  suvarfile = cfg.getParameter<string>("uvarFile");
  sparameterfile = cfg.getParameter<string>("parameterFile");

  outfile        =outpath+outfile;
  outfile2       =outpath+outfile2;
  struefile      =outpath+struefile;
  smisalignedfile=outpath+smisalignedfile;
  salignedfile   =outpath+salignedfile;
  siterationfile =outpath+siterationfile;
  suvarfile      =outpath+suvarfile;
  sparameterfile =outpath+sparameterfile;

  // parameters for APE
  apeparam=cfg.getParameter<string>("apeParam");
  vector<double> vapesp = cfg.getParameter< vector<double> >("apeSPar");
  apesp[0]=vapesp[0];
  apesp[1]=vapesp[1];
  apesp[2]=vapesp[2];
  vector<double> vaperp = cfg.getParameter< vector<double> >("apeRPar");
  aperp[0]=vaperp[0];
  aperp[1]=vaperp[1];
  aperp[2]=vaperp[2];

  theMinimumNumberOfHits = cfg.getParameter<int>("minimumNumberOfHits");
  theMaxRelParameterError = cfg.getParameter<double>("maxRelParameterError");

  // for collector mode (parallel processing)
  isCollector=cfg.getParameter<bool>("collectorActive");  
  theCollectorNJobs=cfg.getParameter<int>("collectorNJobs");
  theCollectorPath=cfg.getParameter<string>("collectorPath");
  if (isCollector) edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm] Collector mode";

  theEventPrescale = cfg.getParameter<int>("eventPrescale");
  theCurrentPrescale = theEventPrescale;

  edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm] constructed.";

}

// Call at beginning of job ---------------------------------------------------

void 
CSA06AlignmentAlgorithm::initialize( const edm::EventSetup& setup, 
                                     AlignableTracker* tracker, AlignableMuon* muon, 
                                     AlignmentParameterStore* store )
{
  edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm] Initializing...";

  allTracks = 0;
  // accessor Det->AlignableDet
  if ( !muon )
    theAlignableDetAccessor = new AlignableNavigator(tracker);
  else if ( !tracker )
    theAlignableDetAccessor = new AlignableNavigator(muon);
  else 
    theAlignableDetAccessor = new AlignableNavigator(tracker,muon);

  // set alignmentParameterStore
  theAlignmentParameterStore=store;

  // get alignables
  theAlignables = theAlignmentParameterStore->alignables();

  // iterate over all alignables and attach user variables
  for( vector<Alignable*>::const_iterator it=theAlignables.begin(); 
       it!=theAlignables.end(); it++ )
    {
      AlignmentParameters* ap = (*it)->alignmentParameters();
      int npar=ap->numSelected();
      CSA06UserVariables* userpar = new CSA06UserVariables(npar);
      ap->setUserVariables(userpar);
    }

  // try to read in alignment parameters from a previous iteration
  AlignablePositions theAlignablePositionsFromFile =
    theIO.readAlignableAbsolutePositions(theAlignables,
   (char*)salignedfile.c_str(),-1,ioerr);

  int numAlignablesFromFile = theAlignablePositionsFromFile.size();

  if (numAlignablesFromFile==0) { // file not there: first iteration 

    // set iteration number to 1
    if (isCollector) theIteration=0;
    else theIteration=1;
    edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm] File not found => iteration "<<theIteration;

    // get true (de-misaligned positions) and write to root file
    // hardcoded iteration=1
    theIO.writeAlignableOriginalPositions(theAlignables,
      (char*)struefile.c_str(),1,false,ioerr);

    // get misaligned positions and write to root file
    // hardcoded iteration=1
    theIO.writeAlignableAbsolutePositions(theAlignables,
      (char*)smisalignedfile.c_str(),1,false,ioerr);

  }

  else { // there have been previous iterations

    edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm] Alignables Read " 
      << numAlignablesFromFile;

    // get iteration number from file     
    theIteration = readIterationFile(siterationfile);

    // increase iteration
    theIteration++;
    edm::LogWarning("Alignment") <<"[CSA06AlignmentAlgorithm] Iteration increased by one!";

    // now apply psotions of file from prev iteration
    edm::LogWarning("Alignment") <<"[CSA06AlignmentAlgorithm] Apply positions from file ...";
    theAlignmentParameterStore->applyAlignableAbsolutePositions(theAlignables, 
      theAlignablePositionsFromFile,ioerr);

  }

  edm::LogWarning("Alignment") <<"[CSA06AlignmentAlgorithm] Current Iteration number: " 
    << theIteration;

  // set alignment position error 
  setAlignmentPositionError();

  // book root trees
  bookRoot();

  // run collector job if we are in parallel mode
  if (isCollector) collector();
  
  //for TrackAngle
  edm::ESHandle<TrackerGeometry> trackerGeom;
  setup.get<TrackerDigiGeometryRecord>().get( trackerGeom );
  theAngleFinder = new TrackLocalAngle( trackerGeom.product() ); 
}

// Call at end of job ---------------------------------------------------------

void CSA06AlignmentAlgorithm::terminate(void)
{

  edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm] Terminating";

  // write user variables
  CSA06UserVariablesIORoot CSA06IO;
  CSA06IO.writeCSA06UserVariables (theAlignables,(char*)suvarfile.c_str(),
    theIteration,false,ioerr);

  // now calculate alignment corrections ...
  int ialigned=0;
  // iterate over alignment parameters
  for(vector<Alignable*>::const_iterator
    it=theAlignables.begin(); it!=theAlignables.end(); it++) {
    Alignable* ali=(*it);
    // Alignment parameters
    AlignmentParameters* par = ali->alignmentParameters();
    // try to calculate parameters
    bool test = calcParameters(ali);
    // if successful, apply parameters
    if (test) { 
      edm::LogInfo("Alignment") << "now apply params";
      theAlignmentParameterStore->applyParameters(ali);
      // set these parameters 'valid'
      ali->alignmentParameters()->setValid(true);
      // increase counter
      ialigned++;
    }
    else par->setValid(false);
  }
  edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm::terminate] Aligned units: " << ialigned;

  // fill alignable wise root tree
  fillRoot();

  edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm] Writing aligned parameters to file: " << theAlignables.size();

  // write new absolute positions to disk
  theIO.writeAlignableAbsolutePositions(theAlignables,
    (char*)salignedfile.c_str(),theIteration,false,ioerr);

  // write alignment parameters to disk
  theIO.writeAlignmentParameters(theAlignables, 
    (char*)sparameterfile.c_str(),theIteration,false,ioerr);

  // write iteration number to file
  writeIterationFile(siterationfile,theIteration);


  // write out trees and close root file

  // eventwise tree
  theFile->cd();
  theTree->Write();
  theFile->Close();
  delete theFile;

  // alignable-wise tree is only filled once
  //if ((!isCollector && theIteration==1)||
      //    ( isCollector && theIteration==2)) { 
  if (theIteration==1) { // only for 1st iteration
    theFile2->cd();
    theTree2->Write(); 
    theFile2->Close();
    delete theFile2;
  }  

}

// Run the algorithm on trajectories and tracks -------------------------------

void CSA06AlignmentAlgorithm::run( const edm::EventSetup& setup, const ConstTrajTrackPairCollection& tracks , const edm::SimTrackContainer& simcoll )
{
  if (isCollector) return;

  TrajectoryStateCombiner tsoscomb;
  TrackerAlignableId id;

  int itr=0;
  int itrsim=0;
  int ahit=0;
  m_Ntracks=0;
  m_NtracksSim=0;
  m_allHits=0;

  theFile->cd();

  float ptsim = -9999.;
  float etasim = -9999.;
  float phisim = -9999.;

   // loop over sim tracks 
  if (&simcoll) {
    for(edm::SimTrackContainer::const_iterator trackCI = simcoll.begin(); 
	trackCI != simcoll.end(); trackCI++) {
      
      ptsim    = trackCI->momentum().perp();
      etasim   = trackCI->momentum().eta();
      phisim   = trackCI->momentum().phi();
      
      // fill sim track parameters in root tree
      if (itrsim<MAXSIM) {
	m_PtSim[itrsim]=ptsim;
	m_EtaSim[itrsim]=etasim;
	m_PhiSim[itrsim]=phisim;
	itrsim++;
	m_NtracksSim=itrsim;
      }
    }

  }

  float pt = -9999.;
  float eta = -9999.;
  float phi = -9999.;
  float chi2n = -9999.;
  int nhit = -9999; 

  if ( tracks.size() > 0 ) {
    
    // loop over tracks  
    for( ConstTrajTrackPairCollection::const_iterator it=tracks.begin();
	 it!=tracks.end();it++) {
      
      const Trajectory* traj = (*it).first;
      const reco::Track* track = (*it).second;
      
      pt    = track->pt();
      eta   = track->eta();
      phi   = track->phi();
      chi2n = track->normalizedChi2();
      nhit  = track->numberOfValidHits(); 
      
      if (verbose) edm::LogInfo("Alignment") << "New track pt,eta,phi,chi2n,hits: " << pt <<","<< eta <<","<< phi <<","<< chi2n << ","<<nhit;
      
      // fill track parameters in root tree
      if (itr<MAXREC) {
	m_Nhits[itr]=nhit;
	m_Pt[itr]=pt;
	m_Eta[itr]=eta;
	m_Phi[itr]=phi;
	m_Chi2n[itr]=chi2n;
	itr++;
	m_Ntracks=itr;
      }
      
      vector<const TransientTrackingRecHit*> hitvec;
      vector<TrajectoryStateOnSurface> tsosvec;
      vector<TrajectoryMeasurement> tmvec;
    
      // loop over measurements	
      vector<TrajectoryMeasurement> measurements = traj->measurements();
      for (vector<TrajectoryMeasurement>::iterator im=measurements.begin();
	   im!=measurements.end(); im++) {
	TrajectoryMeasurement meas = *im;
	const TransientTrackingRecHit* hit = &(*meas.recHit());
       
	if (hit->isValid()) {

	  // this is the updated state (including the current hit)
	  //TrajectoryStateOnSurface tsos=meas.updatedState();
	  // combine fwd and bwd predicted state to get state 
	  // which excludes current hit
	  TrajectoryStateOnSurface tsosc = 
	    tsoscomb.combine(meas.forwardPredictedState(),
			     meas.backwardPredictedState());
	  hitvec.push_back(hit);
	  tmvec.push_back(meas);
	  tsosvec.push_back(tsosc);

	}
      }
      
      // transform RecHit vector to AlignableDet vector
      vector <AlignableDet*> alidetvec = 
	theAlignableDetAccessor->alignableDetsFromHits(hitvec);
      
     // get concatenated alignment parameters for list of alignables
      CompositeAlignmentParameters aap = 
	theAlignmentParameterStore->selectParameters(alidetvec);
      
      vector<TrajectoryStateOnSurface>::const_iterator itsos = tsosvec.begin();
      vector<const TransientTrackingRecHit*>::const_iterator ihit = hitvec.begin();
      vector<TrajectoryMeasurement>::const_iterator imea = tmvec.begin();
      
      // loop over vectors(hit,tsos)
      while (itsos != tsosvec.end()) 
	{
     
	  // get trajectory impact point
	  LocalPoint alvec = (*itsos).localPosition();
	  AlgebraicVector pos(2);
	  pos[0]=alvec.x(); // local x
	  pos[1]=alvec.y(); // local y
	  
	  // get impact point covariance
	  AlgebraicSymMatrix ipcovmat(2);
	  ipcovmat[0][0] = (*itsos).localError().positionError().xx();
	  ipcovmat[1][1] = (*itsos).localError().positionError().yy();
	  ipcovmat[0][1] = (*itsos).localError().positionError().xy();
	  
	  // get hit local position and covariance
	  AlgebraicVector coor(2);
	  coor[0] = (*ihit)->localPosition().x();
	  coor[1] = (*ihit)->localPosition().y();
	  
	  AlgebraicSymMatrix covmat(2);
	  covmat[0][0] = (*ihit)->localPositionError().xx();
	  covmat[1][1] = (*ihit)->localPositionError().yy();
	  covmat[0][1] = (*ihit)->localPositionError().xy();
	  
	  // add hit and impact point covariance matrices
	  covmat = covmat + ipcovmat;
	  
	  // get AlignableDet for this hit
	  const GeomDet* det=(*ihit)->det();
	  AlignableDet* alidet = 
	    theAlignableDetAccessor->alignableDetFromGeomDet(det);
	  
	  // get relevant Alignable
	  Alignable* ali=aap.alignableFromAlignableDet(alidet);

	  if (ali!=0) {
	    // get Alignment Parameters
	    AlignmentParameters* params = ali->alignmentParameters();
	    // get derivatives
	    AlgebraicMatrix derivs=params->selectedDerivatives(*itsos,alidet);
	    
	    // invert covariance matrix
	    int ierr; 
	    covmat.invert(ierr);
	    if (ierr != 0) { 
	      edm::LogError("Alignment") << "Matrix inversion failed!"; 
	      return; 
	    }
	    
	    // calculate user parameters
	    int npar=derivs.num_row();
	    AlgebraicSymMatrix thisjtvj(npar);
	    AlgebraicVector thisjtve(npar);
	    thisjtvj=covmat.similarity(derivs);
	    thisjtve=derivs * covmat * (pos-coor);
	    
	    // access user variables (via AlignmentParameters)
	    CSA06UserVariables* uservar =
	      dynamic_cast<CSA06UserVariables*>(params->userVariables());
	    uservar->jtvj += thisjtvj;
	    uservar->jtve += thisjtve;
	    uservar->nhit ++;

            std::pair<int,int> typeAndLay = id.typeAndLayerFromGeomDet( *det );
            TrajectoryMeasurementResidual *TMR = new TrajectoryMeasurementResidual(*imea);
            
            // fill hit parameters in root tree
	    if (ahit<MAXHIT) {
	      m_hType[ahit] = typeAndLay.first;
              m_hLayer[ahit] = typeAndLay.second;
	      m_hOwnerTrack[ahit] = allTracks;
              m_hR[ahit] = (*ihit)->globalPosition().perp();
	      m_hPhi[ahit] = (*ihit)->globalPosition().phi();
              m_hZ[ahit] = (*ihit)->globalPosition().z();
	      m_hLocalX[ahit] = (*ihit)->localPosition().x();
	      m_hLocalY[ahit] = (*ihit)->localPosition().y();
              m_hLocalZ[ahit] = (*ihit)->localPosition().z();
	      std::pair<float,float> monoStereoAng = theAngleFinder->findtrackangle(*imea);
              m_hLocalAngle[ahit] = monoStereoAng.first; 
              m_UresLoc[ahit] = coor[0] - pos[0];
              m_VresLoc[ahit] = coor[1] - pos[1];
              // m_UresLoc[ahit] = TMR->localXResidual();
              // m_VresLoc[ahit] = TMR->localYResidual();
              m_UerrLoc[ahit] = TMR->localXError();
              m_VerrLoc[ahit] = TMR->localYError();
	      m_Xres[ahit] = TMR->measurementXResidual();
              m_Yres[ahit] = TMR->measurementYResidual();
              m_Xerr[ahit] = TMR->measurementXError();
              m_Yerr[ahit] = TMR->measurementYError(); 
	      ahit++;
	      m_allHits=ahit;
	    } 
	  }
	  
	  itsos++;
	  ihit++;
          imea++;
	} 
    }
    allTracks++;
  } // end of track loop

  // fill eventwise root tree (with prescale defined in pset)
  theCurrentPrescale--;
  if (theCurrentPrescale<=0) {
    theTree->Fill();
    theCurrentPrescale=theEventPrescale;    
  }

}

// ----------------------------------------------------------------------------

int CSA06AlignmentAlgorithm::readIterationFile(string filename)
{
  int result;

  ifstream inIterFile((char*)filename.c_str(), ios::in);
  if (!inIterFile) {
    edm::LogError("Alignment") << "[CSA06AlignmentAlgorithm::readIterationFile] ERROR! "
      << "Unable to open Iteration file";
    result = -1;
  }
  else {
    inIterFile >> result;
    edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm::readIterationFile] "
         << "Read last iteration number from file: " << result;
  }
  inIterFile.close();

  return result;
}

// ----------------------------------------------------------------------------

void CSA06AlignmentAlgorithm::writeIterationFile(string filename,int iter)
{
  ofstream outIterFile((char*)(filename.c_str()), ios::out);
  if (!outIterFile) {
    edm::LogError("Alignment") << "[CSA06AlignmentAlgorithm::writeIterationFile] ERROR: Unable to write Iteration file";
  }
  else {
     outIterFile << iter;
     edm::LogWarning("Alignment") <<"[CSA06AlignmentAlgorithm::writeIterationFile] writing iteration number to file: " << iter;
  }
  outIterFile.close();
}


// ----------------------------------------------------------------------------
// set alignment position error

void CSA06AlignmentAlgorithm::setAlignmentPositionError(void)
{


  // Check if user wants to override APE
  if ( apeparam == "none" )
    {
      edm::LogWarning("Alignment") <<"[CSA06AlignmentAlgorithm::setAlignmentPositionError] No APE applied";
      return; // NO APE APPLIED
    }

  
  edm::LogWarning("Alignment") <<"[CSA06AlignmentAlgorithm::setAlignmentPositionError] Apply APE!";

  // Printout for debug
  for ( int i=0; i<21; ++i ) {
    double apelinstest=calcAPE(apesp,i,"linear");
    double apeexpstest=calcAPE(apesp,i,"exponential");
    double apelinrtest=calcAPE(aperp,i,"linear");
    double apeexprtest=calcAPE(aperp,i,"exponential");
    printf("APE: iter slin sexp rlin rexp: %5d %12.5f %12.5f %12.5f %12.5f\n",
      i,apelinstest,apeexpstest,apelinrtest,apeexprtest);
  }

  // set APE
  double apeshift=calcAPE(apesp,theIteration,apeparam);
  double aperot  =calcAPE(aperp,theIteration,apeparam);
  theAlignmentParameterStore->setAlignmentPositionError( theAlignables, apeshift, aperot );

}

// ----------------------------------------------------------------------------
// calculate APE

double 
CSA06AlignmentAlgorithm::calcAPE(double* par, int iter,std::string param)
{
  double diter=(double)iter;

  if (param == "linear") {
    return max(0.,par[0]+((par[1]-par[0])/par[2])*diter);
  }
  else if (param == "exponential") {
    return max(0.,par[0]*(exp(-pow(diter,par[1])/par[2])));
  }
  else { 
     edm::LogInfo("Alignment") << "Unknown param: " << param;
    return 0.;
  }

}


// ----------------------------------------------------------------------------
// book root trees

void CSA06AlignmentAlgorithm::bookRoot(void)
{
  // create ROOT files
  theFile = new TFile(outfile.c_str(),"update");
  theFile->cd();
  
  // book event-wise ROOT Tree

  TString tname="T1";
  char iterString[5];
  sprintf(iterString, "%i",theIteration);
  tname.Append("_");
  tname.Append(iterString);

  theTree  = new TTree(tname,"Eventwise tree");

   //theTree->Branch("Run",     &m_Run,     "Run/I");
  //theTree->Branch("Event",   &m_Event,   "Event/I");
  theTree->Branch("NtracksSim", &m_NtracksSim, "NtracksSim/I");      
  theTree->Branch("PtSim",       m_PtSim,      "PtSim[NtracksSim]/F");
  theTree->Branch("EtaSim",      m_EtaSim,     "EtaSim[NtracksSim]/F");
  theTree->Branch("PhiSim",      m_PhiSim,     "PhiSim[NtracksSim]/F");
  theTree->Branch("Ntracks", &m_Ntracks, "Ntracks/I");
  theTree->Branch("Nhits",    m_Nhits,   "Nhits[Ntracks]/I");       
  theTree->Branch("Pt",       m_Pt,      "Pt[Ntracks]/F");
  theTree->Branch("Eta",      m_Eta,     "Eta[Ntracks]/F");
  theTree->Branch("Phi",      m_Phi,     "Phi[Ntracks]/F");
  theTree->Branch("Chi2n",    m_Chi2n,   "Chi2n[Ntracks]/F");
  theTree->Branch("allHits", &m_allHits, "allHits/I");
  theTree->Branch("hType",    m_hType,   "hType[allHits]/I");
  theTree->Branch("hLayer",   m_hLayer,  "hLayer[allHits]/I");
  theTree->Branch("hOwnerTrack",m_hOwnerTrack,"hOwnerTrack[allHits]/I"); 
  theTree->Branch("hR",       m_hR,      "hR[allHits]/F");
  theTree->Branch("hPhi",     m_hPhi,    "hPhi[allHits]/F");
  theTree->Branch("hZ",       m_hZ,      "hZ[allHits]/F");
  theTree->Branch("hLocalX",  m_hLocalX, "hLocalX[allHits]/F");
  theTree->Branch("hLocalY",  m_hLocalY, "hLocalY[allHits]/F");
  theTree->Branch("hLocalZ",  m_hLocalZ, "hLocalZ[allHits]/F");
  theTree->Branch("hLocalAngle",m_hLocalAngle,"hLocalAngle[allHits]/F");
  theTree->Branch("uresLoc",  m_UresLoc, "uresLoc[allHits]/F");
  theTree->Branch("vresLoc",  m_VresLoc, "vresLoc[allHits]/F");
  theTree->Branch("uerrLoc",  m_UerrLoc, "uerrLoc[allHits]/F");
  theTree->Branch("verrLoc",  m_VerrLoc, "verrLoc[allHits]/F");
  // theTree->Branch("uresCSA06",m_UresCSA06, "uresCSA06[allHits]/F");
  // theTree->Branch("vresCSA06",m_VresCSA06, "vresCSA06[allHits]/F");
  theTree->Branch("xres",     m_Xres,    "xres[allHits]/F");
  theTree->Branch("yres",     m_Yres,    "yres[allHits]/F");
  theTree->Branch("xerr",     m_Xerr,    "xerr[allHits]/F");
  theTree->Branch("yerr",     m_Yerr,    "yerr[allHits]/F"); 

  // book Alignable-wise ROOT Tree

  theFile2 = new TFile(outfile2.c_str(),"update");
  theFile2->cd();

  theTree2 = new TTree("T2","Alignablewise tree");

  theTree2->Branch("Nhit",   &m2_Nhit,    "Nhit/I");
  theTree2->Branch("Type",   &m2_Type,    "Type/I");
  theTree2->Branch("Layer",  &m2_Layer,   "Layer/I");
  theTree2->Branch("Xpos",   &m2_Xpos,    "Xpos/F");
  theTree2->Branch("Ypos",   &m2_Ypos,    "Ypos/F");
  theTree2->Branch("Zpos",   &m2_Zpos,    "Zpos/F");
  theTree2->Branch("Eta",    &m2_Eta,     "Eta/F");
  theTree2->Branch("Phi",    &m2_Phi,     "Phi/F");
  theTree2->Branch("Id",     &m2_Id,      "Id/I");
  theTree2->Branch("ObjId",  &m2_ObjId,   "ObjId/I");

  edm::LogWarning("Alignment") << "[CSA06AlignmentAlgorithm::bookRoot] Root trees booked.";

}

// ----------------------------------------------------------------------------
// fill alignable-wise root tree

void CSA06AlignmentAlgorithm::fillRoot(void)
{
  TrackerAlignableId id;

  theFile2->cd();

  int naligned=0;

  for(vector<Alignable*>::const_iterator
    it=theAlignables.begin(); it!=theAlignables.end(); it++) {
    Alignable* ali=(*it);
    AlignmentParameters* dap = ali->alignmentParameters();

    // consider only those parameters classified as 'valid'
    if (dap->isValid()) {

      printf("------------------------------------------------------------------------\n");
      printf(" ALIGNABLE: %6d \n",naligned);

      // get number of hits from user variable
      CSA06UserVariables* uservar =
        dynamic_cast<CSA06UserVariables*>(dap->userVariables());
      m2_Nhit  = uservar->nhit;

      // get type/layer
      std::pair<int,int> tl=theAlignmentParameterStore->typeAndLayer(ali);
      m2_Type=tl.first;
      m2_Layer=tl.second;

      // get identifier (as for IO)
      m2_Id    = id.alignableId(ali);
      m2_ObjId = id.alignableTypeId(ali);

      // get position
      GlobalPoint pos=ali->surface().position();
      m2_Xpos=pos.x();
      m2_Ypos=pos.y();
      m2_Zpos=pos.z();
      m2_Eta=pos.eta();
      m2_Phi=pos.phi();

      printf("hits: %4d type: %4d layer %4d id %4d objId: %4d \n",
	     m2_Nhit,m2_Type,m2_Layer,m2_Id,m2_ObjId);
      printf("x,y,z: %12.5f %12.5f %12.5f eta,phi: %12.5f %12.5f \n",
             m2_Xpos,m2_Ypos,m2_Zpos,m2_Eta,m2_Phi);

      AlgebraicVector pars=dap->parameters();
      printf("params:  %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n",
	     pars[0],pars[1],pars[2],pars[3],pars[4],pars[5]);

      naligned++;
      theTree2->Fill();
    }
  }
}

// ----------------------------------------------------------------------------

bool CSA06AlignmentAlgorithm::calcParameters(Alignable* ali)
{

  // Alignment parameters
  AlignmentParameters* par = ali->alignmentParameters();
  // access user variables
  CSA06UserVariables* uservar =
    dynamic_cast<CSA06UserVariables*>(par->userVariables());
  int nhit = uservar->nhit;

  if (nhit < theMinimumNumberOfHits) {
    par->setValid(false);
    return false;
  }

  AlgebraicSymMatrix jtvj = uservar->jtvj;
  AlgebraicVector jtve = uservar->jtve;
  int ierr;
  AlgebraicSymMatrix jtvjinv=jtvj.inverse(ierr);
  if (ierr !=0) { 
    edm::LogError("Alignment") << "Matrix inversion failed!"; 
    return false;
  }

  // these are the alignment corrections+covariance (for selected params)
  AlgebraicVector params = - (jtvjinv * jtve);
  AlgebraicSymMatrix cov = jtvjinv;

  edm::LogInfo("Alignment") << "parameters " << params;

  // errors of parameters
  int npar=params.num_row();    
  AlgebraicVector paramerr(npar);
  AlgebraicVector relerr(npar);
  for (int i=0;i<npar;i++) {
     if (abs(cov[i][i])>0) paramerr[i]=sqrt(abs(cov[i][i]));
     else paramerr[i]=params[i];
     relerr[i] = abs(paramerr[i]/params[i]);
     if (relerr[i] >= theMaxRelParameterError) { 
       params[i]=0; 
       paramerr[i]=0; 
     }
  }

  // store alignment parameters
  AlignmentParameters* parnew = par->cloneFromSelected(params,cov);
  ali->setAlignmentParameters(parnew);
  parnew->setValid(true);
  return true;
}

//-----------------------------------------------------------------------------

void CSA06AlignmentAlgorithm::collector(void)
{
  edm::LogWarning("Alignment") <<"[CSA06AlignmentAlgorithm::collector] called for iteration " << theIteration <<endl;

  CSA06UserVariablesIORoot CSA06IO;

  for (int ijob=1; ijob<=theCollectorNJobs; ijob++) {

    edm::LogWarning("Alignment") <<"reading uservar for job " << ijob;

    stringstream ss;
    string str;
    ss << ijob;
    ss >> str;
    string uvfile = theCollectorPath+"/job"+str+"/IOUserVariables.root";

    vector<AlignmentUserVariables*> uvarvec = 
      CSA06IO.readCSA06UserVariables (theAlignables,(char*)uvfile.c_str(),
      theIteration,ioerr);

    if (ioerr!=0) { 
      edm::LogWarning("Alignment") <<"[CSA06AlignmentAlgorithm::collector] could not read user variable files!";
      return;
    }

    // add
    vector<AlignmentUserVariables*> uvarvecadd;
    vector<AlignmentUserVariables*>::const_iterator iuvarnew=uvarvec.begin(); 
    for (vector<Alignable*>::const_iterator it=theAlignables.begin(); 
     it!=theAlignables.end(); it++) {
     Alignable* ali = *it;
     AlignmentParameters* ap = ali->alignmentParameters();

     CSA06UserVariables* uvarold = 
      dynamic_cast<CSA06UserVariables*>(ap->userVariables());
     CSA06UserVariables* uvarnew = 
      dynamic_cast<CSA06UserVariables*>(*iuvarnew);

     CSA06UserVariables* uvar = uvarold->clone();
     if (uvarnew!=0) {
       uvar->nhit=(uvarold->nhit)+(uvarnew->nhit);
       uvar->jtvj=(uvarold->jtvj)+(uvarnew->jtvj);
       uvar->jtve=(uvarold->jtve)+(uvarnew->jtve);
       delete uvarnew;
     }

     uvarvecadd.push_back(uvar);
     iuvarnew++;
    }

   theAlignmentParameterStore->attachUserVariables(theAlignables,
      uvarvecadd,ioerr);

  }

}

