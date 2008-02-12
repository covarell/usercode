/// \file AlignmentProducer.cc
///
///  \author    : Frederic Ronga
///  Revision   : $Revision: 1.13 $
///  last update: $Date: 2007/08/29 02:24:34 $
///  by         : $Author: ratnik $

#include "Alignment/CommonAlignmentProducer/plugins/AlignmentProducer.h"

#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

// System include files
#include <memory>
#include <sstream>

// Framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Conditions database
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeomBuilderFromGeometricDet.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonNumberingRecord.h"
#include "Geometry/MuonNumbering/interface/MuonDDDConstants.h"
#include "Geometry/DTGeometryBuilder/src/DTGeometryBuilderFromDDD.h"
#include "Geometry/CSCGeometryBuilder/src/CSCGeometryBuilderFromDDD.h"
#include "Geometry/TrackingGeometryAligner/interface/GeometryAligner.h"
#include "CondFormats/DataRecord/interface/TrackerAlignmentRcd.h"
#include "CondFormats/DataRecord/interface/TrackerAlignmentErrorRcd.h"
#include "CondFormats/DataRecord/interface/DTAlignmentRcd.h"
#include "CondFormats/DataRecord/interface/DTAlignmentErrorRcd.h"
#include "CondFormats/DataRecord/interface/CSCAlignmentRcd.h"
#include "CondFormats/DataRecord/interface/CSCAlignmentErrorRcd.h"
#include "CondFormats/DataRecord/interface/TrackerSurveyRcd.h"
#include "CondFormats/DataRecord/interface/TrackerSurveyErrorRcd.h"
#include "CondFormats/DataRecord/interface/DTSurveyRcd.h"
#include "CondFormats/DataRecord/interface/DTSurveyErrorRcd.h"
#include "CondFormats/DataRecord/interface/CSCSurveyRcd.h"
#include "CondFormats/DataRecord/interface/CSCSurveyErrorRcd.h"

// Tracking 	 
#include "TrackingTools/PatternTools/interface/Trajectory.h" 

// Alignment
#include "CondFormats/Alignment/interface/Alignments.h"
#include "CondFormats/Alignment/interface/AlignmentErrors.h"
#include "CondFormats/Alignment/interface/SurveyErrors.h"
#include "Alignment/TrackerAlignment/interface/TrackerScenarioBuilder.h"
#include "Alignment/MuonAlignment/interface/MuonScenarioBuilder.h"
#include "Alignment/CommonAlignment/interface/SurveyDet.h"
#include "Alignment/CommonAlignmentParametrization/interface/RigidBodyAlignmentParameters.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmPluginFactory.h"
#include "Alignment/CommonAlignmentMonitor/interface/AlignmentMonitorPluginFactory.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterSelector.h"

//_____________________________________________________________________________
AlignmentProducer::AlignmentProducer(const edm::ParameterSet& iConfig) :
  theAlignmentAlgo(0), theAlignmentParameterStore(0),
  theAlignableTracker(0), theAlignableMuon(0),
  nevent_(0), theParameterSet(iConfig),
  theMaxLoops( iConfig.getUntrackedParameter<unsigned int>("maxLoops",0) ),
  stNFixAlignables_(iConfig.getParameter<int>("nFixAlignables") ),
  stRandomShift_(iConfig.getParameter<double>("randomShift")),
  stRandomRotation_(iConfig.getParameter<double>("randomRotation")),
  applyDbAlignment_( iConfig.getUntrackedParameter<bool>("applyDbAlignment",false) ),
  doMisalignmentScenario_(iConfig.getParameter<bool>("doMisalignmentScenario")),
  saveToDB_(iConfig.getParameter<bool>("saveToDB")),
  doTracker_( iConfig.getUntrackedParameter<bool>("doTracker") ),
  doMuon_( iConfig.getUntrackedParameter<bool>("doMuon") ),
  useSurvey_( iConfig.getParameter<bool>("useSurvey") ),
  isData_( iConfig.getUntrackedParameter<bool>("isData") )
{

  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::AlignmentProducer";

  // Tell the framework what data is being produced
  if (doTracker_) {
     setWhatProduced(this, &AlignmentProducer::produceTracker);
  }
  if (doMuon_) {
     setWhatProduced(this, &AlignmentProducer::produceDT);
     setWhatProduced(this, &AlignmentProducer::produceCSC);
  }

  // Create the alignment algorithm
  edm::ParameterSet algoConfig = iConfig.getParameter<edm::ParameterSet>( "algoConfig" );
  std::string algoName = algoConfig.getParameter<std::string>("algoName");
  theAlignmentAlgo = AlignmentAlgorithmPluginFactory::get( )->create( algoName, algoConfig  );

  // Check if found
  if ( !theAlignmentAlgo )
	throw cms::Exception("BadConfig") << "Couldn't find algorithm called " << algoName;

  edm::ParameterSet monitorConfig = iConfig.getParameter<edm::ParameterSet>( "monitorConfig" );
  std::vector<std::string> monitors = monitorConfig.getUntrackedParameter<std::vector<std::string> >( "monitors" );

  for (std::vector<std::string>::const_iterator miter = monitors.begin();  miter != monitors.end();  ++miter) {
//     AlignmentMonitorBase* newMonitor = dynamic_cast<AlignmentMonitorBase*>(
//	AlignmentMonitorPluginFactory::getMonitor(*miter, monitorConfig.getParameter<edm::ParameterSet>(*miter)));
     AlignmentMonitorBase* newMonitor = AlignmentMonitorPluginFactory::get()->create(*miter, monitorConfig.getUntrackedParameter<edm::ParameterSet>(*miter));

     if (!newMonitor) throw cms::Exception("BadConfig") << "Couldn't find monitor named " << *miter;

     theMonitors.push_back(newMonitor);
  }

}


//_____________________________________________________________________________
// Delete new objects
AlignmentProducer::~AlignmentProducer()
{

  delete theAlignmentParameterStore;
  delete theAlignableTracker;
  delete theAlignableMuon;

}


//_____________________________________________________________________________
// Produce tracker geometry
boost::shared_ptr<TrackerGeometry>
AlignmentProducer::produceTracker( const TrackerDigiGeometryRecord& iRecord )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::produceTracker";
  return theTracker;
}

//_____________________________________________________________________________
// Produce muonDT geometry
boost::shared_ptr<DTGeometry>
AlignmentProducer::produceDT( const MuonGeometryRecord& iRecord )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::produceDT";
  return theMuonDT;
}

//_____________________________________________________________________________
// Produce muonCSC geometry
boost::shared_ptr<CSCGeometry>
AlignmentProducer::produceCSC( const MuonGeometryRecord& iRecord )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::produceCSC";
  return theMuonCSC;  
}


//_____________________________________________________________________________
// Initialize algorithm
void AlignmentProducer::beginOfJob( const edm::EventSetup& iSetup )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob";

  nevent_ = 0;
  GeometryAligner aligner;

  // Create the geometries from the ideal geometries (first time only)
  this->createGeometries_( iSetup );
  
  // Retrieve and apply alignments, if requested (requires DB setup)
  if ( applyDbAlignment_ ) {
    if ( doTracker_ ) {
      edm::ESHandle<Alignments> alignments;
      iSetup.get<TrackerAlignmentRcd>().get( alignments );
      edm::ESHandle<AlignmentErrors> alignmentErrors;
      iSetup.get<TrackerAlignmentErrorRcd>().get( alignmentErrors );
      aligner.applyAlignments<TrackerGeometry>( &(*theTracker), &(*alignments), &(*alignmentErrors) );
    }
    if ( doMuon_ ) {
      edm::ESHandle<Alignments> dtAlignments;
      iSetup.get<DTAlignmentRcd>().get( dtAlignments );
      edm::ESHandle<AlignmentErrors> dtAlignmentErrors;
      iSetup.get<DTAlignmentErrorRcd>().get( dtAlignmentErrors );
      aligner.applyAlignments<DTGeometry>( &(*theMuonDT), &(*dtAlignments), &(*dtAlignmentErrors) );

      edm::ESHandle<Alignments> cscAlignments;
      iSetup.get<CSCAlignmentRcd>().get( cscAlignments );
      edm::ESHandle<AlignmentErrors> cscAlignmentErrors;
      iSetup.get<CSCAlignmentErrorRcd>().get( cscAlignmentErrors );
      aligner.applyAlignments<CSCGeometry>( &(*theMuonCSC), &(*cscAlignments), &(*cscAlignmentErrors) );
    }
  }

  // Create alignable tracker and muon 
  if (doTracker_)
  {
    theAlignableTracker = new AlignableTracker( &(*theGeometricDet), &(*theTracker) );

    if (useSurvey_)
    {
      edm::ESHandle<Alignments> surveys;
      edm::ESHandle<SurveyErrors> surveyErrors;

      iSetup.get<TrackerSurveyRcd>().get(surveys);
      iSetup.get<TrackerSurveyErrorRcd>().get(surveyErrors);

      theSurveyIndex  = 0;
      theSurveyValues = &*surveys;
      theSurveyErrors = &*surveyErrors;
      addSurveyInfo_(theAlignableTracker);
    }
  }

  if (doMuon_)
  {
     theAlignableMuon = new AlignableMuon( &(*theMuonDT), &(*theMuonCSC) );

     if (useSurvey_)
     {
	edm::ESHandle<Alignments> dtSurveys;
	edm::ESHandle<SurveyErrors> dtSurveyErrors;
	edm::ESHandle<Alignments> cscSurveys;
	edm::ESHandle<SurveyErrors> cscSurveyErrors;
	
	iSetup.get<DTSurveyRcd>().get(dtSurveys);
	iSetup.get<DTSurveyErrorRcd>().get(dtSurveyErrors);
	iSetup.get<CSCSurveyRcd>().get(cscSurveys);
	iSetup.get<CSCSurveyErrorRcd>().get(cscSurveyErrors);

	theSurveyIndex  = 0;
	theSurveyValues = &*dtSurveys;
	theSurveyErrors = &*dtSurveyErrors;
	std::vector<Alignable*> barrels = theAlignableMuon->DTBarrel();
	for (std::vector<Alignable*>::const_iterator iter = barrels.begin();  iter != barrels.end();  ++iter) {
	   addSurveyInfo_(*iter);
	}

	theSurveyIndex  = 0;
	theSurveyValues = &*cscSurveys;
	theSurveyErrors = &*cscSurveyErrors;
	std::vector<Alignable*> endcaps = theAlignableMuon->CSCEndcaps();
	for (std::vector<Alignable*>::const_iterator iter = endcaps.begin();  iter != endcaps.end();  ++iter) {
	   addSurveyInfo_(*iter);
	}
     }
  }

  // Create alignment parameter builder
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                            << "Creating AlignmentParameterBuilder";
  edm::ParameterSet aliParamBuildCfg = 
    theParameterSet.getParameter<edm::ParameterSet>("ParameterBuilder");
  AlignmentParameterBuilder alignmentParameterBuilder(theAlignableTracker,
                                                      theAlignableMuon,
                                                      aliParamBuildCfg );
  // Fix alignables if requested
  if (stNFixAlignables_>0) alignmentParameterBuilder.fixAlignables(stNFixAlignables_);

  // Get list of alignables
  Alignables theAlignables = alignmentParameterBuilder.alignables();
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                            << "got " << theAlignables.size() << " alignables";

  // Create AlignmentParameterStore 
  edm::ParameterSet aliParamStoreCfg = 
    theParameterSet.getParameter<edm::ParameterSet>("ParameterStore");
  theAlignmentParameterStore = new AlignmentParameterStore(theAlignables, aliParamStoreCfg);
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                            << "AlignmentParameterStore created!";

  // Apply misalignment scenario to alignable tracker and muon if requested
  // WARNING: this assumes scenarioConfig can be passed to both muon and tracker
  if (doMisalignmentScenario_ && (doTracker_ || doMuon_)) {
    edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                              << "Applying misalignment scenario to "
                              << (doTracker_ ? "tracker" : "")
                              << (doMuon_    ? (doTracker_ ? " and muon" : "muon") : ".");
    edm::ParameterSet scenarioConfig 
      = theParameterSet.getParameter<edm::ParameterSet>( "MisalignmentScenario" );
    if (doTracker_) {
      TrackerScenarioBuilder scenarioBuilder( theAlignableTracker );
      scenarioBuilder.applyScenario( scenarioConfig );
    }
    if (doMuon_) {
      MuonScenarioBuilder muonScenarioBuilder( theAlignableMuon );
      muonScenarioBuilder.applyScenario( scenarioConfig );
    }
  } else {
    edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::beginOfJob" 
                              << "NOT applying misalignment scenario!";
  }

  // Apply simple misalignment
  const std::string sParSel(theParameterSet.getParameter<std::string>("parameterSelectorSimple"));
  this->simpleMisalignment_(theAlignables, sParSel, stRandomShift_, stRandomRotation_, true);

  // Initialize alignment algorithm
  theAlignmentAlgo->initialize( iSetup, theAlignableTracker,
                                theAlignableMuon, theAlignmentParameterStore );

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();  monitor != theMonitors.end();  ++monitor) {
     (*monitor)->beginOfJob(theAlignableTracker, theAlignableMuon, theAlignmentParameterStore);
  }
}

//_____________________________________________________________________________
// Terminate algorithm
void AlignmentProducer::endOfJob()
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::endOfJob";

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();  monitor != theMonitors.end();  ++monitor) {
     (*monitor)->endOfJob();
  }

  // Save alignments to database
  if (saveToDB_) {
    edm::LogInfo("Alignment") << "Writing Alignments to DB...";
    // Call service
    edm::Service<cond::service::PoolDBOutputService> poolDbService;
    if( !poolDbService.isAvailable() ) // Die if not available
      throw cms::Exception("NotAvailable") << "PoolDBOutputService not available";
    
    if ( doTracker_ ) {
       // Get alignments+errors
       Alignments* alignments = theAlignableTracker->alignments();
       AlignmentErrors* alignmentErrors = theAlignableTracker->alignmentErrors();

       // Store
       std::string alignRecordName( "TrackerAlignmentRcd" );
       std::string errorRecordName( "TrackerAlignmentErrorRcd" );

       if ( poolDbService->isNewTagRequest(alignRecordName) )
	  poolDbService->createNewIOV<Alignments>( alignments, poolDbService->endOfTime(), 
						   alignRecordName );
       else
	  poolDbService->appendSinceTime<Alignments>( alignments, poolDbService->currentTime(), 
						      alignRecordName );
       if ( poolDbService->isNewTagRequest(errorRecordName) )
	  poolDbService->createNewIOV<AlignmentErrors>( alignmentErrors, poolDbService->endOfTime(), 
							errorRecordName );
       else
	  poolDbService->appendSinceTime<AlignmentErrors>( alignmentErrors, poolDbService->currentTime(), 
							   errorRecordName );
    }
 
    if ( doMuon_ ) {
       // Get alignments+errors
       Alignments*      dtAlignments       = theAlignableMuon->dtAlignments();
       AlignmentErrors* dtAlignmentErrors  = theAlignableMuon->dtAlignmentErrors();
       Alignments*      cscAlignments      = theAlignableMuon->cscAlignments();
       AlignmentErrors* cscAlignmentErrors = theAlignableMuon->cscAlignmentErrors();

       std::string dtAlignRecordName( "DTAlignmentRcd" );
       std::string dtErrorRecordName( "DTAlignmentErrorRcd" );
       std::string cscAlignRecordName( "CSCAlignmentRcd" );
       std::string cscErrorRecordName( "CSCAlignmentErrorRcd" );

       if (poolDbService->isNewTagRequest(dtAlignRecordName)) {
	  poolDbService->createNewIOV<Alignments>( &(*dtAlignments), poolDbService->endOfTime(), dtAlignRecordName);
       }
       else {
	  poolDbService->appendSinceTime<Alignments>( &(*dtAlignments), poolDbService->currentTime(), dtAlignRecordName);
       }
       if (poolDbService->isNewTagRequest(dtErrorRecordName)) {
	  poolDbService->createNewIOV<AlignmentErrors>( &(*dtAlignmentErrors), poolDbService->endOfTime(), dtErrorRecordName);
       }
       else {
	  poolDbService->appendSinceTime<AlignmentErrors>( &(*dtAlignmentErrors), poolDbService->currentTime(), dtErrorRecordName);
       }
       if (poolDbService->isNewTagRequest(cscAlignRecordName)) {
	  poolDbService->createNewIOV<Alignments>( &(*cscAlignments), poolDbService->endOfTime(), cscAlignRecordName);
       }
       else {
	  poolDbService->appendSinceTime<Alignments>( &(*cscAlignments), poolDbService->currentTime(), cscAlignRecordName);
       }
       if (poolDbService->isNewTagRequest(cscErrorRecordName)) {
	  poolDbService->createNewIOV<AlignmentErrors>( &(*cscAlignmentErrors), poolDbService->endOfTime(), cscErrorRecordName);
       }
       else {
	  poolDbService->appendSinceTime<AlignmentErrors>( &(*cscAlignmentErrors), poolDbService->currentTime(), cscErrorRecordName);
       }
    }
  }
}

//_____________________________________________________________________________
// Called at beginning of loop
void AlignmentProducer::startingNewLoop(unsigned int iLoop )
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::startingNewLoop" 
                            << "Starting loop number " << iLoop;

  theAlignmentAlgo->startNewLoop();

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();  monitor != theMonitors.end();  ++monitor) {
     (*monitor)->startingNewLoop();
  }

  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::startingNewLoop" 
                            << "Now physically apply alignments to  geometry...";


  // Propagate changes to reconstruction geometry (from initialisation or iteration)
  GeometryAligner aligner;
  if ( doTracker_ ) {
    std::auto_ptr<Alignments> alignments(theAlignableTracker->alignments());
    std::auto_ptr<AlignmentErrors> alignmentErrors(theAlignableTracker->alignmentErrors());
    aligner.applyAlignments<TrackerGeometry>( &(*theTracker),&(*alignments),&(*alignmentErrors));
  }
  if ( doMuon_ ) {
    std::auto_ptr<Alignments>      dtAlignments(       theAlignableMuon->dtAlignments());
    std::auto_ptr<AlignmentErrors> dtAlignmentErrors(  theAlignableMuon->dtAlignmentErrors());
    std::auto_ptr<Alignments>      cscAlignments(      theAlignableMuon->cscAlignments());
    std::auto_ptr<AlignmentErrors> cscAlignmentErrors( theAlignableMuon->cscAlignmentErrors());

    aligner.applyAlignments<DTGeometry>( &(*theMuonDT), &(*dtAlignments), &(*dtAlignmentErrors) );
    aligner.applyAlignments<CSCGeometry>( &(*theMuonCSC), &(*cscAlignments), &(*cscAlignmentErrors) );
  }
}


//_____________________________________________________________________________
// Called at end of loop
edm::EDLooper::Status 
AlignmentProducer::endOfLoop(const edm::EventSetup& iSetup, unsigned int iLoop)
{
  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::endOfLoop" 
                            << "Ending loop " << iLoop;

  edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::endOfLoop" 
                            << "Terminating algorithm.";
  theAlignmentAlgo->terminate();

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();  monitor != theMonitors.end();  ++monitor) {
     (*monitor)->endOfLoop(iSetup);
  }

  if ( iLoop == theMaxLoops-1 || iLoop >= theMaxLoops ) return kStop;
  else return kContinue;
}

//_____________________________________________________________________________
// Called at each event
edm::EDLooper::Status 
AlignmentProducer::duringLoop( const edm::Event& event, 
  const edm::EventSetup& setup )
{
  nevent_++;

  for ( int i=10; i<10000000; i*=10 )
    if ( nevent_<10*i && (nevent_%i)==0 )
      edm::LogInfo("Alignment") << "@SUB=AlignmentProducer::duringLoop" 
                                << "Events processed: " << nevent_;

  // Retrieve trajectories and tracks from the event
  edm::InputTag tjTag = theParameterSet.getParameter<edm::InputTag>("tjTkAssociationMapTag");
  edm::Handle<TrajTrackAssociationCollection> m_TrajTracksMap;
  event.getByLabel( tjTag, m_TrajTracksMap );
  edm::InputTag simTag = theParameterSet.getParameter<edm::InputTag>("simtrackTag");
  edm::Handle<edm::SimTrackContainer> m_SimCollection;
  if (!isData_) event.getByLabel( simTag, m_SimCollection );

  // Form pairs of trajectories and tracks
  ConstTrajTrackPairCollection trajTracks;
  for ( TrajTrackAssociationCollection::const_iterator iPair = m_TrajTracksMap->begin();
        iPair != m_TrajTracksMap->end(); iPair++ )
    trajTracks.push_back( ConstTrajTrackPair( &(*(*iPair).key), &(*(*iPair).val) ) );

  // Run the alignment algorithm
  theAlignmentAlgo->run(  setup, trajTracks, *(m_SimCollection.product()) );

  for (std::vector<AlignmentMonitorBase*>::const_iterator monitor = theMonitors.begin();  monitor != theMonitors.end();  ++monitor) {
     (*monitor)->duringLoop(setup, trajTracks);
  }

  return kContinue;
}

// ----------------------------------------------------------------------------

void AlignmentProducer::simpleMisalignment_(const Alignables &alivec, const std::string &selection, 
                                            float shift, float rot, bool local)
{

  std::ostringstream output; // collecting output

  if (shift > 0. || rot > 0.) {
    output << "Adding random flat shift of max size " << shift
           << " and adding random flat rotation of max size " << rot <<" to ";

    std::vector<bool> commSel(0);
    if (selection != "-1") {
      AlignmentParameterSelector aSelector(0,0); // no alignable needed here...
      const std::vector<char> cSel(aSelector.convertParamSel(selection));
      if (cSel.size() < RigidBodyAlignmentParameters::N_PARAM) {
	throw cms::Exception("BadConfig") 
	  << "[AlignmentProducer::simpleMisalignment_]\n"
	  << "Expect selection string '" << selection << "' to be at least of length " 
	  << RigidBodyAlignmentParameters::N_PARAM << " or to be '-1'.\n"
	  << "(Most probably you have to adjust the parameter 'parameterSelectorSimple'.)";
      }
      for (std::vector<char>::const_iterator cIter = cSel.begin(); cIter != cSel.end(); ++cIter) {
        commSel.push_back(*cIter == '0' ? false : true);
      }
      output << "parameters defined by (" << selection 
             << "), representing (x,y,z,alpha,beta,gamma),";
    } else {
      output << "the active parameters of each alignable,";
    }
    output << " in " << (local ? "local" : "global") << " frame.";

    for (std::vector<Alignable*>::const_iterator it = alivec.begin(); it != alivec.end(); ++it) {
      Alignable* ali=(*it);
      std::vector<bool> mysel(commSel.empty() ? ali->alignmentParameters()->selector() : commSel);
      
      if (std::abs(shift)>0.00001) {
	double s0 = 0., s1 = 0., s2 = 0.;
        if (mysel[RigidBodyAlignmentParameters::dx]) s0 = shift * double(random()%1000-500)/500.;
        if (mysel[RigidBodyAlignmentParameters::dy]) s1 = shift * double(random()%1000-500)/500.;
        if (mysel[RigidBodyAlignmentParameters::dz]) s2 = shift * double(random()%1000-500)/500.;
        
        if (local) ali->move( ali->surface().toGlobal(align::LocalVector(s0,s1,s2)) );
	else       ali->move( align::GlobalVector(s0,s1,s2) );

      //AlignmentPositionError ape(dx,dy,dz);
      //ali->addAlignmentPositionError(ape);
      }

      if (std::abs(rot)>0.00001) {
	align::EulerAngles r(3);
        if (mysel[RigidBodyAlignmentParameters::dalpha]) r(1)=rot*double(random()%1000-500)/500.;
        if (mysel[RigidBodyAlignmentParameters::dbeta])  r(2)=rot*double(random()%1000-500)/500.;
        if (mysel[RigidBodyAlignmentParameters::dgamma]) r(3)=rot*double(random()%1000-500)/500.;

        const align::RotationType mrot = align::toMatrix(r);
        if (local) ali->rotateInLocalFrame(mrot);
        else       ali->rotateInGlobalFrame(mrot);
        
      //ali->addAlignmentPositionErrorFromRotation(mrot);
      }
    } // end loop on alignables
  } else {
    output << "No simple misalignment added!";
  }
  edm::LogInfo("Alignment")  << "@SUB=AlignmentProducer::simpleMisalignment_" << output.str();
}


//__________________________________________________________________________________________________
void AlignmentProducer::createGeometries_( const edm::EventSetup& iSetup )
{
   edm::ESHandle<DDCompactView> cpv;
   iSetup.get<IdealGeometryRecord>().get( cpv );

   if (doTracker_) {
      iSetup.get<IdealGeometryRecord>().get( theGeometricDet );
      TrackerGeomBuilderFromGeometricDet trackerBuilder;
      theTracker = boost::shared_ptr<TrackerGeometry>( trackerBuilder.build(&(*theGeometricDet)) );
   }

   if (doMuon_) {
      edm::ESHandle<MuonDDDConstants> mdc;
      iSetup.get<MuonNumberingRecord>().get(mdc);
      DTGeometryBuilderFromDDD DTGeometryBuilder;
      CSCGeometryBuilderFromDDD CSCGeometryBuilder;
      theMuonDT = boost::shared_ptr<DTGeometry>(DTGeometryBuilder.build(&(*cpv), *mdc));

      //theMuonCSC = boost::shared_ptr<CSCGeometry>(CSCGeometryBuilder.build(&(*cpv), *mdc));
      theMuonCSC = boost::shared_ptr<CSCGeometry>( new CSCGeometry );
      CSCGeometryBuilder.build( theMuonCSC, &(*cpv), *mdc );
   }
}

void AlignmentProducer::addSurveyInfo_(Alignable* ali)
{
  const std::vector<Alignable*>& comp = ali->components();

  unsigned int nComp = comp.size();

  for (unsigned int i = 0; i < nComp; ++i) addSurveyInfo_(comp[i]);

  const SurveyError& error = theSurveyErrors->m_surveyErrors[theSurveyIndex];

  if ( ali->geomDetId().rawId() != error.rawId() ||
       ali->alignableObjectId() != error.structureType() )
  {
    throw cms::Exception("DatabaseError")
      << "Error reading survey info from DB. Mismatched id!";
  }

  const CLHEP::Hep3Vector&  pos = theSurveyValues->m_align[theSurveyIndex].translation();
  const CLHEP::HepRotation& rot = theSurveyValues->m_align[theSurveyIndex].rotation();

  AlignableSurface surf( align::PositionType( pos.x(), pos.y(), pos.z() ),
			 align::RotationType( rot.xx(), rot.xy(), rot.xz(),
					      rot.yx(), rot.yy(), rot.yz(),
					      rot.zx(), rot.zy(), rot.zz() ) );

  surf.setWidth( ali->surface().width() );
  surf.setLength( ali->surface().length() );

  ali->setSurvey( new SurveyDet( surf, error.matrix() ) );

  ++theSurveyIndex;
}

DEFINE_ANOTHER_FWK_LOOPER( AlignmentProducer );
