// -*- C++ -*-
//
// Package:    EvtGenInterface
// Class:      EvtGenProducer
// 
/**\class EvtGenProducer GeneratorInterface/EvtGenInterface/src/EvtGenProducer.cc


 Description: EvtGen interface - decays B mesons (left stable by Pythia) by EvtGen

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Nello Nappi
//         Created:  Fri May 11 15:19:32 CEST 2007
// $Id: EvtGenProducer.cc,v 1.2 2007/07/03 13:29:15 covarell Exp $
//
//
#include "FWCore/PluginManager/interface/PluginManager.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "HepMC/IO_HEPEVT.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "Utilities/General/interface/FileInPath.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/Random.h"

#include "GeneratorInterface/EvtGenInterface/interface/EvtGenProducer.h"
#include "GeneratorInterface/EvtGenInterface/interface/myEvtRandomEngine.h"


#include <iostream>

EvtGenProducer::EvtGenProducer(edm::ParameterSet const & p)
{   

// create random engine and initialize seed using Random Number Generator Service 
// as configured in the configuration file

  edm::Service<edm::RandomNumberGenerator> rngen;
  
  if ( ! rngen.isAvailable()) {                //check available, else generate exception
    throw cms::Exception("Configuration")
      << "The EvtGenProducer module requires the RandomNumberGeneratorService\n"
      "which is not present in the configuration file.  You must add the service\n"
      "in the configuration file if you  want to run EvtGenProducer";
    }
  CLHEP::HepRandomEngine& m_engine = rngen->getEngine();
  m_flat = new CLHEP::RandFlat(m_engine, 0., 1.);
  //  HepRandom::setTheEngine(m_engine);   not needed ????
  myEvtRandomEngine* the_engine = new myEvtRandomEngine(&m_engine); 

  // Get data from parameter set

  std::string decay_table = p.getParameter<std::string>("decay_table");
  std::string pdt = p.getParameter<std::string>("particle_property_file");
  // any number of alias names for forced decays can be specified using dynamic std vector of strings
  std::vector<std::string> forced_names = p.getParameter< std::vector<std::string> >("list_forced_decays");
    
  produces<edm::HepMCProduct>();   // declare 
    
  m_EvtGen = new EvtGen (decay_table.c_str(),pdt.c_str(),the_engine);  // 4 th parameter should be rad cor
  // EvtPythia::pythiaInit(0);      // Patrick Robbe's advice

  std::vector<std::string>::const_iterator i;
  nforced=0;
  for (i=forced_names.begin(); i!=forced_names.end(); ++i)   // i will point to strings containing
                                                            // names of particles with forced decay
    {
      nforced++;
      EvtId found = EvtPDL::getId(*i);        // EvtPDL::getID finds the EvtId corresponding to name
      if (found.getId()==-1)
	{
	  throw cms::Exception("Configuration")
	    << "name in part list for forced decays not found: " << *i; // ??? OK ????
	}
      if (found.getId()==found.getAlias())
	{
	  throw cms::Exception("Configuration")
	    << "name in part list for forced decays is not an alias: " << *i; // ??? OK ????
	}
      forced_Evt.push_back(found);                   // forced_Evt is the list of EvtId's
      forced_Hep.push_back(EvtPDL::getStdHep(found));// forced_Hep is the list of stdhep codes
    }

}

EvtGenProducer::~EvtGenProducer() 
{ 
  //this is causing a seg fault when an exception occurs while constructing
  // an HcalSD.  Need to check for memory problems. 
  // if (m_runManager!=0) delete m_runManager; 
}

void EvtGenProducer::beginJob(const edm::EventSetup & es)
{
  ntotal = 0;
  std::cout << " EvtGenProducer starting ... " << std::endl;
  /*    StaticRandomEngineSetUnset random(m_engine);

    std::cout << " EvtGenProducer initializing " << std::endl;
    m_runManager->initG4(es);
  */
}
 
void EvtGenProducer::endJob()
{ std::cout << " EvtGenProducer terminating ... " << std::endl; }
 
void EvtGenProducer::produce(edm::Event & e, const edm::EventSetup & es)
{
  int idHep,ipart,status;
  EvtId idEvt;
  std::auto_ptr<edm::HepMCProduct> new_product( new edm::HepMCProduct() ); 

  edm::Handle< edm::HepMCProduct > EvtHandle ;
  e.getByLabel( "source", EvtHandle ) ;
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
  /*
  StaticRandomEngineSetUnset random(m_engine);  
  */

  // Loop through existing particles to find undecayed B's
  int nlist = 0; 
  HepMC::GenParticle *listp[10]; 
  int index[10];       // list of candidates to be forced 
  for (HepMC::GenEvent::particle_const_iterator p=Evt->particles_begin(); p != Evt->particles_end(); ++p)
    {
      status = (*p)->status();
 
      if(status!=2)                 // only not decayed (status = 2) particles
	{ 
	  idHep = (*p)->pdg_id();

	  int do_force=0;
	  for(int i=0;i<nforced;i++)           // First check if part with forced decay
	    {                                  // In that case do not decay immediately 
	      if(idHep==forced_Hep[i])         // (only 1 per event will be forced)	 
		{                              // Fill list
		  if(nlist<10)                 // not nice ... but 10 is a reasonable maximum?
		    {
		      listp[nlist]=*p;
		      index[nlist++]=i;
		    }
		  else
		    {
		      throw cms::Exception("runtime")
			<< "more than 10 candidates to be forced "; 
		    }
		  do_force=1;
		}
	    }
	  if(do_force==0)         // particles with decays not forced are decayed immediately 
	    {
	      idEvt = EvtPDL::evtIdFromStdHep(idHep);
	      ipart = idEvt.getId();
	      if(ipart==-1)continue;                          // particle not known to EvtGen       
	      if(EvtDecayTable::getNMode(ipart)==0)continue;  // particles stable for EvtGen
	      decay(*p,idEvt);                                 // generate decay
	    }
	}
    }
  if(nlist!=0)
    {
      // decide randomly which one to decay as alias
      int which = (int)(nlist*m_flat->fire()); 
      if(which==nlist)which=nlist-1;
	  for(int k=0;k<nlist;k++)
	    {
	      if(k==which)
		{
		  decay(listp[k],forced_Evt[k]);           // decay as alias
		}
	      else
		{
		  int id_non_alias=forced_Evt[k].getId();
		  EvtId non_alias(id_non_alias,id_non_alias); // create new EvtId with id = alias
		  decay(listp[k],non_alias);                    // decay as standard (non alias
		}
	    }
    }

    HepMC::IO_HEPEVT conv;
    //HepMC::GenEvent* evt = conv.getGenEventfromHEPEVT();
    HepMC::GenEvent* evt = conv.read_next_event();
    evt->set_signal_process_id( Evt->signal_process_id() );
    evt->set_event_scale( Evt->event_scale() );
    evt->set_event_number( Evt->event_number() );

    if (evt) new_product->addHepMCData( evt );

    e.put( new_product );

}


void EvtGenProducer::decay(HepMC::GenParticle* partHep, EvtId idEvt)
{
  // Set spin type
  EvtSpinType::spintype stype = EvtPDL::getSpinType(idEvt);
  EvtParticle* partEvt;
    switch (stype){
    case EvtSpinType::SCALAR: 
      partEvt = new EvtScalarParticle();
      break;
    case EvtSpinType::STRING:
      partEvt = new EvtStringParticle();    // ?????
      break;
    case EvtSpinType::DIRAC: 
      partEvt = new EvtDiracParticle();
      break;
    case EvtSpinType::VECTOR:
      partEvt = new EvtVectorParticle();
      break;
    case EvtSpinType::RARITASCHWINGER:
      partEvt = new EvtRaritaSchwingerParticle();
      break;
    case EvtSpinType::TENSOR:
      partEvt = new EvtTensorParticle();
      break;
    case EvtSpinType::SPIN5HALF: case EvtSpinType::SPIN3: case EvtSpinType::SPIN7HALF: case EvtSpinType::SPIN4:
      partEvt = new EvtHighSpinParticle();
      break;
    default:
      std::cout << "Unknown spintype in EvtSpinType!" << std::endl;   
      return;
    }

    // Generate decay
    EvtVector4R momEvt; // translate particle 4 momentum from Hep to Evt format
    HepMC::FourVector momHep = partHep->momentum();
    momEvt.set(momHep.t(),momHep.x(),momHep.y(),momHep.z());
    partEvt->init(idEvt,momEvt);
    partEvt->setDiagonalSpinDensity();        // unpolarized ??? 
    partEvt->decay();                    

    if (ntotal % 2000 == 0) {
      partEvt->printParticle();                    // DEBUG
      partEvt->printTree();
      std::cout << "--------------------------------------------" << std::endl;
    }
    ntotal++;
    partEvt->deleteTree();

    // Store particle in HEPEVT format
    static EvtStdHep evtstdhep;
    static EvtSecondary evtsecondary;
    EvtId        list_of_stable[10];
    EvtParticle* stable_parent[10];
        
    list_of_stable[0]=EvtId(-1,-1);
    stable_parent[0]=0;
    
    evtstdhep.init();
    evtsecondary.init();
    partEvt->makeStdHep(evtstdhep,evtsecondary,list_of_stable);
}        

/* StaticRandomEngineSetUnset::StaticRandomEngineSetUnset() {

    using namespace edm;
    edm::Service<RandomNumberGenerator> rng;

    if ( ! rng.isAvailable()) {
        throw cms::Exception("Configuration")
            << "The EvtGenProducer module requires the RandomNumberGeneratorService\n"
               "which is not present in the configuration file.  You must add the service\n"
               "in the configuration file if you want to run EvtGenProducer";
    }
    m_currentEngine = &(rng->getEngine());

    m_previousEngine = HepRandom::getTheEngine();
    HepRandom::setTheEngine(m_currentEngine);
}

StaticRandomEngineSetUnset::StaticRandomEngineSetUnset(CLHEP::HepRandomEngine * engine) {

    m_currentEngine = engine;

    m_previousEngine = HepRandom::getTheEngine();
    HepRandom::setTheEngine(m_currentEngine);
}

StaticRandomEngineSetUnset::~StaticRandomEngineSetUnset() {
    HepRandom::setTheEngine(m_previousEngine);
}

CLHEP::HepRandomEngine*
StaticRandomEngineSetUnset::getEngine() const { return m_currentEngine; } */

// DEFINE_FWK_MODULE(EvtGenProducer);
 
