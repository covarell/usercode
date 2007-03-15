/**
 * \file AlignmentParameterStore.cc
 *
 *  $Revision: 1.9 $
 *  $Date: 2007/02/12 16:06:19 $
 *  (last update by $Author: flucke $)
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "Alignment/CommonAlignment/interface/AlignableDet.h"
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"

#include "Alignment/CommonAlignmentParametrization/interface/AlignmentTransformations.h"
#include "Alignment/CommonAlignmentParametrization/interface/RigidBodyAlignmentParameters.h"

// This class's header
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentParameterStore.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentCorrelationsStore.h"
#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentExtendedCorrelationsStore.h"


//__________________________________________________________________________________________________
AlignmentParameterStore::AlignmentParameterStore( const Alignables &alis,
						  const edm::ParameterSet& config ) :
  theAlignables(alis)
{
  if (config.getUntrackedParameter<bool>("UseExtendedCorrelations")) {
    theCorrelationsStore = new AlignmentExtendedCorrelationsStore
      (config.getParameter<edm::ParameterSet>("ExtendedCorrelationsConfig"));
  } else {
    theCorrelationsStore = new AlignmentCorrelationsStore();
  }

  theTrackerAlignableId = new TrackerAlignableId;

  edm::LogInfo("Alignment") << "@SUB=AlignmentParameterStore"
                            << "Created with " << theAlignables.size() << " alignables.";
}

//__________________________________________________________________________________________________
AlignmentParameterStore::~AlignmentParameterStore()
{
  delete theCorrelationsStore;
  delete theTrackerAlignableId;
}

//__________________________________________________________________________________________________
CompositeAlignmentParameters
AlignmentParameterStore::selectParameters( const std::vector<AlignableDet*>& alignabledets ) const
{

  std::vector<Alignable*> alignables;
  std::map <AlignableDet*,Alignable*> alidettoalimap;
  std::map <Alignable*,int> aliposmap;
  std::map <Alignable*,int> alilenmap;
  int nparam=0;

  // iterate over AlignableDet's
  std::vector<AlignableDet*>::const_iterator iad;
  for( iad = alignabledets.begin(); iad != alignabledets.end(); ++iad ) 
  {
    Alignable* ali = alignableFromAlignableDet( *iad );
    if ( ali ) 
    {
      alidettoalimap[ *iad ] = ali; // Add to map
      // Check if Alignable already there, insert into vector if not
      if ( find(alignables.begin(),alignables.end(),ali) == alignables.end() ) 
      {
	alignables.push_back(ali);
	AlignmentParameters* ap = ali->alignmentParameters();
	nparam += ap->numSelected();
      }
    }
  }

  AlgebraicVector* selpar = new AlgebraicVector( nparam, 0 );
  AlgebraicSymMatrix* selcov = new AlgebraicSymMatrix( nparam, 0 );

  // Fill in parameters and corresponding covariance matricess
  int ipos = 1; // NOTE: .sub indices start from 1
  std::vector<Alignable*>::const_iterator it1;
  for( it1 = alignables.begin(); it1 != alignables.end(); ++it1 ) 
  {
    AlignmentParameters* ap = (*it1)->alignmentParameters();
    selpar->sub( ipos, ap->selectedParameters() );
    selcov->sub( ipos, ap->selectedCovariance() );
    int npar = ap->numSelected();
    aliposmap[*it1]=ipos;
    alilenmap[*it1]=npar;
    ipos +=npar;
  }

  // Fill in the correlations. Has to be an extra loop, because the
  // AlignmentExtendedCorrelationsStore (if used) needs the
  // alignables' covariance matrices already present.
  ipos = 1;
  for( it1 = alignables.begin(); it1 != alignables.end(); ++it1 ) 
  {
    int jpos=1;

    // Look for correlations between alignables
    std::vector<Alignable*>::const_iterator it2;
    for( it2 = alignables.begin(); it2 != it1; ++it2 ) 
    {
      theCorrelationsStore->correlations( *it1, *it2, *selcov, ipos-1, jpos-1 );
      jpos += (*it2)->alignmentParameters()->numSelected();
    }

    ipos += (*it1)->alignmentParameters()->numSelected();
  }

  AlignmentParametersData::DataContainer data( new AlignmentParametersData( selpar, selcov ) );
  CompositeAlignmentParameters aap( data, alignables, alidettoalimap, aliposmap, alilenmap );

  return aap;
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::updateParameters( const CompositeAlignmentParameters& aap )
{

  std::vector<Alignable*> alignables = aap.components();
  const AlgebraicVector& parameters = aap.parameters();
  const AlgebraicSymMatrix& covariance = aap.covariance();

  int ipos = 1; // NOTE: .sub indices start from 1

  // Loop over alignables
  for( std::vector<Alignable*>::const_iterator it=alignables.begin(); it != alignables.end(); ++it ) 
  {
    // Update parameters and local covariance   
    AlignmentParameters* ap = (*it)->alignmentParameters();
    int nsel = ap->numSelected();
    AlgebraicVector subvec = parameters.sub( ipos, ipos+nsel-1 );
    AlgebraicSymMatrix subcov = covariance.sub( ipos, ipos+nsel-1 );
    AlignmentParameters* apnew = ap->cloneFromSelected( subvec, subcov );
    (*it)->setAlignmentParameters( apnew );
	  
    // Now update correlations between detectors
    int jpos = 1;
    for( std::vector<Alignable*>::const_iterator it2 = alignables.begin(); it2 != it; ++it2 ) 
    {
      theCorrelationsStore->setCorrelations( *it, *it2, covariance, ipos-1, jpos-1 );
      jpos += (*it2)->alignmentParameters()->numSelected();
    }

    ipos+=nsel;
  }

}


//__________________________________________________________________________________________________
std::vector<Alignable*> AlignmentParameterStore::validAlignables(void) const
{ 
  std::vector<Alignable*> result;
  for (std::vector<Alignable*>::const_iterator iali = theAlignables.begin();
       iali != theAlignables.end(); ++iali)
    if ( (*iali)->alignmentParameters()->isValid() ) result.push_back(*iali);

  LogDebug("Alignment") << "@SUB=AlignmentParameterStore::validAlignables"
                        << "Valid alignables: " << result.size()
                        << "out of " << theAlignables.size();
  return result;
}

//__________________________________________________________________________________________________
Alignable* AlignmentParameterStore::alignableFromAlignableDet( AlignableDet* alignableDet ) const
{
  Alignable *mother = alignableDet;
  while (mother) {
    if (mother->alignmentParameters()) return mother;
    mother = mother->mother();
  }

  return 0;
}

//__________________________________________________________________________________________________
void AlignmentParameterStore::applyParameters(void)
{
  std::vector<Alignable*>::const_iterator iali;
  for ( iali = theAlignables.begin(); iali != theAlignables.end(); ++iali) 
    applyParameters( *iali );
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::applyParameters(Alignable* alignable)
{

  // Get alignment parameters
  RigidBodyAlignmentParameters* ap = 
    dynamic_cast<RigidBodyAlignmentParameters*>( alignable->alignmentParameters() );

  if ( !ap )
    throw cms::Exception("BadAlignable") 
      << "applyParameters: provided alignable does not have rigid body alignment parameters";

  // Translation in local frame
  AlgebraicVector shift = ap->translation();

  // Translation local->global
  LocalPoint l0 = Local3DPoint( 0.0,  0.0, 0.0);
  LocalPoint l1 = Local3DPoint(shift[0], shift[1], shift[2]);
  GlobalPoint g0 = alignable->surface().toGlobal( l0 );
  GlobalPoint g1 = alignable->surface().toGlobal( l1 );
  GlobalVector dg = g1-g0;
  alignable->move(dg);

  // Rotation in local frame
  AlgebraicVector rota = ap->rotation();
  if ( fabs(rota[0]) > 1e-5 || fabs(rota[1]) > 1e-5 || fabs(rota[2]) > 1e-5 ) 
  {
    AlignmentTransformations alignTransform;
    Surface::RotationType rot = alignTransform.rotationType( alignTransform.rotMatrix3(rota) );
    Surface::RotationType rot2 =
      alignTransform.localToGlobalMatrix( rot, alignable->globalRotation() );
    alignable->rotateInGlobalFrame(rot2);
  }
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::resetParameters(void)
{
  // Erase contents of correlation map
  theCorrelationsStore->resetCorrelations();

  // Iterate over alignables in the store and reset parameters
  std::vector<Alignable*>::const_iterator iali;
  for ( iali = theAlignables.begin(); iali != theAlignables.end(); ++iali )
    resetParameters( *iali );
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::resetParameters( Alignable* ali )
{
  if ( ali ) 
  {
    // Get alignment parameters for this alignable
    AlignmentParameters* ap = ali->alignmentParameters();
    if ( ap ) 
    {
      int npar=ap->numSelected();
          
      AlgebraicVector par(npar,0);
      AlgebraicSymMatrix cov(npar,0);
      AlignmentParameters* apnew = ap->cloneFromSelected(par,cov);
      ali->setAlignmentParameters(apnew);
      apnew->setValid(false);
    }
    else 
      edm::LogError("BadArgument") << "@SUB=AlignmentParameterStore::resetParameters"
				   << "alignable has no alignment parameter";
  }
  else
    edm::LogError("BadArgument") << "@SUB=AlignmentParameterStore::resetParameters"
                                 << "argument is NULL";
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::acquireRelativeParameters(void)
{

  AlignmentTransformations alignTransform;
  std::vector<Alignable*>::const_iterator iali;
  for ( iali = theAlignables.begin(); iali != theAlignables.end(); ++iali ) 
  {
    RigidBodyAlignmentParameters* ap = 
      dynamic_cast<RigidBodyAlignmentParameters*>( (*iali)->alignmentParameters() );

    if ( !ap )
      throw cms::Exception("BadAlignable") 
	<< "acquireRelativeParameters: "
	<< "provided alignable does not have rigid body alignment parameters";

    AlgebraicVector par( ap->size(),0 );
    AlgebraicSymMatrix cov( ap->size(), 0 );
	  
    // Get displacement and transform global->local
    LocalVector dloc = (*iali)->surface().toLocal( (*iali)->displacement() );
    par[0]=dloc.x();
    par[1]=dloc.y();
    par[2]=dloc.z();
	  
    // Global rel rotation
    Surface::RotationType rot = (*iali)->rotation();
    // Global abs rotation
    Surface::RotationType detrot = (*iali)->surface().rotation();

    // Global euler angles
    AlgebraicVector euglob = alignTransform.eulerAngles( rot,0 );

    // Transform to local euler angles
    AlgebraicVector euloc = alignTransform.globalToLocalEulerAngles( euglob, detrot );
    par[3]=euloc[0];
    par[4]=euloc[1];
    par[5]=euloc[2];
	  
    // Clone parameters
    RigidBodyAlignmentParameters* apnew = ap->clone(par,cov);
	  
    (*iali)->setAlignmentParameters(apnew);
  }
}


//__________________________________________________________________________________________________
// Get type/layer from Alignable
// type: -6   -5   -4   -3   -2    -1     1     2    3    4    5    6
//      TEC- TOB- TID- TIB- PxEC- PxBR- PxBr+ PxEC+ TIB+ TID+ TOB+ TEC+
// Layers start from zero
std::pair<int,int> AlignmentParameterStore::typeAndLayer(const Alignable* ali) const
{
  return theTrackerAlignableId->typeAndLayerFromAlignable( ali );
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::
applyAlignableAbsolutePositions( const Alignables& alivec, 
                                 const AlignablePositions& newpos, 
                                 int& ierr )
{
  unsigned int nappl=0;
  ierr=0;

  // Iterate over list of alignables
  for ( Alignables::const_iterator iali = alivec.begin(); iali != alivec.end(); ++iali ) 
  {
    Alignable* ali = *iali;
    unsigned int detId = theTrackerAlignableId->alignableId(ali);
    int typeId = theTrackerAlignableId->alignableTypeId(ali);

    // Find corresponding entry in AlignablePositions
    bool found=false;
    for ( AlignablePositions::const_iterator ipos = newpos.begin(); ipos != newpos.end(); ++ipos ) 
      if ( detId == ipos->id() && typeId == ipos->objId() ) 
	if ( found )
	  edm::LogError("DuplicatePosition")
	    << "New positions for alignable found more than once!";
	else
	{
	  // New position/rotation
	  GlobalPoint pnew = ipos->pos();
	  Surface::RotationType rnew = ipos->rot();
	  // Current position / rotation
	  GlobalPoint pold = ali->surface().position();
	  Surface::RotationType rold = ali->surface().rotation();
				
	  // shift needed to move from current to new position
	  GlobalVector shift = pnew - pold;
	  ali->move( shift );
	  LogDebug("NewPosition") << "moving by" << shift;
				
	  // Delta-rotation needed to rotate from current to new rotation
	  int ierr;
	  AlignmentTransformations alignTransform;
	  Surface::RotationType rot = 
	    alignTransform.rotationType(alignTransform.algebraicMatrix(rold).inverse(ierr)) * rnew;
	  if ( ierr )
	    edm::LogError("InversionError") << "Matrix inversion failed: not rotating";
	  else
	    { 
	      // 'Repair' matrix for rounding errors 
	      Surface::RotationType rotfixed = alignTransform.rectify(rot);
	      ali->rotateInGlobalFrame(rotfixed);
	      AlgebraicMatrix mrot = alignTransform.algebraicMatrix( rotfixed );
	      LogDebug("NewRotation") << "rotating by: " << mrot;
	    }
				
	  // add position error
	  // AlignmentPositionError ape(shift.x(),shift.y(),shift.z());
	  // (*iali)->addAlignmentPositionError(ape);
	  // (*iali)->addAlignmentPositionErrorFromRotation(rot);
				
	  found=true;
	  ++nappl;
	}
  }

  if ( nappl< newpos.size() )
    edm::LogError("Mismatch") << "Applied only " << nappl << " new positions" 
			      << " out of " << newpos.size();

  LogDebug("NewPositions") << "Applied new positions for " << nappl
                           << " out of " << alivec.size() <<" alignables.";

}


//__________________________________________________________________________________________________
void AlignmentParameterStore::
applyAlignableRelativePositions( const Alignables& alivec, const AlignableShifts& shifts, int& ierr )
{

  unsigned int nappl=0;
  ierr=0;

  // Iterate over list of alignables
  for ( Alignables::const_iterator iali = alivec.begin(); iali != alivec.end(); ++iali) 
  {
    unsigned int detId = theTrackerAlignableId->alignableId( *iali );
    int typeId=theTrackerAlignableId->alignableTypeId( *iali );

    // Find corresponding entry in AlignableShifts
    bool found = false;
    for ( AlignableShifts::const_iterator ipos = shifts.begin(); ipos != shifts.end(); ++ipos ) 
    {
      if ( detId == ipos->id() && typeId == ipos->objId() ) 
	if ( found )
	  edm::LogError("DuplicatePosition")
	    << "New positions for alignable found more than once!";
	else
	{
	  // New position/rotation shift
	  GlobalVector pnew = ipos->pos();
	  Surface::RotationType rnew = ipos->rot();

	  (*iali)->move(pnew);
	  (*iali)->rotateInGlobalFrame(rnew);
				
	  // Add position error
	  //AlignmentPositionError ape(pnew.x(),pnew.y(),pnew.z());
	  //(*iali)->addAlignmentPositionError(ape);
	  //(*iali)->addAlignmentPositionErrorFromRotation(rnew);

	  found=true;
	  ++nappl;
	}
    }
  }
  
  if ( nappl < shifts.size() )
    edm::LogError("Mismatch") << "Applied only " << nappl << " new positions" 
			      << " out of " << shifts.size();

  LogDebug("NewPositions") << "Applied new positions for " << nappl << " alignables.";
}



//__________________________________________________________________________________________________
void AlignmentParameterStore::attachAlignmentParameters( const Parameters& parvec, int& ierr )
{
  attachAlignmentParameters( theAlignables, parvec, ierr);
}



//__________________________________________________________________________________________________
void AlignmentParameterStore::attachAlignmentParameters( const Alignables& alivec, 
                                                         const Parameters& parvec, int& ierr )
{
  int ipass = 0;
  int ifail = 0;
  ierr = 0;

  // Iterate over alignables
  for ( Alignables::const_iterator iali = alivec.begin(); iali != alivec.end(); ++iali ) 
  {
    // Iterate over Parameters
    bool found=false;
    for ( Parameters::const_iterator ipar = parvec.begin(); ipar != parvec.end(); ++ipar) 
    {
      // Get new alignment parameters
      RigidBodyAlignmentParameters* ap = dynamic_cast<RigidBodyAlignmentParameters*>(*ipar); 

      // Check if parameters belong to alignable 
      if ( ap->alignable() == (*iali) )
      {
	if (!found) 
	{
          (*iali)->setAlignmentParameters(ap);
          ++ipass;
          found=true;
        } 
        else edm::LogError("DuplicateParameters") << "More than one parameters for Alignable";
      }
    }
    if (!found) ++ifail;
  }
  if (ifail>0) ierr=-1;
  
  LogDebug("attachAlignmentParameters") << " Parameters, Alignables: " << parvec.size() << ","
                                        << alivec.size() << "\n pass,fail: " << ipass << ","<< ifail;
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::attachCorrelations( const Correlations& cormap, 
                                                  bool overwrite, int& ierr )
{
  attachCorrelations( theAlignables, cormap, overwrite, ierr );
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::attachCorrelations( const Alignables& alivec, 
                                                  const Correlations& cormap, 
                                                  bool overwrite, int& ierr )
{
  ierr=0;
  int icount=0;

  // Iterate over correlations
  for ( Correlations::const_iterator icor = cormap.begin(); icor!=cormap.end(); ++icor ) 
  {
    AlgebraicMatrix mat=(*icor).second;
    Alignable* ali1 = (*icor).first.first;
    Alignable* ali2 = (*icor).first.second;

    // Check if alignables exist
    if ( find( alivec.begin(), alivec.end(), ali1 ) != alivec.end() && 
         find( alivec.begin(), alivec.end(), ali2 ) != alivec.end() )
    {
      // Check if correlations already existing between these alignables
      if ( !theCorrelationsStore->correlationsAvailable(ali1,ali2) || (overwrite) ) 
       {
         theCorrelationsStore->setCorrelations(ali1,ali2,mat);
         ++icount;
       }
      else edm::LogInfo("AlreadyExists") << "Correlation existing and not overwritten";
    }
    else edm::LogInfo("IgnoreCorrelation") << "Ignoring correlation with no alignables!";
  }

  LogDebug( "attachCorrelations" ) << " Alignables,Correlations: " << alivec.size() <<","<< cormap.size() 
                                   << "\n applied: " << icount ;

}


//__________________________________________________________________________________________________
void AlignmentParameterStore::
attachUserVariables( const Alignables& alivec,
                     const std::vector<AlignmentUserVariables*>& uvarvec, int& ierr )
{
  ierr=0;

  LogDebug("DumpArguments") << "size of alivec:   "  << alivec.size()
                            << "\nsize of uvarvec: " << uvarvec.size();

  std::vector<AlignmentUserVariables*>::const_iterator iuvar=uvarvec.begin();

  for ( Alignables::const_iterator iali=alivec.begin(); iali!=alivec.end(); ++iali, ++iuvar ) 
  {
    AlignmentParameters* ap = (*iali)->alignmentParameters();
    AlignmentUserVariables* uvarnew = (*iuvar);
    ap->setUserVariables(uvarnew);
  }
}


//__________________________________________________________________________________________________
void AlignmentParameterStore::setAlignmentPositionError( const Alignables& alivec, 
                                                         double valshift, double valrot )
{
  bool first=true;
  for ( Alignables::const_iterator iali = alivec.begin(); iali != alivec.end(); ++iali ) 
  {

    // First reset APE	 
    AlignmentPositionError nulApe(0,0,0);	 
    (*iali)->setAlignmentPositionError(nulApe);

    // Set APE from displacement
    AlignmentPositionError ape(valshift,valshift,valshift);
    if ( valshift > 0. ) (*iali)->addAlignmentPositionError(ape);
    else (*iali)->setAlignmentPositionError(ape);
    if (first) LogDebug("StoreAPE") << "Store APE from shift: " << valshift;

    // Set APE from rotation
    AlignmentTransformations alignTransform;
    AlgebraicVector r(3);
    r[0]=valrot; r[1]=valrot; r[2]=valrot;
    Surface::RotationType aperot = alignTransform.rotationType( alignTransform.rotMatrix3(r) );
    (*iali)->addAlignmentPositionErrorFromRotation(aperot);
    if (first) LogDebug("StoreAPE") << "Store APE from rotation: " << valrot;

    first=false;
  }
}
