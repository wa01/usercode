#include "Workspace/TrackProducerFromSeed/interface/TrackProducerFromSeedAlgorithm.h"

#include "DataFormats/Common/interface/OrphanHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoTracker/TrackProducer/interface/TrackingRecHitLessFromGlobalPosition.h"

#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "Utilities/General/interface/CMSexception.h"

// #include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "TrackingTools/TrackFitters/interface/RecHitSorter.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

template <> bool
TrackProducerFromSeedAlgorithm<reco::Track>::buildTrack (const Propagator * thePropagator,
						   AlgoProductCollection& algoResults,
						   TransientTrackingRecHit::RecHitContainer& hits,
						   TrajectoryStateOnSurface& theTSOS,
						   const TrajectorySeed& seed,
						   float ndof,
						   const reco::BeamSpot& bs,
 						   SeedRef seedRef,
						   int qualityMask)						 
{
//   std::cout << "MTPA::buildTrack start" << std::endl;
  //variable declarations
  reco::Track * theTrack;
  Trajectory * theTraj; 
  PropagationDirection seedDir = seed.direction();
  SetPropagationDirection setDir(*thePropagator,anyDirection);
 
  std::vector<Trajectory> trajVec(1, Trajectory(seed, thePropagator->propagationDirection()));
  Trajectory& myTraj = trajVec.front();
  myTraj.reserve(hits.size());
  for ( TransientTrackingRecHit::RecHitContainer::const_iterator ihit = hits.begin(); ihit != hits.end(); ++ihit ) {

    const TransientTrackingRecHit & hit = (**ihit);
    if ( !hit.isValid() && hit.surface() == 0) {
      continue;
    }

    TrajectoryStateOnSurface predTsos = thePropagator->propagate(theTSOS,*(hit.surface()));
    if ( !predTsos.isValid() )  continue;

    myTraj.push(TrajectoryMeasurement(predTsos,predTsos,*ihit,0,0));
  }
  if ( myTraj.empty() ) {
    edm::LogError("EmptyTrajectory") << " from " << hits.size() << " hits";
    return false;
  }

  LogDebug("TrackProducer") <<" FITTER FOUND "<< trajVec.size() << " TRAJECTORIES" <<"\n";
  TrajectoryStateOnSurface innertsos;
  
  if (trajVec.size() != 0){
    theTraj = new Trajectory( trajVec.front() );
    theTraj->setSeedRef(seedRef);
    
    if (theTraj->direction() == alongMomentum) {
      innertsos = theTraj->firstMeasurement().updatedState();
    } else { 
      innertsos = theTraj->lastMeasurement().updatedState();
    }
    
    ndof = 0;
    TransientTrackingRecHit::RecHitContainer validHits;
    theTraj->validRecHits(validHits);
    for (TransientTrackingRecHit::RecHitContainer::iterator h=validHits.begin();h!=validHits.end();++h)
      ndof = ndof + ((*h)->dimension())*((*h)->weight());
    if (theTSOS.magneticField()->inTesla(GlobalPoint(0,0,0)).mag2()<DBL_MIN) ndof = ndof - 4;
    else ndof = ndof - 5;
    
    TSCBLBuilderNoMaterial tscblBuilder;
    //    const FreeTrajectoryState & stateForProjectionToBeamLine=*innertsos.freeState();
    const TrajectoryStateOnSurface & stateForProjectionToBeamLineOnSurface = theTraj->closestMeasurement(GlobalPoint(bs.x0(),bs.y0(),bs.z0())).updatedState();
    if (!stateForProjectionToBeamLineOnSurface.isValid()){
      edm::LogError("CannotPropagateToBeamLine")<<"the state on the closest measurement isnot valid. skipping track.";
      delete theTraj;
      return false;
    }
    const FreeTrajectoryState & stateForProjectionToBeamLine=*stateForProjectionToBeamLineOnSurface.freeState();

    LogDebug("TrackProducer") << "stateForProjectionToBeamLine=" << stateForProjectionToBeamLine;

    TrajectoryStateClosestToBeamLine tscbl = tscblBuilder(stateForProjectionToBeamLine,bs);

    if (tscbl.isValid()==false) {
        delete theTraj;
        return false;
    }

    GlobalPoint v = tscbl.trackStateAtPCA().position();
    math::XYZPoint  pos( v.x(), v.y(), v.z() );
    GlobalVector p = tscbl.trackStateAtPCA().momentum();
    math::XYZVector mom( p.x(), p.y(), p.z() );

    LogDebug("TrackProducer") << "pos=" << v << " mom=" << p << " pt=" << p.perp() << " mag=" << p.mag();

    theTrack = new reco::Track(theTraj->chiSquared(),
			       int(ndof),//FIXME fix weight() in TrackingRecHit
			       pos, mom, tscbl.trackStateAtPCA().charge(), 
			       tscbl.trackStateAtPCA().curvilinearError(),
			       algo_);
   
    theTrack->setQualityMask(qualityMask);
     theTrack->setNLoops(0);

    LogDebug("TrackProducer") << "theTrack->pt()=" << theTrack->pt();

    LogDebug("TrackProducer") <<"track done\n";

    AlgoProduct aProduct(theTraj,std::make_pair(theTrack,seedDir));
    algoResults.push_back(aProduct);
    
//     std::cout << "MTPA::buildTrack end 1" << std::endl;
    return true;
  } 
  else  return false;
}
