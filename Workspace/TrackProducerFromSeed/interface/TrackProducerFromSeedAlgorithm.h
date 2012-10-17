#ifndef TrackProducerFromSeedAlgorithm_h
#define TrackProducerFromSeedAlgorithm_h

/** \class TrackProducerFromSeedAlgorithm
 *  This class calls the Final Fit and builds the Tracks then produced by the TrackProducer or by the TrackRefitter
 *
 *  $Date: 2012/03/13 16:42:53 $
 *  $Revision: 1.27 $
 *  \author cerati
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/PatternTools/interface/TrackConstraintAssociation.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

class MagneticField;
class TrackingGeometry;
class TrajectoryFitter;
class Propagator;
class Trajectory;
class TrajectoryStateOnSurface;
class TransientTrackingRecHitBuilder;


template <class T>
class TrackProducerFromSeedAlgorithm {
public:
  typedef std::vector<T> TrackCollection;
  typedef std::pair<Trajectory*, std::pair<T*,PropagationDirection> > AlgoProduct; 
  typedef std::vector< AlgoProduct >  AlgoProductCollection;
  typedef edm::RefToBase<TrajectorySeed> SeedRef;
 public:

  /// Constructor
  TrackProducerFromSeedAlgorithm(const edm::ParameterSet& conf) : 
    conf_(conf),
    algoName_(conf_.getParameter<std::string>( "AlgorithmName" )),
    algo_(reco::TrackBase::algoByName(algoName_)),
    reMatchSplitHits_(false)
      {
	if (conf_.exists("reMatchSplitHits"))
	  reMatchSplitHits_=conf_.getParameter<bool>("reMatchSplitHits");
      }

  /// Destructor
  ~TrackProducerFromSeedAlgorithm() {}
  
  /// Run the Final Fit taking TrackCandidates as input
  void runWithSeed(const TrackingGeometry *, 
		   const MagneticField *, 
		   const TrajectorySeedCollection&,
		   const Propagator *,
		   const TransientTrackingRecHitBuilder*,
		   const reco::BeamSpot&,
		   AlgoProductCollection &);


  /// Construct Tracks to be put in the event
  bool buildTrack(const Propagator *,
		  AlgoProductCollection& ,
		  TransientTrackingRecHit::RecHitContainer&,
		  TrajectoryStateOnSurface& ,
		  const TrajectorySeed&,		  
		  float,
		  const reco::BeamSpot&,
 		  SeedRef seedRef = SeedRef(),
		  int qualityMask=0);

 private:
  edm::ParameterSet conf_;  
  std::string algoName_;
  reco::TrackBase::TrackAlgorithm algo_;
  bool reMatchSplitHits_;

};

#include "Workspace/TrackProducerFromSeed/interface/TrackProducerFromSeedAlgorithm.icc"

template <> bool
TrackProducerFromSeedAlgorithm<reco::Track>::buildTrack(const Propagator *,
						  AlgoProductCollection& ,
						  TransientTrackingRecHit::RecHitContainer&,
						  TrajectoryStateOnSurface& ,
						  const TrajectorySeed&,
						  float,
						  const reco::BeamSpot&,
						  SeedRef seedRef,
						  int qualityMask);



#endif
