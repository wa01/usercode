#ifndef KfTrackProducerFromSeedBase_h
#define KfTrackProducerFromSeedBase_h

/** \class KfTrackProducerFromSeedBase
 *  Produce Tracks from TrackCandidates
 *
 *  $Date: 2010/09/29 12:36:04 $
 *  $Revision: 1.4 $
 *  \author cerati
 */

#include "Workspace/TrackProducerFromSeed/interface/TrackProducerFromSeedBase.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

class Trajectory;

class KfTrackProducerFromSeedBase : public TrackProducerFromSeedBase<reco::Track> {
public:

  /// Constructor
  explicit KfTrackProducerFromSeedBase(bool trajectoryInEvent, bool split) :
    TrackProducerFromSeedBase<reco::Track>(trajectoryInEvent),useSplitting(split) {}

  /// Put produced collections in the event
  virtual void putInEvt(edm::Event&,
			const Propagator* prop,
			const MeasurementTracker* measTk,
			std::auto_ptr<TrackingRecHitCollection>&,
			std::auto_ptr<reco::TrackCollection>&,
			std::auto_ptr<reco::TrackExtraCollection>&,
			std::auto_ptr<std::vector<Trajectory> >&,
			AlgoProductCollection&);


  //  void setSecondHitPattern(Trajectory* traj, reco::Track& track);
 private:
  bool useSplitting;

};

#endif
