#ifndef TrackProducerFromSeedBase_h
#define TrackProducerFromSeedBase_h

/** \class TrackProducerFromSeedBase
 *  Base Class To Produce Tracks
 *
 *  $Date: 2010/09/29 12:38:37 $
 *  $Revision: 1.19 $
 *  \author cerati
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackerRecHit2D/interface/ClusterRemovalInfo.h"
#include <RecoTracker/MeasurementDet/interface/MeasurementTracker.h>

class Propagator;
class TrajectoryStateUpdator;
class MeasurementEstimator;
class TrackerGeometry;
class TrajectoryFitter;
class TransientTrackingRecHitBuilder;
class NavigationSchool;

template <class T>
class TrackProducerFromSeedBase {
public:
  typedef std::vector<T> TrackCollection;
  typedef std::pair<Trajectory*, std::pair<T*,PropagationDirection> > AlgoProduct;
  typedef std::vector< AlgoProduct >  AlgoProductCollection;
public:
  /// Constructor
  TrackProducerFromSeedBase(bool trajectoryInEvent = false):
     trajectoryInEvent_(trajectoryInEvent) {}

  /// Destructor
  virtual ~TrackProducerFromSeedBase();
  
  /// Get needed services from the Event Setup
  virtual void getFromES(const edm::EventSetup&,
			 edm::ESHandle<TrackerGeometry>& ,
			 edm::ESHandle<MagneticField>& ,
// 			 edm::ESHandle<TrajectoryFitter>& ,
			 edm::ESHandle<Propagator>& ,
			 edm::ESHandle<MeasurementTracker>& ,
			 edm::ESHandle<TransientTrackingRecHitBuilder>& );

  virtual void getFromEvt(edm::Event&, edm::Handle<TrajectorySeedCollection>&, reco::BeamSpot&);

  /// Method where the procduction take place. To be implemented in concrete classes
  virtual void produce(edm::Event&, const edm::EventSetup&) = 0;

  /// Set parameter set
  void setConf(edm::ParameterSet conf){conf_=conf;}

  /// set label of source collection
  void setSrc(edm::InputTag src, edm::InputTag bsSrc){src_=src;bsSrc_=bsSrc;}

  /// set the aliases of produced collections
  void setAlias(std::string alias){
    alias.erase(alias.size()-6,alias.size());
    alias_=alias;
  }

  void setSecondHitPattern(Trajectory* traj, T& track, 
			   const Propagator* prop, const MeasurementTracker* measTk );

  const edm::ParameterSet& getConf() const {return conf_;}
 private:
  edm::ParameterSet conf_;
  edm::InputTag src_;
 protected:
  std::string alias_;
  bool trajectoryInEvent_;
  edm::OrphanHandle<TrackCollection> rTracks_;
  edm::InputTag bsSrc_;

  edm::ESHandle<NavigationSchool> theSchool;

};

#include "Workspace/TrackProducerFromSeed/interface/TrackProducerFromSeedBase.icc"

#endif
