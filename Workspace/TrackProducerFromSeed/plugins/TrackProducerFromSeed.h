#ifndef TrackProducerFromSeed_h
#define TrackProducerFromSeed_h

#include "Workspace/TrackProducerFromSeed/interface/KfTrackProducerFromSeedBase.h"
#include "Workspace/TrackProducerFromSeed/interface/TrackProducerFromSeedAlgorithm.h"

// #include "TrackingTools/TransientTrack/interface/TransientTrack.h"

class TrackProducerFromSeed : public KfTrackProducerFromSeedBase, public edm::EDProducer {
public:

  /// Constructor
  explicit TrackProducerFromSeed(const edm::ParameterSet& iConfig);

  /// Implementation of produce method
  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  TrackProducerFromSeedAlgorithm<reco::Track> theAlgo;

};

#endif
