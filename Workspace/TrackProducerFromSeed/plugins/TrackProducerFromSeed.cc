#include "Workspace/TrackProducerFromSeed/plugins/TrackProducerFromSeed.h"
// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

TrackProducerFromSeed::TrackProducerFromSeed(const edm::ParameterSet& iConfig):
  KfTrackProducerFromSeedBase(iConfig.getParameter<bool>("TrajectoryInEvent"),
		      iConfig.getParameter<bool>("useHitsSplitting")),
  theAlgo(iConfig)
{
  setConf(iConfig);
  setSrc( iConfig.getParameter<edm::InputTag>( "src" ), iConfig.getParameter<edm::InputTag>( "beamSpot" ));
  setAlias( iConfig.getParameter<std::string>( "@module_label" ) );

  //register your products
  produces<reco::TrackCollection>().setBranchAlias( alias_ + "Tracks" );
  produces<reco::TrackExtraCollection>().setBranchAlias( alias_ + "TrackExtras" );
  produces<TrackingRecHitCollection>().setBranchAlias( alias_ + "RecHits" );
  produces<std::vector<Trajectory> >() ;
  produces<TrajTrackAssociationCollection>();

}


void TrackProducerFromSeed::produce(edm::Event& theEvent, const edm::EventSetup& setup)
{
  LogDebug("TrackProducerFromSeed") << "Analyzing event number: " << theEvent.id() << "\n";
  //
  // create empty output collections
  //
  std::auto_ptr<TrackingRecHitCollection>    outputRHColl (new TrackingRecHitCollection);
  std::auto_ptr<reco::TrackCollection>       outputTColl(new reco::TrackCollection);
  std::auto_ptr<reco::TrackExtraCollection>  outputTEColl(new reco::TrackExtraCollection);
  std::auto_ptr<std::vector<Trajectory> >    outputTrajectoryColl(new std::vector<Trajectory>);

  //
  //declare and get stuff to be retrieved from ES
  //
  edm::ESHandle<TrackerGeometry> theG;
  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<TrajectoryFitter> theFitter;
  edm::ESHandle<Propagator> thePropagator;
  edm::ESHandle<MeasurementTracker>  theMeasTk;
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  getFromES(setup,theG,theMF,thePropagator,theMeasTk,theBuilder);

  //
  //declare and get TrackColection to be retrieved from the event
  //
  AlgoProductCollection algoResults;
  edm::Handle<TrajectorySeedCollection> theTSCollection;
  reco::BeamSpot bs;
  getFromEvt(theEvent,theTSCollection,bs);
  //protect against missing product  
  if (theTSCollection.failedToGet()){
    edm::LogError("TrackProducerFromSeed") <<"could not get the TrackCandidateCollection.";} 
  else{
    LogDebug("TrackProducerFromSeed") << "run the algorithm" << "\n";
    try{  
      theAlgo.runWithSeed(theG.product(), theMF.product(), *theTSCollection, 
			  thePropagator.product(), theBuilder.product(), bs, algoResults);
    } catch (cms::Exception &e){ edm::LogError("TrackProducerFromSeed") << "cms::Exception caught during theAlgo.runWithCandidate." << "\n" << e << "\n"; throw;}
  }
  
  //put everything in the event
  putInEvt(theEvent, thePropagator.product(),theMeasTk.product(), outputRHColl, outputTColl, outputTEColl, outputTrajectoryColl, algoResults);
  LogDebug("TrackProducerFromSeed") << "end" << "\n";
}




