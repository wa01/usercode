#include "HLTrigger/Bphysics/interface/RecoChargedCandidatesToTracks.h"


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// system include files
#include <memory>


//
// constructors and destructor
//
RecoChargedCandidatesToTracks::RecoChargedCandidatesToTracks(const edm::ParameterSet& iConfig) :
  candTag_(iConfig.getParameter<edm::InputTag>("candTag"))
{
   //register your products
   produces<reco::TrackCollection>();
}


RecoChargedCandidatesToTracks::~RecoChargedCandidatesToTracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
RecoChargedCandidatesToTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // JPsi filter object collection
  //
  edm::Handle<trigger::TriggerFilterObjectWithRefs> jpsiHandle;
  iEvent.getByLabel(candTag_,jpsiHandle);

  std::vector<reco::RecoChargedCandidateRef> trackRefs;
  jpsiHandle->getObjects(trigger::TriggerTrack,trackRefs);

  std::vector<reco::RecoChargedCandidateRef> uniqueRefs;
  uniqueRefs.reserve(trackRefs.size());
  for ( unsigned int i=0; i<trackRefs.size(); ++i ) {
    if ( std::find(uniqueRefs.begin(),uniqueRefs.end(),trackRefs[i])==uniqueRefs.end() )
      uniqueRefs.push_back(trackRefs[i]);
  }
  //
  // product
  //
  std::auto_ptr<reco::TrackCollection> product(new reco::TrackCollection);
  product->reserve(uniqueRefs.size());

//   std::cout << "Copying " << uniqueRefs.size() << " tracks out of " << trackRefs.size()
// 	    << " candidates" << std::endl;

  for ( unsigned int i=0; i<uniqueRefs.size(); ++i ) {
    product->push_back(*uniqueRefs[i]->track());
  }

  iEvent.put(product);
}

// ------------ method called once each job just before starting event loop  ------------
void 
RecoChargedCandidatesToTracks::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecoChargedCandidatesToTracks::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoChargedCandidatesToTracks);
