#include "HLTrigger/Bphysics/interface/QuarkoniaTrackSelector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// system include files
#include <memory>
#include <iostream>
#include <sstream>

//
// constructors and destructor
//
QuarkoniaTrackSelector::QuarkoniaTrackSelector(const edm::ParameterSet& iConfig) :
  muonTag_(iConfig.getParameter<edm::InputTag>("muonCandidates")),
  trackTag_(iConfig.getParameter<edm::InputTag>("trackCandidates")),
  minMasses_(iConfig.getParameter< std::vector<double> >("MinMasses")),
  maxMasses_(iConfig.getParameter< std::vector<double> >("MaxMasses")),
  checkCharge_(iConfig.getParameter<bool>("checkCharge"))
{
  LogDebug("HLTTrackFilter") << "instantiated with parameters\n"
			     << "  muonTag  = " << muonTag_ << "\n"
			     << "  trackTag = " << trackTag_ << "\n"
// 			     << "  MinMass  = " << minMass_ << "\n"
// 			     << "  MaxMass  = " << maxMass_ << "\n"
			     << "  checkCharge  = " << checkCharge_;
  //register your products
  produces<reco::TrackCollection>();
  //
  // verify mass windows
  //
  bool massesValid = minMasses_.size()==maxMasses_.size();
  if ( massesValid ) {
    for ( unsigned int i=0; i<minMasses_.size(); ++i ) {
      if ( minMasses_[i]<0 || maxMasses_[i]<0 || 
	   minMasses_[i]>maxMasses_[i] )  massesValid = false;
    }
  }
  if ( !massesValid ) {
    edm::LogError("QuarkoniaTrackSelector") << "Inconsistency in definition of mass windows, "
					    << "no track will be selected";
    minMasses_.clear();
    maxMasses_.clear();
  }
}


QuarkoniaTrackSelector::~QuarkoniaTrackSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
QuarkoniaTrackSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // the product
  //
  std::auto_ptr<reco::TrackCollection> product(new reco::TrackCollection);
  //
  // Muons
  //
  edm::Handle<reco::RecoChargedCandidateCollection> muonHandle;
  iEvent.getByLabel(muonTag_,muonHandle);
  //
  // Tracks
  //
  edm::Handle<reco::RecoChargedCandidateCollection> trackHandle;
  iEvent.getByLabel(trackTag_,trackHandle);
  //
  // Verification
  //
  if ( !muonHandle.isValid() || !trackHandle.isValid() || minMasses_.empty() ) {
    iEvent.put(product);
    return;
  }
  //
  // combinations
  //
  unsigned int nQ(0);
  unsigned int nComb(0);
  std::vector<reco::TrackRef> selectedTracks;
  selectedTracks.reserve(muonHandle->size());
  reco::Particle::LorentzVector p4Muon;
  reco::Particle::LorentzVector p4JPsi;
  for ( unsigned int im=0; im<muonHandle->size(); ++im ) {
    const reco::RecoChargedCandidate& muon = (*muonHandle)[im];
    int qMuon = muon.charge();
    p4Muon = muon.p4();
    for ( unsigned int it=0; it<trackHandle->size(); ++it ) {
      const reco::RecoChargedCandidate& track = (*trackHandle)[it];
      if ( checkCharge_ && track.charge()!=-qMuon )  continue;
      ++nQ;
      double mass = (p4Muon+track.p4()).mass();
      for ( unsigned int j=0; j<minMasses_.size(); ++j ) {
	if ( mass>minMasses_[j] && mass<maxMasses_[j] ) {
	  ++nComb;
	  if ( find(selectedTracks.begin(),selectedTracks.end(),track.track())==
	       selectedTracks.end() )  selectedTracks.push_back(track.track());
	  break;
	}
      }
    }
  }


  for ( unsigned int i=0; i<selectedTracks.size(); ++i )  product->push_back(*selectedTracks[i]);

  if ( edm::isDebugEnabled() ) {
    std::ostringstream stream;
    stream << "Total number of combinations = " << muonHandle->size()*trackHandle->size()
	   << " , after charge " << nQ << " , after mass " << nComb << std::endl;
    stream << "Selected " << product->size() << " tracks with # / q / pt / eta\n";
    for ( unsigned int i=0; i<product->size(); ++i ) {
      const reco::Track& track = (*product)[i];
      stream << "  " << i << " " << track.charge() << " "
	     << track.pt() << " " << track.eta() << "\n";
    }
    LogDebug("QuarkoniaTrackSelector") << stream.str();
  }

  iEvent.put(product);
}


//define this as a plug-in
DEFINE_FWK_MODULE(QuarkoniaTrackSelector);
