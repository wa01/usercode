#include "HLTrigger/Bphysics/interface/HLTCompositeMuonTrackProducer.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/RecoCandidate/interface/DiRefChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/DiRefChargedCandidateFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// system include files
#include <memory>
#include <iostream>
#include <sstream>

//
// constructors and destructor
//
HLTCompositeMuonTrackProducer::HLTCompositeMuonTrackProducer(const edm::ParameterSet& iConfig) :
  muonTag_(iConfig.getParameter<edm::InputTag>("muonCandidates")),
  trackTag_(iConfig.getParameter<edm::InputTag>("trackCandidates")),
  minMass_(iConfig.getParameter<double>("MinMass")),
  maxMass_(iConfig.getParameter<double>("MaxMass")),
  checkCharge_(iConfig.getParameter<bool>("checkCharge"))
{
  LogDebug("HLTCompositeMuonTrackProducer") << "instantiated with parameters\n"
					    << "  muonTag  = " << muonTag_ << "\n"
					    << "  trackTag = " << trackTag_ << "\n"
					    << "  MinMass  = " << minMass_ << "\n"
					    << "  MaxMass  = " << maxMass_ << "\n"
					    << "  checkCharge = " << checkCharge_;
  //register your products
  produces<reco::DiRefChargedCandidateCollection>();
}


HLTCompositeMuonTrackProducer::~HLTCompositeMuonTrackProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
HLTCompositeMuonTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 // The filter object
  std::auto_ptr<reco::DiRefChargedCandidateCollection>
    product(new reco::DiRefChargedCandidateCollection);
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


  unsigned int nQ(0);
  unsigned int nSel(0);
  reco::Particle::LorentzVector p4Muon;
  int qJPsi;
  reco::Particle::LorentzVector p4JPsi;
  for ( unsigned int im=0; im<muonHandle->size(); ++im ) {
    reco::RecoChargedCandidateRef muonRef(muonHandle,im);
    int qMuon = muonRef->charge();
    p4Muon = muonRef->p4();
    for ( unsigned int it=0; it<trackHandle->size(); ++it ) {
      reco::RecoChargedCandidateRef trackRef(trackHandle,it);
      if ( checkCharge_ && trackRef->charge()!=-qMuon )  continue;
      ++nQ;
      qJPsi = qMuon + trackRef->charge();
      p4JPsi = p4Muon + trackRef->p4();
      double mass = p4JPsi.mass();
      if ( mass>minMass_ && mass<maxMass_ ) {
	++nSel;
 	reco::DiRefChargedCandidate jpsiCand(qJPsi,p4JPsi);
 	jpsiCand.setFirstDaughter(muonRef);
 	jpsiCand.setSecondDaughter(trackRef);
	product->push_back(jpsiCand);
      }
    }
  }


  LogDebug("HLTCompositeMuonTrackProducer") << "Total number of combinations = " 
					    << muonHandle->size()*trackHandle->size()
					    << " , after charge " << nQ
					    << " , after mass " << nSel;
  if ( edm::isDebugEnabled() ) {
    std::ostringstream stream;
    stream << "Found " << nSel << " composite candidates with # / mass / q / pt / eta" << std::endl;
    for ( unsigned int i=0; i<product->size(); ++i ) {
      stream << "  " << i << " " 
	     << (*product)[i].p4().M() << " "
	     << (*product)[i].charge() << " "
	     << (*product)[i].p4().P() << " "
	     << (*product)[i].p4().Eta() << "\n";
    }
    LogDebug("HLTCompositeMuonTrackProducer") << stream.str();
  }

  iEvent.put(product);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HLTCompositeMuonTrackProducer);
