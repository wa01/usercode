#include "HLTrigger/Bphysics/interface/HLTOniaMuonTrackMassFilter.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// system include files
#include <memory>
#include <iostream>
#include <sstream>

//
// constructors and destructor
//
HLTOniaMuonTrackMassFilter::HLTOniaMuonTrackMassFilter(const edm::ParameterSet& iConfig) :
  muonTag_(iConfig.getParameter<edm::InputTag>("muonCandidates")),
  trackTag_(iConfig.getParameter<edm::InputTag>("trackCandidates")),
  saveTag_(iConfig.getParameter<bool>("saveTag")),
  minMass_(iConfig.getParameter<double>("MinMass")),
  maxMass_(iConfig.getParameter<double>("MaxMass")),
  checkCharge_(iConfig.getParameter<bool>("checkCharge"))
{
  LogDebug("HLTTrackFilter") << "instantiated with parameters\n"
			     << "  muonTag  = " << muonTag_ << "\n"
			     << "  trackTag = " << trackTag_ << "\n"
			     << "  saveTag = " << saveTag_ << "\n"
			     << "  MinMass  = " << minMass_ << "\n"
			     << "  MaxMass  = " << maxMass_ << "\n"
			     << "  checkCharge  = " << checkCharge_;
  //register your products
  produces<trigger::TriggerFilterObjectWithRefs>();
}


HLTOniaMuonTrackMassFilter::~HLTOniaMuonTrackMassFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HLTOniaMuonTrackMassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 // The filter object
  std::auto_ptr<trigger::TriggerFilterObjectWithRefs>
    filterproduct (new trigger::TriggerFilterObjectWithRefs(path(),module()));
  if ( saveTag_ ) {
    filterproduct->addCollectionTag(muonTag_);
    filterproduct->addCollectionTag(trackTag_);
  }
  //
  // Muons
  //
  edm::Handle<trigger::TriggerFilterObjectWithRefs> muonHandle;
  iEvent.getByLabel(muonTag_,muonHandle);
  //
  // Tracks
  //
  edm::Handle<trigger::TriggerFilterObjectWithRefs> trackHandle;
  iEvent.getByLabel(trackTag_,trackHandle);

  std::vector<reco::RecoChargedCandidateRef> muonRefs;
  muonHandle->getObjects(trigger::TriggerMuon,muonRefs);

  std::vector<reco::RecoChargedCandidateRef> trackRefs;
  trackHandle->getObjects(trigger::TriggerTrack,trackRefs);

  unsigned int nQ(0);
  unsigned int nSel(0);
  reco::Particle::LorentzVector p4Muon;
  reco::Particle::LorentzVector p4JPsi;
  for ( unsigned int im=0; im<muonRefs.size(); ++im ) {
    int qMuon = muonRefs[im]->charge();
    p4Muon = muonRefs[im]->p4();
    for ( unsigned int it=0; it<trackRefs.size(); ++it ) {
      if ( checkCharge_ && trackRefs[it]->charge()!=-qMuon )  continue;
      ++nQ;
      double mass = (p4Muon+trackRefs[it]->p4()).mass();
      if ( mass>minMass_ && mass<maxMass_ ) {
	++nSel;
// 	reco::CompositeCandidate jpsiCand(qJPsi,p4JPsi);
// 	jpsiCand.addDaughter(*muonRefs[im]);
// 	jpsiCand.addDaughter(*trackRefs[it]);
	filterproduct->addObject(trigger::TriggerMuon,muonRefs[im]);
	filterproduct->addObject(trigger::TriggerTrack,trackRefs[it]);
      }
    }
  }


  if ( edm::isDebugEnabled() ) {
    std::ostringstream stream;
    stream << "Total number of combinations = " << muonRefs.size()*trackRefs.size()
	   << " , after charge " << nQ << " , after mass " << nSel << std::endl;
    stream << "Found " << nSel << " jpsi candidates with # / mass / q / pt / eta" << std::endl;
    std::vector<reco::RecoChargedCandidateRef> muRefs;
    std::vector<reco::RecoChargedCandidateRef> tkRefs;
    filterproduct->getObjects(trigger::TriggerMuon,muRefs);
    filterproduct->getObjects(trigger::TriggerTrack,tkRefs);
    reco::Particle::LorentzVector p4Mu;
    reco::Particle::LorentzVector p4Tk;
    reco::Particle::LorentzVector p4JPsi;
    if ( muRefs.size()==tkRefs.size() ) {
      for ( unsigned int i=0; i<muRefs.size(); ++i ) {
	p4Mu = muRefs[i]->p4();
	p4Tk = tkRefs[i]->p4();
	p4JPsi = p4Mu + p4Tk;
	stream << "  " << i << " " 
	       << p4JPsi.M() << " "
	       << muRefs[i]->charge()+tkRefs[i]->charge() << " "
	       << p4JPsi.P() << " "
	       << p4JPsi.Eta() << "\n";
      }
      LogDebug("HLTOniaMuonTrackMassFilter") << stream.str();
    }
    else {
      LogDebug("HLTOniaMuonTrackMassFilter") << "different sizes for muon and track containers!!!";
    }
  }

  iEvent.put(filterproduct);

  return nSel>0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HLTOniaMuonTrackMassFilter);
