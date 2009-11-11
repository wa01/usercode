#include "HLTrigger/Bphysics/interface/HLTOniaCandidateFilter.h"

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
HLTOniaCandidateFilter::HLTOniaCandidateFilter(const edm::ParameterSet& iConfig) :
  candTag_(iConfig.getParameter<edm::InputTag>("candTag")),
  saveTag_(iConfig.getParameter<bool>("saveTag")),
  minMass_(iConfig.getParameter<double>("MinMass")),
  maxMass_(iConfig.getParameter<double>("MaxMass")),
  checkCharge_(iConfig.getParameter<bool>("checkCharge"))
{
  LogDebug("HLTOniaCandidateFilter") << "instantiated with parameters\n"
				     << "  candTag  = " << candTag_ << "\n"
				     << "  saveTag = " << saveTag_ << "\n"
				     << "  MinMass  = " << minMass_ << "\n"
				     << "  MaxMass  = " << maxMass_ << "\n"
				     << "  checkCharge = " << checkCharge_;
  //register your products
  produces<trigger::TriggerFilterObjectWithRefs>();
}


HLTOniaCandidateFilter::~HLTOniaCandidateFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HLTOniaCandidateFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 // The filter object
  std::auto_ptr<trigger::TriggerFilterObjectWithRefs>
    filterproduct (new trigger::TriggerFilterObjectWithRefs(path(),module()));
  if ( saveTag_ ) {
    filterproduct->addCollectionTag(candTag_);
  }
  //
  // Onia candidates
  //
  edm::Handle<reco::DiRefChargedCandidateCollection> candHandle;
  iEvent.getByLabel(candTag_,candHandle);

  unsigned int nQ(0);
  unsigned int nSel(0);
  reco::Particle::LorentzVector p4Muon;
  reco::Particle::LorentzVector p4JPsi;
  for ( unsigned int i=0; i<candHandle->size(); ++i ) {
    const reco::DiRefChargedCandidate& candidate = (*candHandle)[i];
    if ( checkCharge_ && candidate.charge()!=0 )  continue;
      ++nQ;
      double mass = candidate.p4().mass();
      if ( mass>minMass_ && mass<maxMass_ ) {
	++nSel;
	filterproduct->addObject(trigger::TriggerMuon,candidate.firstDaughterRef());
	filterproduct->addObject(trigger::TriggerTrack,candidate.secondDaughterRef());
      }
  }


  if ( edm::isDebugEnabled() ) {
    std::ostringstream stream;
    stream << "Total number of candidates = " << candHandle->size()
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
      LogDebug("HLTOniaCandidateFilter") << stream.str();
    }
  }

  iEvent.put(filterproduct);

  return nSel>0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HLTOniaCandidateFilter);
