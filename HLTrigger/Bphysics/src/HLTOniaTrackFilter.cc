#include "HLTrigger/Bphysics/interface/HLTOniaTrackFilter.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include <memory>
#include <iostream>
#include <sstream>

HLTOniaTrackFilter::HLTOniaTrackFilter(const edm::ParameterSet& iConfig) :
  candTag_(iConfig.getParameter<edm::InputTag>("candTag")), 
  saveTag_(iConfig.getParameter<bool>("saveTag")),
  minPt_(iConfig.getParameter<double>("MinPt")),
  minP_(iConfig.getParameter<double>("MinP")),
  maxEta_(iConfig.getParameter<double>("MaxEta"))
{
  LogDebug("HLTOniaTrackFilter") << "instantiated with parameters\n"
			     << "  candTag = " << candTag_ << "\n"
			     << "  MinPt   = " << minPt_ << "\n"
			     << "  MinP    = " << minP_ << "\n"
			     << "  MaxEta  = " << maxEta_ << "\n"
			     << "  saveTag = " << saveTag_;

  //register your products
  produces<trigger::TriggerFilterObjectWithRefs>();
}

bool 
HLTOniaTrackFilter::filter (edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // The Filter object
  std::auto_ptr<trigger::TriggerFilterObjectWithRefs> 
    filterproduct (new trigger::TriggerFilterObjectWithRefs(path(),module()));
  if ( saveTag_ )  filterproduct->addCollectionTag(candTag_);

  // get hold of candidates
  edm::Handle<reco::RecoChargedCandidateCollection> trackHandle;
  iEvent.getByLabel(candTag_,trackHandle);
  LogDebug("HLTOniaTrackFilter") << "Total number of track candidates = " << trackHandle->size();

  // Ref to Candidate object to be recorded in filter object
  edm::Ref<reco::RecoChargedCandidateCollection> trackRef;

  unsigned int n(0);
  for (unsigned int i=0; i<trackHandle->size(); i++)
    {
      trackRef = edm::Ref<reco::RecoChargedCandidateCollection>(trackHandle,i);
      
      if ( trackRef->pt()>minPt_ && trackRef->p()>minP_ && fabs(trackRef->eta())<maxEta_ ) {
	++n;
	filterproduct->addObject(trigger::TriggerTrack,trackRef);
      }
      
    }

  if ( edm::isDebugEnabled() ) {
    std::ostringstream stream;
    stream << "Found " << n << " track candidates with pt / p / eta" << std::endl;
    std::vector<reco::RecoChargedCandidateRef> refs;
    filterproduct->getObjects(trigger::TriggerTrack,refs);
    for ( unsigned int i=0; i<refs.size(); ++i ) {
      stream << "  " << i << " " 
	     << refs[i]->pt() << " "
	     << refs[i]->p() << " "
	     << refs[i]->eta() << std::endl;
    }
    LogDebug("HLTOniaTrackFilter") << stream.str();
  }

  iEvent.put(filterproduct);

  return n>0;

}
	  

//define this as a plug-in
DEFINE_FWK_MODULE(HLTOniaTrackFilter);
