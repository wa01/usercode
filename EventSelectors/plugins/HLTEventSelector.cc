#include "Workspace/EventSelectors/interface/HLTEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

HLTEventSelector::HLTEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset) {
  // input collection
  triggerResults_ = pset.getParameter<edm::InputTag>("triggerResults");
  // trigger path names
  pathNames_ = pset.getParameter< std::vector<std::string> >("pathNames");
}

bool
HLTEventSelector::select (const edm::Event& event) const
{
  //
  // get the trigger results and check validity
  //
  edm::Handle<edm::TriggerResults> hltHandle;
  event.getByLabel(triggerResults_, hltHandle);
  if ( !hltHandle.isValid() ) {
    edm::LogWarning("HLTEventSelector") << "No trigger results for InputTag " << triggerResults_;
    return false;
  }
  //
  // get results
  //
  edm::TriggerNames trgNames;
  trgNames.init(*hltHandle);
  unsigned int trgSize = trgNames.size();
//   static int first(true);
//   if ( first ) {
//     first = false;
//     std::cout << "Trigger menu" << std::endl;
//     for ( unsigned int i=0; i<trgSize; ++i ) {
//       std::cout << trgNames.triggerName(i) << std::endl;
//     }
//   }
//
// example for OR of all specified triggers
//
  for ( std::vector<std::string>::const_iterator i=pathNames_.begin();
	i!=pathNames_.end(); ++i ) {
    // get index
    unsigned int index = trgNames.triggerIndex(*i);
    if ( index==trgSize ) {
      edm::LogWarning("HLTEventSelector") << "Unknown trigger name " << *i;
//       return false;
      continue;
    }
//     if ( !hltHandle->accept(index) )  return false;
    if ( hltHandle->accept(index) )  return true;
  }
//   return true;
  return false;
}
