#include "SusyAnalysis/EventSelector/interface/HLTEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"

HLTEventSelector::HLTEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset) {
  // input collection
  triggerResults_ = pset.getParameter<edm::InputTag>("triggerResults");
  // trigger path names
  pathNames_ = pset.getParameter< std::vector<std::string> >("pathNames");

  // Add all path names as variables
  for ( std::vector<std::string>::const_iterator i=pathNames_.begin();
	i!=pathNames_.end(); ++i )
    {
      defineVariable( *i );
    }

  edm::LogInfo("HLTEventSelector") << "constructed with \n"
				   << "  src = " << triggerResults_ << "\n"
				   << "  #pathnames = " << pathNames_.size();
}

bool
HLTEventSelector::select (const edm::Event& event) const
{
  // Reset cached variables
  resetVariables();

  // Get the trigger results and check validity
  edm::Handle<edm::TriggerResults> hltHandle;
  event.getByLabel(triggerResults_, hltHandle);
  if ( !hltHandle.isValid() ) {
    edm::LogWarning("HLTEventSelector") << "No trigger results for InputTag " << triggerResults_;
    return false;
  }

  // Get results
  edm::TriggerNames trgNames;
  trgNames.init(*hltHandle);
  unsigned int trgSize = trgNames.size();
  bool result(false);

  // Example for OR of all specified triggers
  for ( std::vector<std::string>::const_iterator i=pathNames_.begin();
	i!=pathNames_.end(); ++i ) {
    // Get index
    unsigned int index = trgNames.triggerIndex(*i);
    if ( index==trgSize ) {
      edm::LogWarning("HLTEventSelector") << "Unknown trigger name " << *i;
      continue;
    }
    if ( hltHandle->accept(index) ) {
      LogDebug("HLTEventSelector") << "Event selected by " << *i;
      result = true;
      setVariable( *i, true );
    } else {
      setVariable( *i, false );
    }
  }

  if (!result) LogDebug("HLTEventSelector") << "Event rejected";

  return result;
}
