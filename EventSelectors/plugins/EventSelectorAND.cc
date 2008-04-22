#include "Workspace/EventSelectors/interface/EventSelectorAND.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

EventSelectorAND::EventSelectorAND (const edm::ParameterSet& pset) :
  SusyEventSelector(pset), 
  sequence_(pset.getParameter< std::vector<std::string> >("components"),
	    pset.getParameter<edm::ParameterSet>("_AllFilters")) {
  edm::LogInfo("EventSelectorAND") << "constructed with " << sequence_.size() << " components";
}

bool
EventSelectorAND::select (const edm::Event& event) const
{
  //
  // logical AND of all results
  //
  const std::vector<const SusyEventSelector*>& selectors = sequence_.selectors();
  for ( unsigned int i=0; i<selectors.size(); ++i ) 
    if ( !selectors[i]->select(event) ) {
      LogDebug("EventSelectorAND") << "Event rejected by " << selectors[i]->name();
      return false;
    }
  //
  LogDebug("EventSelectorAND") << "Event accepted";
  return true;
}
