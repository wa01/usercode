#include "Workspace/EventSelectors/interface/EventSelectorOR.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

EventSelectorOR::EventSelectorOR (const edm::ParameterSet& pset) :
  SusyEventSelector(pset),
  sequence_(pset.getParameter< std::vector<std::string> >("components"),
	    pset.getParameter<edm::ParameterSet>("_AllFilters")) {
  edm::LogInfo("EventSelectorOR") << "constructed with " << sequence_.size() << " components";
}

bool
EventSelectorOR::select (const edm::Event& event) const
{
  //
  // logical OR of all results
  //
  const std::vector<const SusyEventSelector*>& selectors = sequence_.selectors();
  for ( unsigned int i=0; i<selectors.size(); ++i ) 
    if ( selectors[i]->select(event) ) {
      LogDebug("EventSelectorOR") << "Event accepted by " << selectors[i]->name();
      return true;
    }
  //
  LogDebug("EventSelectorOR") << "Event rejected";
  return false;
}
