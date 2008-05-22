#include "SusyAnalysis/EventSelector/interface/EventSelectorOR.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//________________________________________________________________________________________
EventSelectorOR::EventSelectorOR (const edm::ParameterSet& pset) :
  CombinedEventSelector(pset)
{
  edm::LogInfo("EventSelectorOR") << "constructed with " << sequence_.size() << " components";
}

//________________________________________________________________________________________
bool
EventSelectorOR::select (const edm::Event& event) const
{

  // logical OR of all results

  bool result(false); // Default: none passed
  const std::vector<const SusyEventSelector*>& selectors = sequence_.selectors();
  for ( unsigned int i=0; i<selectors.size(); ++i ) 
    if ( selectors[i]->select(event) ) {
      LogDebug("EventSelectorOR") << "Event accepted by " << selectors[i]->name();
      result = true;  // passed!
    }

  if ( !result ) LogDebug("EventSelectorOR") << "Event rejected";

  return result;

}
