#include "Workspace/EventSelectors/interface/EventSelectorOR.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

EventSelectorOR::EventSelectorOR (const edm::ParameterSet& pset) :
  SusyEventSelector(pset), sequence_(pset) {}

bool
EventSelectorOR::select (const edm::Event& event) const
{
  //
  // logical OR of all results
  //
  const std::vector<const SusyEventSelector*>& selectors = sequence_.selectors();
  for ( unsigned int i=0; i<selectors.size(); ++i ) 
    if ( selectors[i]->select(event) )  return true;
  //
  return false;
}
