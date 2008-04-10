#include "Workspace/EventSelectors/interface/EventSelectorAND.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

EventSelectorAND::EventSelectorAND (const edm::ParameterSet& pset) :
  SusyEventSelector(pset), sequence_(pset) {}

bool
EventSelectorAND::select (const edm::Event& event) const
{
  //
  // logical AND of all results
  //
  const std::vector<const SusyEventSelector*>& selectors = sequence_.selectors();
  for ( unsigned int i=0; i<selectors.size(); ++i ) 
    if ( !selectors[i]->select(event) )  return false;
  //
  return true;
}
