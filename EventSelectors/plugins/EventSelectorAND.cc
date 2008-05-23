#include "SusyAnalysis/EventSelector/interface/EventSelectorAND.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

//________________________________________________________________________________________
EventSelectorAND::EventSelectorAND (const edm::ParameterSet& pset) :
  CombinedEventSelector(pset)
{
  edm::LogInfo("EventSelectorAND") << "constructed with " << sequence_.size() << " components";
}

//________________________________________________________________________________________
bool
EventSelectorAND::select (const edm::Event& event) const
{
  // logical AND of all results

  bool result(true);  // Default: all passed
  const std::vector<const SusyEventSelector*>& selectors = sequence_.selectors();
  for ( unsigned int i=0; i<selectors.size(); ++i ) 
    if ( !selectors[i]->select(event) ) {
      LogDebug("EventSelectorAND") << "Event rejected by " << selectors[i]->name();
      result = false; // not passed!
    }

  if ( result ) LogDebug("EventSelectorAND") << "Event accepted";

  return result;

}

//________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"
DEFINE_EDM_PLUGIN(EventSelectorFactory, EventSelectorAND, "EventSelectorAND");
