#include "SusyAnalysis/EventSelector/interface/CombinedEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

//________________________________________________________________________________________
CombinedEventSelector::CombinedEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset), 
  sequence_(pset.getParameter< std::vector<std::string> >("components"),
            pset.getParameter<edm::ParameterSet>("_AllFilters"))
{

}

//________________________________________________________________________________________
double 
CombinedEventSelector::value (const std::string& name) const {

  // Decode variable name (of the form "selector:name") and get it from the sequence
  
  size_t index = name.find(":");
  if ( index == std::string::npos ) {
    edm::LogError("CombinedEventSelector") << "unable to get selector name from \"" << name << "\""
                                      << "\n(looking for \"selector:variable\")";
    return susy::DEFAULT_VALUE;
  }
  
  return sequence_.value( name.substr(0,index), name.substr(index+1) );

}
