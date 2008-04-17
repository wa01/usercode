#include "Workspace/EventSelectors/interface/SelectorSequence.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "Workspace/EventSelectors/interface/EventSelectorFactory.h"

#include <iostream>

SelectorSequence::SelectorSequence (const edm::ParameterSet& iConfig) 
{
  // retrieve parameter sets defining the selector modules
  std::vector<edm::ParameterSet> selectorConfigs =
    iConfig.getParameter< std::vector<edm::ParameterSet> >("Selectors");
  //
  // instantiate selectors according to parameter
  //
  for ( std::vector<edm::ParameterSet>::const_iterator i=selectorConfigs.begin();
	i!=selectorConfigs.end(); ++i ) {
    std::string selectorType = i->getParameter<std::string>("selector");
    const SusyEventSelector* selector = EventSelectorFactory::get()->create(selectorType,*i);
    selectors_.push_back(selector);
  }
  
}


SelectorSequence::~SelectorSequence()
{
  //
  // delete selectors
  // 
  for ( std::vector<const SusyEventSelector*>::const_iterator i=selectors_.begin();
	i!=selectors_.end(); ++i )  delete *i;
}

const std::vector<std::string>&
SelectorSequence::selectorNames () const
{
  if ( selectorNames_.empty() ) {
    selectorNames_.reserve(size());
    for ( size_t i=0; i<size(); ++i )
      selectorNames_.push_back(selectors_[i]->name());
  }
  return selectorNames_;
  
}

std::vector<bool>
SelectorSequence::decisions (const edm::Event& iEvent) const
{
  std::vector<bool> selection(size(),false);
  for ( size_t i=0; i<selectors_.size(); ++i )
    selection[i] = selectors_[i]->select(iEvent);
  return selection;
}

size_t
SelectorSequence::selectorIndex (const std::string& name) const
{
  const std::vector<std::string>& names = selectorNames();
  std::vector<std::string>::const_iterator idx = 
    find(names.begin(),names.end(),name);
  if ( idx==names.end() ) 
    edm::LogError("SelectorSequence") << "undefined selector " << name;
  return idx-names.begin();
}

std::string
SelectorSequence::selectorName (size_t index) const
{
  if ( index<size() ) {
    return selectorNames()[index];
  }
  else {
    edm::LogError("SelectorSequence") << "selector index outside range: " << index;
    return std::string();
  }
}
