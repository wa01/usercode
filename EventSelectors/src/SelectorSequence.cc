#include "Workspace/EventSelectors/interface/SelectorSequence.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "Workspace/EventSelectors/interface/EventSelectorFactory.h"

#include <iostream>

SelectorSequence::SelectorSequence (const edm::ParameterSet& iConfig) 
{
  edm::LogInfo("SelectorSequence") << "created by PSet";
  // retrieve parameter sets defining the selector modules
  edm::ParameterSet filters =
    iConfig.getParameter<edm::ParameterSet>("filters");
  std::vector<std::string> sequence = 
    iConfig.getParameter< std::vector<std::string> >("selectionSequence");
  //
  createSelectors(sequence,filters);
}

SelectorSequence::SelectorSequence (const std::vector<std::string>& sequence,
				    const edm::ParameterSet& filters)
{
  edm::LogInfo("SelectorSequence") << "created by strings";
  createSelectors(sequence,filters);
}

void
SelectorSequence::createSelectors (const std::vector<std::string>& sequence,
				   const edm::ParameterSet& filters)
{
  //
  // go through sequence and instantiate selectors
  //
  for ( std::vector<std::string>::const_iterator i=sequence.begin();
	i!=sequence.end(); ++i ) {
    // retrieve filter definition
    edm::ParameterSet filterPSet = filters.getParameter<edm::ParameterSet>(*i);
    // get selector type
    std::string selectorType = filterPSet.getParameter<std::string>("selector");
    // add name
    filterPSet.addUntrackedParameter<std::string>("name",*i);
    edm::LogInfo("SelectorSequence") << "creating selector of type " << selectorType
				     << " with name " << *i;
    // add full list of filters (for combined selectors)
    filterPSet.addParameter<edm::ParameterSet>("_AllFilters",filters);
    // create selector
    const SusyEventSelector* selector = EventSelectorFactory::get()->create(selectorType,filterPSet);
    selectors_.push_back(selector);
  }
  // prepare cached decision vector
  currentDecisions_.resize(selectors_.size(),false);
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
  if ( newEvent(iEvent) ) {
    // calculate results for all selectors and update cache
    for ( size_t i=0; i<selectors_.size(); ++i )
      currentDecisions_[i] = selectors_[i]->select(iEvent);
    currentEventId_ = iEvent.id();
  }

  return currentDecisions_;
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

bool 
SelectorSequence::globalDecision (const edm::Event& event) const
{
  // make sure that cache is updated
  decisions(event);
  //
  for ( size_t i=0; i<size(); ++i ) {
    if ( !currentDecisions_[i] )  return false;
  }
  return true;
}

bool 
SelectorSequence::decision (const edm::Event& event, size_t index) const
{
  if ( index>=size() ) {
    edm::LogError("SelectorSequence") << "selector index outside range: " << index;
    return false;
  }
  return decisions(event)[index];
}

bool 
SelectorSequence::decision (const edm::Event& event, const std::string& name) const
{
  return decisions(event)[selectorIndex(name)];
}

bool 
SelectorSequence::complementaryDecision (const edm::Event& event, size_t index) const
{
  // make sure that cache is updated
  decisions(event);
  //
  for ( size_t i=0; i<size(); ++i ) {
    // ignore decision of selector "index"
    if ( i==index ) continue;
    // implement AND of other selectors
    if ( !currentDecisions_[i] )  return false;
  }
  return true;
}

bool 
SelectorSequence::complementaryDecision (const edm::Event& event, 
					 const std::string& name) const
{
  return complementaryDecision(event,selectorIndex(name));
}
