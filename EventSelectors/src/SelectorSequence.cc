#include "SusyAnalysis/EventSelector/interface/SelectorSequence.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"

#include <iostream>

//________________________________________________________________________________________
SelectorSequence::SelectorSequence (const edm::ParameterSet& iConfig) 
{
  edm::LogInfo("SelectorSequence") << "created by PSet";
  
  // retrieve parameter sets defining the selector modules
  edm::ParameterSet selectors =
    iConfig.getParameter<edm::ParameterSet>("selectors");
  std::vector<std::string> sequence = 
    iConfig.getParameter< std::vector<std::string> >("selectionSequence");
  //
  createSelectors(sequence,selectors);
}

//________________________________________________________________________________________
SelectorSequence::SelectorSequence (const std::vector<std::string>& sequence,
				    const edm::ParameterSet& selectors)
{
  edm::LogInfo("SelectorSequence") << "created by strings";
  createSelectors(sequence,selectors);
}

//________________________________________________________________________________________
void
SelectorSequence::createSelectors (const std::vector<std::string>& sequence,
				   const edm::ParameterSet& selectors)
{
  //
  // go through sequence and instantiate selectors
  //
  for ( std::vector<std::string>::const_iterator i=sequence.begin();
	i!=sequence.end(); ++i ) {
    // retrieve selector definition
    edm::ParameterSet selectorPSet = selectors.getParameter<edm::ParameterSet>(*i);
    // get selector type
    std::string selectorType = selectorPSet.getParameter<std::string>("selector");
    // add name
    selectorPSet.addUntrackedParameter<std::string>("name",*i);
    edm::LogVerbatim("SelectorSequence") << "creating selector of type " << selectorType
                                         << " with name " << *i;
    // add full list of selectors (for combined selectors)
    selectorPSet.addParameter<edm::ParameterSet>("_AllFilters",selectors);
    // create selector
    const SusyEventSelector* selector = EventSelectorFactory::get()->create(selectorType,selectorPSet);
    selectors_.push_back(selector);
  }
  // prepare cached decision vector
  currentDecisions_.resize(selectors_.size(),false);
}

//________________________________________________________________________________________
SelectorSequence::~SelectorSequence()
{
  //
  // delete selectors
  // 
  for ( std::vector<const SusyEventSelector*>::const_iterator i=selectors_.begin();
	i!=selectors_.end(); ++i )  delete *i;
}

//________________________________________________________________________________________
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

//________________________________________________________________________________________
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

//________________________________________________________________________________________
size_t
SelectorSequence::selectorIndex (const std::string& selectorName) const
{
  const std::vector<std::string>& names = selectorNames();
  std::vector<std::string>::const_iterator idx = 
    find(names.begin(),names.end(),selectorName);
  if ( idx==names.end() ) 
    edm::LogError("SelectorSequence") << "undefined selector " 
				      << selectorName;
  return idx-names.begin();
}

//________________________________________________________________________________________
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

//________________________________________________________________________________________
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

//________________________________________________________________________________________
bool 
SelectorSequence::decision (const edm::Event& event, size_t index) const
{
  if ( index>=size() ) {
    edm::LogError("SelectorSequence") << "selector index outside range: " << index;
    return false;
  }
  return decisions(event)[index];
}

//________________________________________________________________________________________
bool 
SelectorSequence::decision (const edm::Event& event, 
			    const std::string& selectorName) const
{
  return decisions(event)[selectorIndex(selectorName)];
}

//________________________________________________________________________________________
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

//________________________________________________________________________________________
bool 
SelectorSequence::complementaryDecision (const edm::Event& event, 
					 const std::string& selectorName) const
{
  return complementaryDecision(event,selectorIndex(selectorName));
}

//________________________________________________________________________________________
bool 
SelectorSequence::cumulativeDecision (const edm::Event& event, size_t index) const
{
  // check index
  if ( index >= size() ) {
    edm::LogError("SelectorSequence") << "selector index outside range: " << index;
    return false;
  }
  // make sure that cache is updated
  decisions(event);
  //
  for ( size_t i=0; i<=index; ++i ) {
    if ( !currentDecisions_[i] )  return false;
  }
  return true;
}

//________________________________________________________________________________________
bool 
SelectorSequence::cumulativeDecision (const edm::Event& event, 
				      const std::string& selectorName) const
{
  return cumulativeDecision(event,selectorIndex(selectorName));
}

//________________________________________________________________________________________
size_t
SelectorSequence::numberOfVariables () const
{
  size_t result(0);
  for ( std::vector<const SusyEventSelector*>::const_iterator i=selectors_.begin();
	i!=selectors_.end(); ++i ) result += (**i).numberOfVariables();
  return result;
}

//________________________________________________________________________________________
std::vector<std::string>
SelectorSequence::variableNames () const 
{

  std::vector<std::string> result;
  result.reserve(numberOfVariables());

  std::string sName;
  std::vector<std::string> vNames;
  for ( size_t i=0; i<size(); ++i ) {
    sName = selectorNames()[i] + ":";
    vNames = selectors_[i]->variableNames();
    for ( std::vector<std::string>::const_iterator j=vNames.begin();
	  j!=vNames.end(); ++j ) 
      result.push_back( sName+(*j) );
  }
  
  return result;
}

//________________________________________________________________________________________
std::vector<double>
SelectorSequence::values () const 
{
  std::vector<double> result;
  result.reserve(numberOfVariables());

  std::vector<double> values;
  for ( std::vector<const SusyEventSelector*>::const_iterator i=selectors_.begin();
	i!=selectors_.end(); ++i ) {
    values = (**i).values();
    result.insert(result.end(),values.begin(),values.end());
  }
  
  return result;
}

//________________________________________________________________________________________
double
SelectorSequence::value (const std::string& selectorName,
			 const std::string& variableName) const
{
  size_t iSel = selectorIndex(selectorName);
  if ( iSel==selectors_.size() ) {
    edm::LogError("SelectorSequence") << "undefined selector " 
				      << selectorName;
    return susy::DEFAULT_VALUE;
  }
  return selectors_[iSel]->value(variableName);
}
