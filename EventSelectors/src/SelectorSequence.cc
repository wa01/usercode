#include "Workspace/EventSelectors/interface/SelectorSequence.h"

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


std::vector<bool>
SelectorSequence::decisions (const edm::Event& iEvent) const
{
  std::vector<bool> selection(selectors_.size(),false);
  for ( unsigned int i=0; i<selectors_.size(); ++i )
    selection[i] = selectors_[i]->select(iEvent);
  return selection;
}
