//
// Wrapper to use a SusyEventSelector as EDFilter
//
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "Workspace/EventSelectors/interface/EventSelectorFactory.h"

//
// class declaration
//

class EventSelectorFilter : public edm::EDFilter {
public:
  EventSelectorFilter () : selector_(0) {}
  explicit EventSelectorFilter (const edm::ParameterSet& pset);
  ~EventSelectorFilter () {delete selector_;}

private:
  virtual void beginJob (const edm::EventSetup&) {}
  virtual bool filter (edm::Event& iEvent, const edm::EventSetup& iSetup) {
    return selector_ ? selector_->select(iEvent) : false;
  }
  virtual void endJob() {}
      
private:
  SusyEventSelector* selector_;
};

EventSelectorFilter::EventSelectorFilter (const edm::ParameterSet& pset)
{
  // definition of the selector(s)
  edm::ParameterSet filters = pset.getParameter<edm::ParameterSet>("filters");
  // name of the selector
  std::string name = pset.getParameter<std::string>("selector");

  // retrieve configuration of the selector
  edm::ParameterSet filterPSet = filters.getParameter<edm::ParameterSet>(name);
  // get selector type
  std::string selectorType = filterPSet.getParameter<std::string>("selector");
  // add name
  filterPSet.addUntrackedParameter<std::string>("name",name);
  edm::LogInfo("EventSelectorFilter") << "creating selector of type " << selectorType
				      << " with name " << name;
  // add full list of filters (for combined selectors)
  filterPSet.addParameter<edm::ParameterSet>("_AllFilters",filters);
  // create selector
  selector_ = EventSelectorFactory::get()->create(selectorType,filterPSet);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventSelectorFilter);
