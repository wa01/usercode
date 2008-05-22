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
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"

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
  edm::ParameterSet selectors = pset.getParameter<edm::ParameterSet>("selectors");
  // name of the selector
  std::string name = pset.getParameter<std::string>("selector");

  // retrieve configuration of the selector
  edm::ParameterSet selectorPSet = selectors.getParameter<edm::ParameterSet>(name);
  // get selector type
  std::string selectorType = selectorPSet.getParameter<std::string>("selector");
  // add name
  selectorPSet.addUntrackedParameter<std::string>("name",name);
  edm::LogInfo("EventSelectorFilter") << "creating selector of type " << selectorType
				      << " with name " << name;
  // add full list of selectors (for combined selectors)
  selectorPSet.addParameter<edm::ParameterSet>("_AllFilters",selectors);
  // create selector
  selector_ = EventSelectorFactory::get()->create(selectorType,selectorPSet);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventSelectorFilter);
