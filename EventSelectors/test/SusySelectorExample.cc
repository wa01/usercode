#include "Workspace/EventSelectors/test/SusySelectorExample.h"

#include "Workspace/EventSelectors/interface/SusyEventSelector.h"

#include <iostream>

SusySelectorExample::SusySelectorExample (const edm::ParameterSet& iConfig) :
  nrEventTotal_(0), nrEventSelected_(0), selectors_(iConfig) {
  // should change to LogInfo ...
  std::cout << "Selectors defined are:" << std::endl;
  std::vector<std::string> names = selectors_.selectorNames();
  for ( size_t i=0; i<nrOfSelectors(); ++i ) 
    std::cout << "  " << names[i] << std::endl;
}

SusySelectorExample::~SusySelectorExample() {}

void
SusySelectorExample::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // retrieve the decision of each selector module
  //
  std::vector<bool> selection = selectors_.decisions(iEvent);
  //
  // count all events / events passing all selections
  //
  ++nrEventTotal_;
  std::vector<bool> exclusiveSelected(nrOfSelectors(),true);
  for ( size_t i=0; i<nrOfSelectors(); ++i ) {
    if ( !selection[i] )  return;
  }
  ++nrEventSelected_;
}


void 
SusySelectorExample::beginJob(const edm::EventSetup&) {}

void 
SusySelectorExample::endJob() 
{
  std::cout << "Total number of events = " << nrEventTotal_
	    << " ; selected = " << nrEventSelected_ << std::endl;
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SusySelectorExample);
