#include "Workspace/EventSelectors/test/SusySelectorExample.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Workspace/EventSelectors/interface/SusyEventSelector.h"

#include <iostream>
#include <sstream>

SusySelectorExample::SusySelectorExample (const edm::ParameterSet& iConfig) :
  nrEventTotal_(0), nrEventSelected_(0), 
  selectors_(iConfig.getParameter<edm::ParameterSet>("selections")) {
  // should change to LogInfo ...
  std::ostringstream msg;
  msg << "Selectors are" << std::endl;
  std::vector<std::string> names = selectors_.selectorNames();
  for ( size_t i=0; i<nrOfSelectors(); ++i ) msg << "  " << names[i];
  edm::LogInfo("SusySelectorExample") << msg.str();
}

SusySelectorExample::~SusySelectorExample() {}

void
SusySelectorExample::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // retrieve the decision of each selector module
  //
  std::vector<bool> decisions = selectors_.decisions(iEvent);
  //
  // count all events / events passing all selections
  //
  ++nrEventTotal_;
  std::ostringstream dbg;
  dbg << "selector decisions " << std::endl
      << "  name, 2xdecision, 2xcompl. decision, cumul. decision" << std::endl;
  bool dec(true);
  for ( size_t i=0; i<nrOfSelectors(); ++i ) {
    std::string name = selectors_.selectorName(i);
    dec = dec && decisions[i];
    dbg << " " << name
	<< " " << selectors_.decision(iEvent,i)
	<< " " << selectors_.decision(iEvent,name)
	<< " " << selectors_.complementaryDecision(iEvent,i)
	<< " " << selectors_.complementaryDecision(iEvent,name)
	<< " " << selectors_.cumulativeDecision(iEvent,i)
	<< " " << selectors_.cumulativeDecision(iEvent,name)
	<< " " << dec << std::endl;
  }
  dbg << " global decision = " << selectors_.globalDecision(iEvent);
  LogTrace("SusySelectorExample") << "SusySelectorExample: " << dbg.str();
  if ( !selectors_.globalDecision(iEvent) )  return;
  ++nrEventSelected_;
}


void 
SusySelectorExample::beginJob(const edm::EventSetup&) {}

void 
SusySelectorExample::endJob() 
{
  std::ostringstream msg;
  msg << "Total number of events = " << nrEventTotal_
      << " ; selected = " << nrEventSelected_;
  edm::LogInfo("SusySelectorExample") << msg.str();
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SusySelectorExample);
