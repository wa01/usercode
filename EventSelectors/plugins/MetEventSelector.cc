#include "SusyAnalysis/EventSelector/interface/MetEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include <vector>

MetEventSelector::MetEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset) {
  // input collection
  metTag_ = pset.getParameter<edm::InputTag>("metTag");
  // lower cut on MET
  minMet_ = pset.getParameter<double>("minMET");

  // definition of the variable to be cached
  defineVariable("MET");

  edm::LogInfo("MetEventSelector") << "constructed with \n"
				   << "  metTag = " << metTag_ << "\n"
				   << "  minMET = " << minMet_;
}

bool
MetEventSelector::select (const edm::Event& event) const
{
  // reset cached variables
  resetVariables();
  //
  // get the MET result
  //
  edm::Handle< std::vector<pat::MET> > metHandle;
  event.getByLabel(metTag_, metHandle);
  if ( !metHandle.isValid() ) {
    edm::LogWarning("MetEventSelector") << "No Met results for InputTag " << metTag_;
    return false;
  }
  //
  // sanity check on collection
  //
  if ( metHandle->size()!=1 ) {
    edm::LogWarning("MetEventSelector") << "MET collection size is " 
					<< metHandle->size() << " instead of 1";
    return false;
  }
  //
  // apply cut
  //
  setVariable(0,metHandle->front().et());
  LogTrace("MetEventSelector") << "result = " << (metHandle->front().et()>minMet_);
  return metHandle->front().et()>minMet_;
}

//________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"
DEFINE_EDM_PLUGIN(EventSelectorFactory, MetEventSelector, "MetEventSelector");

