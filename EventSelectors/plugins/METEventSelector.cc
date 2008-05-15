#include "Workspace/EventSelectors/interface/METEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include <vector>

METEventSelector::METEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset) {
  // input collection
  metTag_ = pset.getParameter<edm::InputTag>("src");
  // lower cut on MET
  minMet_ = pset.getParameter<double>("minMET");

  // definition of the variable to be cached
  defineVariable("MET");

  edm::LogInfo("METEventSelector") << "constructed with \n"
				   << "  src = " << metTag_ << "\n"
				   << "  minMET = " << minMet_;
}

bool
METEventSelector::select (const edm::Event& event) const
{
  // reset cached variables
  resetVariables();
  //
  // get the MET result
  //
  edm::Handle< std::vector<pat::MET> > metHandle;
  event.getByLabel(metTag_, metHandle);
  if ( !metHandle.isValid() ) {
    edm::LogWarning("METEventSelector") << "No Met results for InputTag " << metTag_;
    return false;
  }
  //
  // sanity check on collection
  //
  if ( metHandle->size()!=1 ) {
    edm::LogWarning("METEventSelector") << "MET collection size is " 
					<< metHandle->size() << " instead of 1";
    return false;
  }
  //
  // apply cut
  //
  setVariable(0,metHandle->front().et());
  LogTrace("METEventSelector") << " METEventSelector result = " << (metHandle->front().et()>minMet_);
  return metHandle->front().et()>minMet_;
}
