#include "Workspace/EventSelectors/interface/JetEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <vector>

JetEventSelector::JetEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset) {
  // input collection
  jetTag_ = pset.getParameter<edm::InputTag>("src");
  // lower cuts on jet Et (defines also min. nr. of jets)
  minEt_ = pset.getParameter< std::vector<double> >("minEt");
  // upper cuts on jet |eta| (defines also min. nr. of jets)
  maxEta_ = pset.getParameter< std::vector<double> >("maxEta");
}

bool
JetEventSelector::select (const edm::Event& event) const
{
  //
  if ( minEt_.size() != maxEta_.size() ) {
    edm::LogError("JetEventSelector") << "Inconsistent length of vector of cut values";
    return false;
  }
  // get the jets
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  event.getByLabel(jetTag_, jetHandle);
  if ( !jetHandle.isValid() ) {
    edm::LogWarning("JetEventSelector") << "No Jet results for InputTag " << jetTag_;
    return false;
  }
  //
  // check number of jets
  //
  if ( jetHandle->size()<minEt_.size() )  return false;
//   std::cout << "Jet Et =";
//   for ( unsigned int i=0; i<jetHandle->size(); ++i )  std::cout << " " << (*jetHandle)[i].et();
//   std::cout << std::endl;
  //
  // check cuts (assume that jets are sorted by Et)
  //
  for ( unsigned int i=0; i<minEt_.size(); ++i ) {
    if ( (*jetHandle)[i].et()<minEt_[i] ||
	 fabs((*jetHandle)[i].eta())>maxEta_[i] )  return false;
  }
  return true;
}
