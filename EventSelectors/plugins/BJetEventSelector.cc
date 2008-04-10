#include "Workspace/EventSelectors/interface/BJetEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <vector>

BJetEventSelector::BJetEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset) {
  // input collection
  jetTag_ = pset.getParameter<edm::InputTag>("src");
  // name of tagger
  tagLabel_ = pset.getParameter<std::string>("tagLabel");
  // lower cuts on discriminator (defines also min. nr. of jets)
  minTag_ = pset.getParameter< std::vector<double> >("minTag");
  
  edm::LogInfo("BJetEventSelector") << "constructed with \n"
				    << "  src = " << jetTag_ << "\n"
				    << "  tagger = " << tagLabel_ << "\n"
				    << "  min #jets = " << minTag_.size();
}

bool
BJetEventSelector::select (const edm::Event& event) const
{
  // get the jets
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  event.getByLabel(jetTag_, jetHandle);
  if ( !jetHandle.isValid() ) {
    edm::LogWarning("BJetEventSelector") << "No Jet results for InputTag " << jetTag_;
    return false;
  }
  //
  // check number of jets
  //
  if ( jetHandle->size()<minTag_.size() )  return false;
  //
  // sort discriminator value
  //
  std::vector<float> discriminators;
  discriminators.reserve(jetHandle->size());
  for ( unsigned int i=0; i<jetHandle->size(); ++i ) 
    discriminators.push_back((*jetHandle)[i].bDiscriminator(tagLabel_));
  std::sort(discriminators.begin(),discriminators.end());
  //
  // apply cuts
  //
  for ( unsigned int i=0; i<minTag_.size(); ++i ) {
    if ( discriminators[i]<minTag_[i] ) {
      LogDebug("BJetEventSelector") << "failed at jet " << (i+1);
      return false;
    }
  }
  LogDebug("BJetEventSelector") << "all jets passed";
  return true;
}
