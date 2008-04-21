#include "Workspace/EventSelectors/interface/JetEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>

JetEventSelector::JetEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset) {
  // input collection
  jetTag_ = pset.getParameter<edm::InputTag>("src");
  // jet correction
  correction_ = pat::Jet::correctionType(pset.getParameter<std::string>("correction"));
  // lower cuts on jet Et (defines also min. nr. of jets)
  minEt_ = pset.getParameter< std::vector<double> >("minEt");
  // upper cuts on jet |eta| (defines also min. nr. of jets)
  maxEta_ = pset.getParameter< std::vector<double> >("maxEta");
  // upper cuts on jet EM fraction (defines also min. nr. of jets)
  maxFem_ = pset.getParameter< std::vector<double> >("maxEMFraction");

  edm::LogInfo("JetEventSelector") << "constructed with \n"
				   << "  src = " << jetTag_ << "\n"
				   << "  min #jets = " << minEt_.size();
}

bool
JetEventSelector::select (const edm::Event& event) const
{
  //
  if ( minEt_.size()!=maxEta_.size() ||
       maxFem_.size()!=maxEta_.size() ) {
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
  // sort jets by corrected eta
  //
  std::vector<float> correctedEts;
  correctedEts.reserve(jetHandle->size());
  for ( size_t i=0; i<jetHandle->size(); ++i ) {
    const pat::Jet& jet = (*jetHandle)[i];
    float et = jet.et();
    if ( correction_ != pat::Jet::DefaultCorrection ) {
      et *= jet.correctionFactor(pat::Jet::NoCorrection);
      if ( correction_ != pat::Jet::NoCorrection )  
	et *= jet.correctionFactor(correction_);
    }
    correctedEts.push_back(et);
//     std::cout << "Applying jet correction original / nocorr / " << correction_
// 	      << " / new = " << jet.et() 
// 	      << " " << jet.correctionFactor(pat::Jet::NoCorrection)
// 	      << " " << jet.correctionFactor(correction_) 
// 	      << " " << correctedEts.back() << std::endl;
  }
//   std::cout << "Ets before sorting =";
//   for ( size_t i=0; i<jetHandle->size(); ++i )  std::cout << " " << (*jetHandle)[i].et();
//   std::cout << std::endl;
//   std::vector<size_t> sortedIndices = IndexSorter< std::vector<float> >(correctedEts,true)();
//   std::cout << "Ets after sorting =";
//   for ( size_t i=0; i<jetHandle->size(); ++i )  std::cout << " " << correctedEts[sortedIndices[i]];
//   std::cout << std::endl;
  //
  // check cuts (assume that jets are sorted by Et)
  //
  for ( unsigned int i=0; i<minEt_.size(); ++i ) {
    if ( (*jetHandle)[i].et()<minEt_[i] ||
	 fabs((*jetHandle)[i].eta())>maxEta_[i] ||
	 (*jetHandle)[i].emEnergyFraction()>maxFem_[i] ) {
      LogDebug("JetEventSelector") << "failed at jet " << (i+1);
      return false;
    }
  }
  LogDebug("JetEventSelector") << "all jets passed";
  return true;
}
