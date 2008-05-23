#include "SusyAnalysis/EventSelector/interface/JetEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <sstream>

JetEventSelector::JetEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset) {
  // input collection
  jetTag_ = pset.getParameter<edm::InputTag>("jetTag");
  // jet correction
  correction_ = pat::Jet::correctionType(pset.getParameter<std::string>("correction"));
  // lower cuts on jet Et (defines also min. nr. of jets)
  minEt_ = pset.getParameter< std::vector<double> >("minEt");
  // upper cuts on jet |eta| (defines also min. nr. of jets)
  maxEta_ = pset.getParameter< std::vector<double> >("maxEta");
  // upper cuts on jet EM fraction (defines also min. nr. of jets)
  maxFem_ = pset.getParameter< std::vector<double> >("maxEMFraction");

  /// definition of variables to be cached
  defineVariable("NumberOfJets");
  for ( size_t i=0; i<minEt_.size(); ++i ) {
    std::ostringstream strEt;
    strEt << "Jet" << i << "Et";
    defineVariable(strEt.str());
    std::ostringstream strEta;
    strEta << "Jet" << i << "Eta";
    defineVariable(strEta.str());
    std::ostringstream strFem;
    strFem << "Jet" << i << "EMfraction";
    defineVariable(strFem.str());
  }

  edm::LogInfo("JetEventSelector") << "constructed with \n"
				   << "  jetTag    = " << jetTag_ << "\n"
				   << "  min #jets = " << minEt_.size();
}

bool
JetEventSelector::select (const edm::Event& event) const
{
  // reset cached variables
  resetVariables();
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
  setVariable(0,jetHandle->size());
  if ( jetHandle->size()<minEt_.size() )  return false;
  //
  // sort jets by corrected Et
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
  }
  //
  // check cuts (assume that jets are sorted by Et)
  //
  bool result(true);
  for ( unsigned int i=0; i<minEt_.size(); ++i ) {
    if ( (*jetHandle)[i].et()<minEt_[i] ||
	 fabs((*jetHandle)[i].eta())>maxEta_[i] ||
	 (*jetHandle)[i].emEnergyFraction()>maxFem_[i] ) {
      LogTrace("JetEventSelector") << "JetEventSelector: failed at jet " << (i+1);
      result = false;
    }
    setVariable(3*i+1,(*jetHandle)[i].et());
    setVariable(3*i+2,(*jetHandle)[i].eta());
    setVariable(3*i+3,(*jetHandle)[i].emEnergyFraction());
  }
  LogTrace("JetEventSelector") << "JetEventSelector: all jets passed";
  return result;
}

//________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"
DEFINE_EDM_PLUGIN(EventSelectorFactory, JetEventSelector, "JetEventSelector");
