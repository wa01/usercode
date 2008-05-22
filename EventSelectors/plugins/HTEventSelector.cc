#include "SusyAnalysis/EventSelector/interface/HTEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

//__________________________________________________________________________________________________
HTEventSelector::HTEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset),
  jetTag_( pset.getParameter<edm::InputTag>("jet") ),
  minHT_ ( pset.getParameter<double>("minHT") )
{ 

  // Store computed HT
  defineVariable("HT");

}

//__________________________________________________________________________________________________
bool
HTEventSelector::select (const edm::Event& event) const
{
  // reset cached variables
  resetVariables();
  // Get the jets
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  event.getByLabel(jetTag_, jetHandle);
  if ( !jetHandle.isValid() ) {
    edm::LogWarning("HTEventSelector") << "No Jet results for InputTag " << jetTag_;
    return false;
  }

  // Sum over jet Ets
  float myHT = 0.0;
  std::vector<pat::Jet>::const_iterator iJet = jetHandle->begin();
  while ( iJet != jetHandle->end() ) {
    myHT += iJet->et();
    ++iJet;
  }

  setVariable("HT",myHT);

  return myHT>minHT_;

}

//__________________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"
DEFINE_EDM_PLUGIN(EventSelectorFactory, HTEventSelector, "HTEventSelector");
