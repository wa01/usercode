#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
 

#include "SusyAnalysis/EventSelector/interface/MeffSelector.h"


//________________________________________________________________________________________
MeffSelector::MeffSelector (const edm::ParameterSet& pset ) :
  SusyEventSelector(pset),
  jetTag_(  pset.getParameter<edm::InputTag>("jetTag") ),
  metTag_(  pset.getParameter<edm::InputTag>("metTag") ),
  minMeff_( pset.getParameter<double>("minMeff")       )
{

  // Define variables to store
  defineVariable("Meff");

}

//________________________________________________________________________________________
bool
MeffSelector::select (const edm::Event& event) const
{

  // Reset cached variables
  resetVariables();

  // Get the MET result
  edm::Handle< std::vector<pat::MET> > metHandle;
  event.getByLabel(metTag_, metHandle);
  if ( !metHandle.isValid() ) {
    edm::LogWarning("MeffEventSelector") << "No Met results for InputTag " << metTag_;
    return false;
  }
  if ( metHandle->size()!=1 ) {
    edm::LogWarning("MeffEventSelector") << "MET collection size is "
                                         << metHandle->size() << " instead of 1";    
    return false;
  }

  // Get the jets
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  event.getByLabel(jetTag_, jetHandle);
  if ( !jetHandle.isValid() ) {
    edm::LogWarning("MetJetEventSelector") << "No Jet results for InputTag " << jetTag_;
    return false;
  }  

  // Compute Meff
  float meff = (*metHandle)[0].et();
  for( std::vector< pat::Jet >::const_iterator iJet = jetHandle->begin();
       iJet != jetHandle->end(); ++iJet )
    {
      meff += meff + iJet->et();
    }

  setVariable("Meff",meff);
  
  // Select
  return meff > minMeff_;

}

//__________________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"
DEFINE_EDM_PLUGIN(EventSelectorFactory, MeffSelector, "MeffSelector");

