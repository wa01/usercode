#include "SusyAnalysis/EventSelector/interface/EMFractionSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <vector>

//__________________________________________________________________________________________________
EMFractionSelector::EMFractionSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset),
  jetTag_( pset.getParameter<edm::InputTag>("src") ),
  minFraction_( pset.getParameter<double>("minFraction") )
{

  // Store computed fraction
  defineVariable("emFraction");

}

//__________________________________________________________________________________________________
bool
EMFractionSelector::select (const edm::Event& event) const
{

  // reset cached variables
  resetVariables();
  // Get the jets
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  event.getByLabel(jetTag_, jetHandle);
  if ( !jetHandle.isValid() ) {
    edm::LogWarning("EMFractionSelector") << "No Jet results for InputTag " << jetTag_;
    return false;
  }

  // Calculate sum of fractions
  float etSum  = 0.;
  float etFrac = 0.;

  for ( std::vector< pat::Jet >::const_iterator iJet = jetHandle->begin();
        iJet != jetHandle->end(); ++iJet )
    {
      etFrac += iJet->et()*(iJet->emEnergyFraction());
      etSum  += iJet->et();
    }
  
  float emFraction = etFrac/etSum;

  setVariable("emFraction",emFraction); // Cache variable

  // Return selection
  return emFraction > minFraction_;

}


//__________________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"
DEFINE_EDM_PLUGIN(EventSelectorFactory, EMFractionSelector, "EMFractionSelector");
