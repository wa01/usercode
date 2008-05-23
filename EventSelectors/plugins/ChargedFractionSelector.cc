#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "SusyAnalysis/EventSelector/interface/ChargedFractionSelector.h"


//________________________________________________________________________________________
ChargedFractionSelector::ChargedFractionSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset),
  jetTag_( pset.getParameter<edm::InputTag>("jetTag") ),
  minFraction_( pset.getParameter<double>("minFraction") ),
  maxEta_     ( pset.getParameter<double>("maxEta")      ),
  minPt_      ( pset.getParameter<double>("minPt")       ),
  minTracks_  ( pset.getParameter<int>("minTracks")      )
{

  // Store computed fraction
  defineVariable("chFraction");

}

//________________________________________________________________________________________
bool
ChargedFractionSelector::select (const edm::Event& event) const
{

  // reset cached variables
  resetVariables();

  // Get the jets
  edm::Handle< std::vector<pat::Jet> > jets;
  event.getByLabel(jetTag_, jets );
  if ( !jets.isValid() ) {
    edm::LogWarning("ChargedFractionSelector") << "No Jet results for InputTag " << jetTag_;
    return false;
  }

  // Calculate sum of fractions
  float etFrac = 0.;
  // Loop on jets
  for ( std::vector< pat::Jet >::const_iterator iJet = jets->begin();
        iJet != jets->end(); ++iJet )
    {
      const reco::TrackRefVector& jetTracks = iJet->associatedTracks();
      // Acceptance cuts
      if (    iJet->pt()        > minPt_
           && fabs(iJet->eta()) < maxEta_ 
           && jetTracks.size()  >= minTracks_ ) {
        // Loop on associated tracks
        for ( reco::TrackRefVector::const_iterator iTrack = jetTracks.begin();
              iTrack!= jetTracks.end(); ++iTrack ) {
          etFrac += (*iTrack)->pt()/iJet->et();
        }
      }
    }
  
  // Average charge fraction: if no jets, keep default (extreme) value
  double chFraction = 0.;
  if ( jets->size()>0 ) {
    chFraction = etFrac/static_cast<float>(jets->size()); // Make sure we convert size to float
    setVariable("chFraction",chFraction); // Cache variable
  }

  // Return selection
  return chFraction > minFraction_;

}


//________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"
DEFINE_EDM_PLUGIN(EventSelectorFactory, ChargedFractionSelector, "ChargedFractionSelector");
