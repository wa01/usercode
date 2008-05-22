#include "SusyAnalysis/AnalysisSkeleton/interface/MetJetEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "PhysicsTools/Utilities/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "TMath.h"
#include <vector>

//________________________________________________________________________________________
MetJetEventSelector::MetJetEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset),
  metTag_( pset.getParameter<edm::InputTag>("met") ),
  jetTag_( pset.getParameter<edm::InputTag>("jet") ),
  metDphiMin_    ( pset.getParameter<double>("metDPhiMin")     ),
  dPhiJet2MetMin_( pset.getParameter<double>("dPhiJet2MetMin") ),
  rDistJetsMin_  ( pset.getParameter<double>("rDistJetsMin")   )
{

  // Define all variables we want to cache (and eventually plot...)
  defineVariable("metDphi");
  defineVariable("numberOfJets");
  defineVariable("dPhiJet2Met");
  defineVariable("R1");
  defineVariable("R2");

}

//________________________________________________________________________________________
bool
MetJetEventSelector::select (const edm::Event& event) const
{

  // Reset cached variables
  resetVariables();

  // Get the MET result
  edm::Handle< std::vector<pat::MET> > metHandle;
  event.getByLabel(metTag_, metHandle);
  if ( !metHandle.isValid() ) {
    edm::LogWarning("MetJetEventSelector") << "No Met results for InputTag " << metTag_;
    return false;
  }
  if ( metHandle->size()!=1 ) {
    edm::LogWarning("MetJetEventSelector") << "MET collection size is "
                                           << metHandle->size() << " instead of 1";    
    return false;
  }
  const std::vector<pat::MET>& mets = (*metHandle); // For simplicity...

  // Get the jets
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  event.getByLabel(jetTag_, jetHandle);
  if ( !jetHandle.isValid() ) {
    edm::LogWarning("MetJetEventSelector") << "No Jet results for InputTag " << jetTag_;
    return false;
  }
  const std::vector<pat::Jet>& jets = (*jetHandle); // For simplicity...

  //
  // Preselection: number of jets (need at least 3(?) to make sense out of these variables)
  // 
  setVariable("numberOfJets",jets.size());
  if ( jets.size()<3 ) return false;

  //
  // Compute variables
  //

  // MET "isolation"
  float metIso = 100.;
  for ( int iJet=0; iJet<4; ++iJet) {
    // FIXME: new deltaPhi function now (>2.1) accepts collection
    double deltaPhiAbs = fabs( reco::deltaPhi(jets[iJet].phi(),mets[0].phi()) );
    if ( metIso > deltaPhiAbs ) metIso = deltaPhiAbs;
  }
  setVariable("metDphi",metIso);

  // MET and leading jet deltaPhi (in degrees!)
  // FIXME: new deltaPhi function now (>2.1) accepts collection
  double dPhiJet2Met = TMath::RadToDeg()*fabs( reco::deltaPhi(jets[1].phi(),mets[0].phi()) ); 
  setVariable("dPhiJet2Met",dPhiJet2Met);

  // R1 & R2
  double dphi1 = reco::deltaPhi(jets[0].phi(),mets[0].phi());
  double dphi2 = reco::deltaPhi(jets[1].phi(),mets[0].phi());
  float R1 = sqrt( TMath::Power(dphi2,2) + TMath::Power(dphi1-TMath::Pi(),2) );
  float R2 = sqrt( TMath::Power(dphi1,2) + TMath::Power(dphi2-TMath::Pi(),2) );
  setVariable("R1",R1);
  setVariable("R2",R2);
  

  //
  // Perform final selection
  //
  if ( metIso      < metDphiMin_     ) return false;
  if ( dPhiJet2Met < dPhiJet2MetMin_ ) return false;
  if ( R1          < rDistJetsMin_   ) return false;
  if ( R2          < rDistJetsMin_   ) return false;

  return true;

}

//________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"

DEFINE_EDM_PLUGIN(EventSelectorFactory, MetJetEventSelector, "MetJetEventSelector");
