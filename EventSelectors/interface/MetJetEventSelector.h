#ifndef SusyAnalysis_AnalysisSkeleton_MetJetEventSelector_h_
#define SusyAnalysis_AnalysisSkeleton_MetJetEventSelector_h_
/// Selection on delta phi MET / jets.
///
/// Performs a number of selections based on MET and jets topology:
/// - at least 3 jets (otherwise the variables are not calculated)
/// - minimum delta-phi between MET and any of the (up to 4) leading jets
/// - minimum delta-phi between MET and second leading jet 
/// - minimum R1 and R2 distances (between MET and 2 leading jets)
///
/// $Id $


// system include files
#include <memory>

// user include files
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class MetJetEventSelector : public SusyEventSelector {
public:
  MetJetEventSelector (const edm::ParameterSet&);
  virtual bool select (const edm::Event&) const;
  virtual ~MetJetEventSelector () {}
private:
  edm::InputTag metTag_; ///< tag for input collection
  edm::InputTag jetTag_; ///< tag for input collection

  float metDphiMin_;     ///< Minimum angle between MET and any leading jet
  float dPhiJet2MetMin_; ///< Minimum angle between MET and 2nd leading jet
  float rDistJetsMin_;   ///< Minimum R1 and R2

};
#endif
