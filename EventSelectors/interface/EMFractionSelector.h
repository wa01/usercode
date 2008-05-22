#ifndef SusyAnalysis_EventSelector_EMFractionSelector_h_
#define SusyAnalysis_EventSelector_EMFractionSelector_h_
///
/// Selection based on EM fraction.
///
/// Calculates the total EM fraction of all jets
///
/// $Id: $

// system include files
#include <memory>

// user include files
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class EMFractionSelector : public SusyEventSelector {
public:
  EMFractionSelector (const edm::ParameterSet&);
  virtual bool select (const edm::Event&) const;
  virtual ~EMFractionSelector () {}

private:
  edm::InputTag jetTag_; ///< tag for input collection
  float minFraction_;    ///< lower EM fraction cut
  float maxEta_;         ///< upper |eta| cut on all jets
  float minPt_;          ///< lower pt cut on all jets

};
#endif
