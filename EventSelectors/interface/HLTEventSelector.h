#ifndef SusyAnalysis_HLTEventSelector_h_
#define SusyAnalysis_HLTEventSelector_h_

/// Trivial example for a HLT selector.
///
/// To be modified for analysis!
///
/// Original author: W. Adam, 10/4/08
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

#include <vector>
#include <string>

class HLTEventSelector : public SusyEventSelector {
public:
  HLTEventSelector (const edm::ParameterSet&);
  virtual bool select (const edm::Event&) const;
  virtual ~HLTEventSelector () {}
private:
  edm::InputTag triggerResults_;        ///< tag for input collection
  std::vector<std::string> pathNames_;  ///< trigger path names (ORed)
};
#endif
