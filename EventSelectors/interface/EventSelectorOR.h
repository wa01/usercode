#ifndef SusyAnalysis_EventSelectorOR_h_
#define SusyAnalysis_EventSelectorOR_h_
///
/// Combination of selectors by logical OR.
///
/// Original author: W. Adam, 10/4/08
///
/// $Id: $

// system include files
#include <memory>

// user include files
#include "SusyAnalysis/EventSelector/interface/CombinedEventSelector.h"
#include "SusyAnalysis/EventSelector/interface/SelectorSequence.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <vector>

class SelectorSequence;

class EventSelectorOR : public CombinedEventSelector {
public:
  EventSelectorOR ();
  EventSelectorOR (const edm::ParameterSet&);
  virtual ~EventSelectorOR () {}

  /// Selection: OR of all selectors
  virtual bool select (const edm::Event&) const;

};
#endif
