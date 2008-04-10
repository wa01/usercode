#ifndef Workspace_EventSelectorAND_h_
#define Workspace_EventSelectorAND_h_
/** Combination of selectors by logical AND.
 */
// Original author: W. Adam, 10/4/08

// system include files
#include <memory>

// user include files
#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "Workspace/EventSelectors/interface/SelectorSequence.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <vector>

class SelectorSequence;

class EventSelectorAND : public SusyEventSelector {
public:
  EventSelectorAND ();
  EventSelectorAND (const edm::ParameterSet&);
  virtual bool select (const edm::Event&) const;
  virtual ~EventSelectorAND () {}
private:
  SelectorSequence sequence_;
};
#endif
