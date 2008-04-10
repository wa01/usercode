#ifndef Workspace_EventSelectorOR_h_
#define Workspace_EventSelectorOR_h_
/** Combination of selectors by logical OR.
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

class EventSelectorOR : public SusyEventSelector {
public:
  EventSelectorOR ();
  EventSelectorOR (const edm::ParameterSet&);
  virtual bool select (const edm::Event&) const;
  virtual ~EventSelectorOR () {}
private:
  SelectorSequence sequence_;
};
#endif
