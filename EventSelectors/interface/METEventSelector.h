#ifndef Workspace_METEventSelector_h_
#define Workspace_METEventSelector_h_
/** Trivial example for a MET selector.
 *  To be modified for analysis!
 */
// Original author: W. Adam, 10/4/08

// system include files
#include <memory>

// user include files
#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <vector>
#include <string>

class METEventSelector : public SusyEventSelector {
public:
  METEventSelector (const edm::ParameterSet&);
  virtual bool select (const edm::Event&) const;
  virtual ~METEventSelector () {}
private:
  edm::InputTag metTag_;      ///< tag for input collection
  float minMet_;              ///< lower cut on MET
};
#endif
