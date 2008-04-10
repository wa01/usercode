#ifndef Workspace_JetEventSelector_h_
#define Workspace_JetEventSelector_h_
/** Trivial example for a Jet selector.
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

class JetEventSelector : public SusyEventSelector {
public:
  JetEventSelector (const edm::ParameterSet&);
  virtual bool select (const edm::Event&) const;
  virtual ~JetEventSelector () {}
private:
  edm::InputTag jetTag_;      ///< tag for input collection
  std::vector<double> minEt_; ///< lower Et cuts (defines also min. #jets)
};
#endif
