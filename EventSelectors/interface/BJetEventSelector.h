#ifndef Workspace_BJetEventSelector_h_
#define Workspace_BJetEventSelector_h_
/** Trivial example for a b-jet selector.
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

class BJetEventSelector : public SusyEventSelector {
public:
  BJetEventSelector (const edm::ParameterSet&);
  virtual bool select (const edm::Event&) const;
  virtual ~BJetEventSelector () {}
private:
  edm::InputTag jetTag_;                ///< tag for input collection
  std::string tagLabel_;                ///< b-tag label
  std::vector<double> minTag_;          ///< lower cuts on discriminator (defines also min. #jets)
};
#endif
