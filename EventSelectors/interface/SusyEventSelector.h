#ifndef Workspace_SusyEventSelector_h_
#define Workspace_SusyEventSelector_h_
/** Base class for event selection modules for SUSY analysis.
 */
// Original author: W. Adam, 10/4/08

// system include files
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class SusyEventSelector {
public:
  SusyEventSelector () {}
  SusyEventSelector (const edm::ParameterSet& iConfig) {
    std::string selector = iConfig.getParameter<std::string>("selector");
    name_ = iConfig.getUntrackedParameter<std::string>("name",selector);
  }
  virtual ~SusyEventSelector () {}
  /// name of the module (from configuration)
  const std::string& name () const {return name_;}
  /// decision of the selector module
  virtual bool select (const edm::Event&) const = 0;

private:
  std::string name_;
};
#endif
