#ifndef Workspace_EventSelectors_h_
#define Workspace_EventSelectors_h_
/** Example code for the use of a SelectorSequence (i.e.,
 *  a series of selector modules) in an EDAnalyzer.
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Workspace/EventSelectors/interface/SelectorSequence.h"

#include <vector>

class SusySelectorExample : public edm::EDAnalyzer {
public:
  explicit SusySelectorExample(const edm::ParameterSet&);
  ~SusySelectorExample();
  
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  inline unsigned int nrOfSelectors () const {return selectors_.size();}

private:
  unsigned int nrEventTotal_;
  unsigned int nrEventSelected_;

  SelectorSequence selectors_;
  std::vector<std::string> selVarNames_;
};
#endif
