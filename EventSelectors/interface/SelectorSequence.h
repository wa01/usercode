#ifndef EventSelectors_SelectorSequence_h_
#define EventSelectors_SelectorSequence_h_

/** Class regrouping a series of selector modules. The modules are instantiated 
 *  according to specifications in the configuration files.
 */
// Original author: W. Adam, 10/4/08

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

class SusyEventSelector;

class SelectorSequence {
public:
  explicit SelectorSequence (const edm::ParameterSet&);
  ~SelectorSequence();
  
  /// access to individual selectors
  const std::vector<const SusyEventSelector*>& selectors () const {
    return selectors_;
  }
  /// number of selectors
  unsigned int size () const {return selectors_.size();}
  /// selector results (in the same order as defined in selectors()
  std::vector<bool> decisions (const edm::Event&) const;
  

private:
  std::vector<const SusyEventSelector*> selectors_;
};

#endif
