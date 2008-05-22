#ifndef SusyAnalysis_CombinedEventSelector_h_
#define SusyAnalysis_CombinedEventSelector_h_
///
/// Combination of selectors by logical AND.
///
/// Original author: W. Adam, 10/4/08
///
/// $Id: $

// system include files
#include <memory>

// user include files
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "SusyAnalysis/EventSelector/interface/SelectorSequence.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <vector>

class SelectorSequence;

class CombinedEventSelector : public SusyEventSelector {
public:
  CombinedEventSelector ();
  CombinedEventSelector (const edm::ParameterSet&);
  // Needs to be implemented in concrete classes
  virtual bool select (const edm::Event&) const = 0;
  virtual ~CombinedEventSelector () {}

  // re-implement access to variables
  /// number of cached variables
  virtual size_t numberOfVariables () const { return sequence_.numberOfVariables(); }
  /// variable names
  virtual std::vector<std::string> variableNames () const { return sequence_.variableNames(); }
  /// variable values
  virtual std::vector<double> values () const { return sequence_.values(); }
  /// variable value by name
  virtual double value (const std::string& name) const;

protected:
  SelectorSequence sequence_; ///< Sequence of combined selectors

};
#endif
