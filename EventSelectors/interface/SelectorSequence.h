#ifndef EventSelector_SelectorSequence_h_
#define EventSelector_SelectorSequence_h_

/** Class regrouping a series of selector modules. The modules are instantiated 
 *  according to specifications in the configuration files.
 */
// Original author: W. Adam, 10/4/08

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <string>

class SusyEventSelector;

class SelectorSequence {
public:
  explicit SelectorSequence (const edm::ParameterSet&);
  explicit SelectorSequence (const std::vector<std::string>& sequence,
			     const edm::ParameterSet& selectors);
  ~SelectorSequence();
  
  /// number of selectors
  inline size_t size () const {return selectors_.size();}
  /// names of selectors
  const std::vector<std::string>& selectorNames () const;
  /// access to individual selectors
  const std::vector<const SusyEventSelector*>& selectors () const {
    return selectors_;
  }
  /// selector index from name 
  size_t selectorIndex (const std::string& selectorName) const;
  /// selector name from index
  std::string selectorName (size_t index) const;

  /// selector results (in the same order as defined in selectors()
  std::vector<bool> decisions (const edm::Event& event) const;

  /// selector result by index
  bool decision (const edm::Event& event, size_t index) const;

  /// selector result by name
  bool decision (const edm::Event& event, 
		 const std::string& selectorName) const;

  /// global decision (AND of all selectors)
  bool globalDecision (const edm::Event& event) const;

  /// complementary selection result by index (i.e., excluding one selector)
  bool complementaryDecision (const edm::Event& event, size_t index) const;

  /// complementary selection result by name (i.e., excluding one selector)
  bool complementaryDecision (const edm::Event& event, 
			      const std::string& selectorName) const;

  /** cumulative selection result by index (AND of the selectors in the list 
   *  from 0 to - and including - "index") */
  bool cumulativeDecision (const edm::Event& event, size_t index) const;

  /** cumulative selection result by name (AND of the selectors in the list 
   *  up to and including the one identified by "name") */
  bool cumulativeDecision (const edm::Event& event, 
			   const std::string& selectorName) const;

  /// number of variables
  size_t numberOfVariables () const;
  /// get all variable names
  std::vector<std::string> variableNames () const;
  /// get all variable values
  std::vector<double> values () const;
  /// get value of a specific variable
  double value (const std::string& selectorName,
		const std::string& variableName) const;

private:
  void createSelectors (const std::vector<std::string>& sequence,
			const edm::ParameterSet& selectors);
  inline bool newEvent (const edm::Event& event) const {return event.id()!=currentEventId_;}
  
private:
  std::vector<const SusyEventSelector*> selectors_;
  mutable std::vector<std::string> selectorNames_;

  mutable edm::EventID currentEventId_;
  mutable std::vector<bool> currentDecisions_;
};

#endif
