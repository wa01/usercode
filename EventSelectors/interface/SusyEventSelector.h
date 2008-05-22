#ifndef SusyAnalysis_SusyEventSelector_h_
#define SusyAnalysis_SusyEventSelector_h_
///
/// Base class for event selection modules for SUSY analysis.
///
/// Original author: W. Adam, 10/4/08
///
/// $Id: $


// system include files
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// default (extreme) value for a variable
namespace susy {
  static const double DEFAULT_VALUE = -1.e30;
}

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

  /// number of cached variables
  virtual size_t numberOfVariables () const {return variableNames_.size();}
  /// variable names
  virtual std::vector<std::string> variableNames () const {return variableNames_;}
  /// variable values
  virtual std::vector<double> values () const {return variableValues_;}
  /// variable value by name
  virtual double value (const std::string& name) const;

protected:
  /// definition of a cached variable (returns index)
  size_t defineVariable (const std::string& name);
  /// reset of all variables
  void resetVariables () const {
    fill(variableValues_.begin(),variableValues_.end(),susy::DEFAULT_VALUE);
  }
  /// setting the value of a variable by index
  void setVariable (size_t index, const double& value) const;
  /// setting the value of a variable by name
  void setVariable (const std::string& name, const double& value) const;

private:
  /// index from name
  size_t variableIndex (const std::string& name) const;

private:
  std::string name_;

  std::vector<std::string> variableNames_;
  mutable std::vector<double> variableValues_;
};
#endif
