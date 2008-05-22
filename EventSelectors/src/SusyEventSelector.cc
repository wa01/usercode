#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"

#include "FWCore/Utilities/interface/Exception.h"

using namespace std;

double
SusyEventSelector::value (const std::string& name) const
{
  return variableValues_[variableIndex(name)];
}

size_t
SusyEventSelector::defineVariable (const std::string& name)
{
  size_t result = variableNames_.size();
  variableNames_.push_back(name);
  variableValues_.push_back(susy::DEFAULT_VALUE); // Initialise to extreme value...
  return result;
}

void
SusyEventSelector::setVariable (size_t index, const double& value) const
{
  if ( index>=variableValues_.size() )
    throw cms::Exception("Error") << "index out of bounds " << index;

  variableValues_[index] = value;
}

void
SusyEventSelector::setVariable (const std::string& name, 
				const double& value) const
{
  variableValues_[variableIndex(name)] = value;
}

size_t
SusyEventSelector::variableIndex (const std::string& name) const
{
  vector<string>::const_iterator i =
    find(variableNames_.begin(),variableNames_.end(),name);
  if ( i==variableNames_.end() )
    throw cms::Exception("Error") << "no variable with name " << name;

  return i-variableNames_.begin();
}
