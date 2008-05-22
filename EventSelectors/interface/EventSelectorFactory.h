#ifndef SusyAnalysis_EventSelector_EventSelectorFactory_H
#define SusyAnalysis_EventSelector_EventSelectorFactory_H

/// Plugin factory for SUSY event selector modules. 
///
/// Original author: W. Adam, 10/4/08
/// $Id: $

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "SusyAnalysis/EventSelector/interface/SusyEventSelector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

typedef edmplugin::PluginFactory< SusyEventSelector* (const edm::ParameterSet&) > EventSelectorFactory;

#endif
