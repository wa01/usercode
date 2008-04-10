#ifndef Workspace_EventSelector_EventSelectorFactory_H
#define Workspace_EventSelector_EventSelectorFactory_H

/** Plugin factory for SUSY event selector modules. */
// Original author: W. Adam, 10/4/08

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "Workspace/EventSelectors/interface/SusyEventSelector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

typedef edmplugin::PluginFactory< SusyEventSelector* (const edm::ParameterSet&) > EventSelectorFactory;

#endif
