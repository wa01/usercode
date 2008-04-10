#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

DEFINE_SEAL_MODULE();

#include "Workspace/EventSelectors/interface/EventSelectorFactory.h"

#include "Workspace/EventSelectors/interface/EventSelectorAND.h"
#include "Workspace/EventSelectors/interface/EventSelectorOR.h"
#include "Workspace/EventSelectors/interface/HLTEventSelector.h"
#include "Workspace/EventSelectors/interface/METEventSelector.h"
#include "Workspace/EventSelectors/interface/JetEventSelector.h"
#include "Workspace/EventSelectors/interface/BJetEventSelector.h"

DEFINE_EDM_PLUGIN(EventSelectorFactory, EventSelectorAND, "EventSelectorAND");
DEFINE_EDM_PLUGIN(EventSelectorFactory, EventSelectorOR, "EventSelectorOR");
DEFINE_EDM_PLUGIN(EventSelectorFactory, HLTEventSelector, "HLTEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, METEventSelector, "METEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, JetEventSelector, "JetEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, BJetEventSelector, "BJetEventSelector");
