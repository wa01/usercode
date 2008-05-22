#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"

#include "SusyAnalysis/EventSelector/interface/EventSelectorAND.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorOR.h"
#include "SusyAnalysis/EventSelector/interface/HLTEventSelector.h"
#include "SusyAnalysis/EventSelector/interface/METEventSelector.h"
#include "SusyAnalysis/EventSelector/interface/JetEventSelector.h"
#include "SusyAnalysis/EventSelector/interface/BJetEventSelector.h"

DEFINE_EDM_PLUGIN(EventSelectorFactory, EventSelectorAND, "EventSelectorAND");
DEFINE_EDM_PLUGIN(EventSelectorFactory, EventSelectorOR, "EventSelectorOR");
DEFINE_EDM_PLUGIN(EventSelectorFactory, HLTEventSelector, "HLTEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, METEventSelector, "METEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, JetEventSelector, "JetEventSelector");
DEFINE_EDM_PLUGIN(EventSelectorFactory, BJetEventSelector, "BJetEventSelector");
