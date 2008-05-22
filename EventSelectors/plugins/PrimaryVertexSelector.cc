#include "SusyAnalysis/EventSelector/interface/PrimaryVertexEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//________________________________________________________________________________________
PrimaryVertexEventSelector::PrimaryVertexEventSelector (const edm::ParameterSet& pset) :
  SusyEventSelector(pset),
  vertexTag_( pset.getParameter<edm::InputTag>("src") )
{ 

  defineVariable("numberOfVertices");

}

//________________________________________________________________________________________
bool
PrimaryVertexEventSelector::select (const edm::Event& event) const
{
  // reset cached variables
  resetVariables();
  // Get the vertices
  edm::Handle< std::vector<reco::Vertex> > vertexHandle;
  event.getByLabel(vertexTag_, vertexHandle);
  if ( !vertexHandle.isValid() ) {
    edm::LogWarning("PrimaryVertexEventSelector") 
      << "No Vertex results for InputTag " << vertexTag_;
    return false;
  }

  setVariable( "numberOfVertices", (*vertexHandle).size() );

  return (*vertexHandle).size()>=1;

}

//________________________________________________________________________________________
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "SusyAnalysis/EventSelector/interface/EventSelectorFactory.h"

DEFINE_EDM_PLUGIN(EventSelectorFactory, PrimaryVertexEventSelector, "PrimaryVertexEventSelector");
