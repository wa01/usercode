// -*- C++ -*-
//
// Package:    GsfWriter
// Class:      GsfWriter
// 
/**\class GsfWriter GsfWriter.cc Workspace/GsfWriter/src/GsfWriter.cc

 Description: dump GSF mixture and simulated values to file.
   Selection (simulated di-electrons with opposite momenta):
   - either 1 or 2 reconstructed GSF tracks 
       (in the 2nd case: deltaR>1.5 between the two tracks)
   - exactly two simulated tracks from main vertex
   - sim / reco according to smallest deltaR (and requiring deltaR<1.5)
   Output: local parameters (q/p, dx/dz, dy/dz, x, y) on a plane
   - containing the simulated vertex and with
   -- local z || to global x or y (choosing the solution with the smaller angle 
        between the simulated transverse momentum and the axis)
   -- local y == global z
   Format of the output file. For each electron the following lines
   - local position (==(0/0/0)), local momentum and charge (simulation)
   - local parameters (simulation)
   - N (number of GSF components)
   - N lines with weight + local parameters 
*/
//
// Original Author:  Wolfgang Adam,40 4-A28,,
//         Created:  Thu Oct 20 15:42:28 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "DataFormats/GeometrySurface/interface/OpenBounds.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GsfTracking/interface/GsfMaterialEffectsUpdator.h"
#include "TrackingTools/GsfTracking/interface/GsfPropagatorWithMaterial.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/Math/interface/deltaR.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"

#include <fstream>
//
// class declaration
//


class GsfWriter : public edm::EDAnalyzer {
public:
  explicit GsfWriter(const edm::ParameterSet&);
  ~GsfWriter() {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) {}
  virtual void endRun(edm::Run const&, edm::EventSetup const&) {}
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

  void dumpTrack (const SimTrack&, const SimVertex&, const reco::TransientTrack&, const Propagator&);

  // ----------member data ---------------------------
private:
  ofstream outFile_;                  // output stream

  edm::InputTag gsfInputTag_;         // tag for GSF tracks
  edm::InputTag simTrackInputTag_;    // tag for simulated tracks
  edm::InputTag simVertexInputTag_;   // tag for simulated vertices

  TH1* hDRGsfGsf_;                    // control histogram (dR between GSF tracks)
  TH1* hDRGsfSim_;                    // control histogram (min dR between GSF track and simulation)
  TH1* hDPt_;                         // control histogram (deltaPt/Pt between GSF track and simulation)
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GsfWriter::GsfWriter(const edm::ParameterSet& iConfig) :
  outFile_(iConfig.getUntrackedParameter<std::string>("file","GsfWriter.dat").c_str()),
  gsfInputTag_(iConfig.getParameter<edm::InputTag>("gsfTracks")),
  simTrackInputTag_(iConfig.getParameter<edm::InputTag>("simTracks")),
  simVertexInputTag_(iConfig.getParameter<edm::InputTag>("simVertices")) {}


//
// member functions
//

// ------------ method called for each event  ------------
void
GsfWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //
   // TransientTrack builder
   //
   ESHandle<TransientTrackBuilder> ttBuilderHandle;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttBuilderHandle);
   const TransientTrackBuilder& ttBuilder = *ttBuilderHandle;
   //
   // Material effects
   //
   edm::ESHandle<GsfMaterialEffectsUpdator> matProducer;
   iSetup.get<TrackingComponentsRecord>().get("ElectronMaterialEffects",matProducer);
   //
   // Geom propagator
   //
   edm::ESHandle<Propagator> geomProducer;
   iSetup.get<TrackingComponentsRecord>().get("bwdAnalyticalPropagator",geomProducer);
   //
   // GSF propagator
   //
   GsfPropagatorWithMaterial propagator(*geomProducer.product(),*matProducer.product());
   //
   // Handles for GSF tracks and simulation
   //
   Handle<reco::GsfTrackCollection> gsfHandle;
   iEvent.getByLabel(gsfInputTag_,gsfHandle);
   //
   Handle<edm::SimTrackContainer> simTrackHandle;
   iEvent.getByLabel(simTrackInputTag_,simTrackHandle);
   //
   Handle<edm::SimVertexContainer> simVertexHandle;
   iEvent.getByLabel(simVertexInputTag_,simVertexHandle);
   //
   // for convenience: vector of pointers to GsfTracks
   //
   std::vector<const reco::GsfTrack*> gsfTrackPtrs;
   const reco::GsfTrackCollection& gsfTracks = *gsfHandle;
   edm::LogInfo("GsfWriter") << "Found " << gsfTracks.size() << " GSF tracks";
   for ( unsigned int i=0; i<gsfTracks.size(); ++i ) {
     const reco::GsfTrack& gsfTrack = gsfTracks[i];
     gsfTrackPtrs.push_back(&gsfTrack);
     edm::LogInfo("GsfWriter") << "  " << i << " " << gsfTrack.pt() << " " << gsfTrack.eta() 
			       << " " << gsfTrack.phi();
   }
   //
   // only accept events with a pair of well-separated tracks or a single track
   //
   if ( gsfTrackPtrs.size()==0 || gsfTrackPtrs.size()>2 ||
	(gsfTrackPtrs.size()==2 &&
	 reco::deltaR<reco::GsfTrack,reco::GsfTrack>(*gsfTrackPtrs[0],*gsfTrackPtrs[1])<1.5 ) ) {
     edm::LogInfo("GsfWriter") << "Event rejected: gsf track count / dR";
     return;
   }
   //
   // create vector of pointers to SimTracks from the first generated vertex
   //   (reject secondaries)
   //   
   std::vector<const SimTrack*> simTrackPtrs;
   const edm::SimTrackContainer& simTracks = *simTrackHandle;
   edm::LogInfo("GsfWriter") << "Found " << simTracks.size() << " SIM tracks";
   for ( unsigned int i=0; i<simTracks.size(); ++i ) {
     const SimTrack& simTrack = simTracks[i];
     if ( simTrack.vertIndex()==0 ) {
       simTrackPtrs.push_back(&simTrack);
       edm::LogInfo("GsfWriter") << "  " << i
		 << " " << simTrack.momentum().pt() 
		 << " " << simTrack.momentum().eta() 
		 << " " << simTrack.momentum().phi();
     }
   }
   //
   // only accept events with two simulated (primary) electrons
   //
   if ( simTrackPtrs.size()!=2 ) {
     edm::LogInfo("GsfWriter") << "Event rejected: #sim electrons != 2";
     return;
   }
   //
   // (trivial) matching of GSF tracks to SimTracks (by deltaR)
   //
   std::vector<const reco::GsfTrack*> matchedGsfPtrs(simTrackPtrs.size(),0);
   for ( unsigned int i=0; i<simTrackPtrs.size(); ++i ) {
     double drMin(1.e30);
     size_t gsfMin(gsfTrackPtrs.size());
     for ( unsigned int j=0; j<gsfTrackPtrs.size(); ++j ) {
       double dr = 
	 reco::deltaR<math::XYZTLorentzVectorD,reco::GsfTrack>(simTrackPtrs[i]->momentum(),*gsfTrackPtrs[j]);
       if ( dr<drMin ) {
	 drMin = dr;
	 gsfMin = j;
       }
     }
     edm::LogInfo("GsfWriter") << "Lowest dR for simTrack " << i << " is GsfTrack " << gsfMin
	       << " with dR = " << drMin;
     if ( gsfMin<gsfTrackPtrs.size() ) {
       if ( drMin<1.5 )  matchedGsfPtrs[i] = gsfTrackPtrs[gsfMin];
       hDRGsfSim_->Fill(drMin);
       hDPt_->Fill(gsfTrackPtrs[gsfMin]->pt()/simTrackPtrs[i]->momentum().pt()-1.);
     }
   }
   //
   // dump information for SimTrack / GsfTrack pairs (after conversion
   //   of the GsfTrack to a TransientTrack for access to TrajectoryStateOnSurface
   //   with the full GSF mixture)
   for ( unsigned int i=0; i<simTrackPtrs.size(); ++i ) {
     if ( matchedGsfPtrs[i] ) 
       dumpTrack(*simTrackPtrs[i],(*simVertexHandle)[0],ttBuilder.build(*matchedGsfPtrs[i]),propagator);
   }

//    const edm::SimVertexContainer& simVertices = *simVertexHandle;
//    edm::LogInfo("GsfWriter") << "Found " << simVertices.size() << " SIM vertices";
//    for ( unsigned int i=0; i<simVertices.size(); ++i ) {
//      const SimVertex& simVertex = simVertices[i];
//      edm::LogInfo("GsfWriter") << "  " << i << " " << simVertex.position();
//    }

}


// ------------ method called once each job just before starting event loop  ------------
void 
GsfWriter::beginJob()
{
  //
  // definition of histograms
  //
  edm::Service<TFileService> fs;
  hDRGsfGsf_ = fs->make<TH1F>("hDRGsfGsf","min. #Delta R gsf track / gsf track",200,0.,5);
  hDRGsfSim_ = fs->make<TH1F>("hDRGsfSim","min. #Delta R sim track / gsf track",200,0.,5);
  hDPt_ = fs->make<TH1F>("hDPt","#Delta p_{T}/p_{T}",100,-2.,2);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GsfWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
GsfWriter::dumpTrack (const SimTrack& simTrack, const SimVertex& simVertex,
		      const reco::TransientTrack& tTrack, const Propagator& propagator)
{
  //
  // SimTrack origin and momentum in global co-ordinates
  //
  GlobalPoint simPos(simVertex.position().X(),
		     simVertex.position().Y(),
		     simVertex.position().Z());
  GlobalVector simMom(simTrack.momentum().x(),
		      simTrack.momentum().y(),
		      simTrack.momentum().z());
  //
  // Plane containing the SimTrack origin and || to global x-z or y-z
  //   (choice according to smallest angle between local z and transverse momentum vector)
  //
  BoundPlane::PositionType planePos(simPos);
  static const GlobalVector xGlob(1.,0.,0.);
  static const GlobalVector yGlob(0.,1.,0.);
  static const GlobalVector zGlob(0.,0.,1.);
  bool horizontal = fabs(simTrack.momentum().x())>fabs(simTrack.momentum().y());
  const GlobalVector& xAxis = horizontal ? yGlob : xGlob;
  BoundPlane::RotationType planeRot(xAxis,zGlob);
  BoundPlane::BoundPlanePointer plane = BoundPlane::build(planePos,planeRot,OpenBounds());
  edm::LogInfo("GsfWriter") 
//     << "Plane axes are:" << std::endl
// 			    << plane->toGlobal(LocalVector(1.,0.,0.)) << std::endl
// 			    << plane->toGlobal(LocalVector(0.,1.,0.)) << std::endl
// 			    << plane->toGlobal(LocalVector(0.,0.,1.)) << std::endl
			    << " Local position : " << plane->toLocal(simPos) << std::endl
			    << " Local momentum : " << plane->toLocal(simMom);
  //
  // Just to be sure: extrapolate SimTrack to plane (should give same result as above)
  FreeTrajectoryState simFTS(simPos,simMom,simTrack.charge(),propagator.magneticField());
  TrajectoryStateOnSurface simState = propagator.propagate(simFTS,*plane);
  if ( !simState.isValid() ) {
    edm::LogError("GsfWriter") << "simState is invalid";
    return;
  }
  edm::LogInfo("GsfWriter") << simState.localPosition() << " " << simState.localMomentum();
  //
  // write simulation parameters to output file
  //
  // position, momentum in global co-ordinates and charge
  outFile_ << simState.localPosition().x() << " "
	   << simState.localPosition().y() << " "
	   << simState.localPosition().z() << " "
	   << simState.localMomentum().perp() << " "
	   << simState.localMomentum().eta() << " "
	   << simState.localMomentum().phi() << " "
	   << simState.charge() << std::endl;
  // local parameters
  AlgebraicVector5 simLocPar(simState.localParameters().vector());
  for ( size_t i=0; i<5; ++i ) outFile_ << " " << simLocPar[i];
  outFile_ << std::endl; 
  //
  // Propagation of the GsfTrack from the innermost measured point to the plane
  //
  TrajectoryStateOnSurface gsfState = 
    propagator.propagate(tTrack.innermostMeasurementState(),*plane);
  if ( !gsfState.isValid() ) {
    edm::LogError("GsfWriter") << "gsfState is invalid";
    return;
  }
  //
  // write GSF parameters to the output file:
  //   first #components, then for each component weight and local parameters
  //
  std::vector<TrajectoryStateOnSurface> components(gsfState.components());
  outFile_ << components.size() << std::endl;
  for ( size_t i=0; i<components.size(); ++i ) {
    outFile_ << " " << components[i].weight() << " ";
    AlgebraicVector5 gsfLocPar(components[i].localParameters().vector());
    for ( size_t i=0; i<5; ++i ) outFile_ << " " << gsfLocPar[i];
    outFile_ << std::endl;
  }

  

}

//define this as a plug-in
DEFINE_FWK_MODULE(GsfWriter);
