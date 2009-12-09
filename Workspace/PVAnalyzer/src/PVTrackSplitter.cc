// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <vector>

//
// class decleration
//

class PVTrackSplitter : public edm::EDFilter {
public:
  explicit PVTrackSplitter(const edm::ParameterSet&);
  ~PVTrackSplitter();
  
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
private:
//   edm::InputTag bsTag_;
  edm::InputTag pvTag_;
  
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
PVTrackSplitter::PVTrackSplitter(const edm::ParameterSet& iConfig) :
//   bsTag_(iConfig.getParameter<edm::InputTag>("beamSpot")),
  pvTag_(iConfig.getParameter<edm::InputTag>("primaryVertices")) {
   //now do what ever initialization is needed
  produces<reco::TrackCollection>("pvSelected");
  produces<reco::TrackCollection>("pvRejected");
}


PVTrackSplitter::~PVTrackSplitter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
bool
PVTrackSplitter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
//   Handle<reco::BeamSpot> beamspot;
//   iEvent.getByLabel(bsTag_,beamspot);
  
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(pvTag_,vertices);

//   bsPar_[0] = beamspot->x0();
//   bsPar_[1] = beamspot->y0();
//   bsPar_[2] = beamspot->z0();

  std::auto_ptr<reco::TrackCollection> acceptedTracks(new reco::TrackCollection());
  std::auto_ptr<reco::TrackCollection> rejectedTracks(new reco::TrackCollection());

  for ( size_t iv=0; iv<vertices->size(); ++iv ) {
    const reco::Vertex& vertex = (*vertices)[iv];
    if ( vertex.isFake() )  continue;

    unsigned int nacc(0);
    unsigned int nrej(0);
    for ( reco::Vertex::trackRef_iterator it=vertex.tracks_begin();
	  it!=vertex.tracks_end(); ++it ) {
      double wgt = vertex.trackWeight(*it);
      if ( wgt>0.5 )  acceptedTracks->push_back(**it);
      else            rejectedTracks->push_back(**it);
      if ( wgt>0.5 )  ++nacc;
      else            ++nrej;
    }
    std::cout << "  vertex " << iv << " acc / rej = " << nacc << " " << nrej << std::endl;
  }

  iEvent.put(acceptedTracks,"pvSelected");
  iEvent.put(rejectedTracks,"pvRejected");

  return vertices->size()>1;
}


// ------------ method called once each job just before starting event loop  ------------
void 
PVTrackSplitter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PVTrackSplitter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PVTrackSplitter);
