// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <vector>
#include <string>

//
// class decleration
//

class PVAnalyzer : public edm::EDAnalyzer {
public:
  explicit PVAnalyzer(const edm::ParameterSet&);
  ~PVAnalyzer();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
private:
  edm::InputTag bsTag_;
  edm::InputTag pvTag_;

  TH1* h_mult_;
  TH1* h_tkmult_;
  TH1* h_wgts_;
  TH1* h_ndf_;
  TH1* h_chi2_;
  TH1* h_vx_;
  TH1* h_vy_;
  TH1* h_vz_;
  TH2* h_vxy_;
  TH2* h_vzx_;
  TH2* h_vzy_;
  TH1* h_evx_;
  TH1* h_evy_;
  TH1* h_evz_;
  TH1* h_dvx_;
  TH1* h_dvy_;
  TH1* h_dzpvpv_;

  TTree* tree_;

  float bsPar_[3];
  unsigned short npv_;
  float dzpvpv_;
  unsigned short pvTkMult_;
  float pvPars_[6];
  float pvQual_[5];
  enum { NTKMAX = 200 };
  int pvNTk_;
  unsigned short test_[NTKMAX];
  unsigned short tkQualI_[3][NTKMAX];
  float tkQualF_[3][NTKMAX];
  float tkPar_[5][NTKMAX];
  float tkErr_[5][NTKMAX];
  
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
PVAnalyzer::PVAnalyzer(const edm::ParameterSet& iConfig) :
  bsTag_(iConfig.getParameter<edm::InputTag>("beamSpot")),
  pvTag_(iConfig.getParameter<edm::InputTag>("primaryVertices")) {
   //now do what ever initialization is needed

}


PVAnalyzer::~PVAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PVAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  Handle<reco::BeamSpot> beamspot;
  iEvent.getByLabel(bsTag_,beamspot);
  
  Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(pvTag_,vertices);

//     if ( !beamspot.isValid() )  std::cout << "no beamspot" << std::endl;
//     if ( !beamspot.isValid() )  return;
//     std::cout << "beamspot x / y = " 
// 	      << beamspot->x0()
// 	      << " " << beamspot->y0() 
// 	      << std::endl;
  bsPar_[0] = beamspot->x0();
  bsPar_[1] = beamspot->y0();
  bsPar_[2] = beamspot->z0();
  reco::Track::Point bsPos(beamspot->x0(),beamspot->y0(),beamspot->z0());
  
  npv_ = 0;
  dzpvpv_ = 999999.;
  pvTkMult_ = 0;
  for ( size_t i=0; i<6; ++i )  pvPars_[i] = 999999.;
  for ( size_t i=0; i<5; ++i )  pvQual_[i] = -1.;
  pvNTk_ = 0;
  for ( size_t iv=0; iv<vertices->size(); ++iv ) {
    const reco::Vertex& vertex = (*vertices)[iv];
    if ( vertex.isFake() )  continue;
    //    if ( vertex.tracksSize()<=2 )  continue;
    
    ++npv_;
    if ( iv==0 ) {
      pvTkMult_ = vertex.tracksSize();
      h_tkmult_->Fill(pvTkMult_);
      pvPars_[0] = vertex.x();
      pvPars_[1] = vertex.y();
      pvPars_[2] = vertex.z();
      pvPars_[3] = vertex.xError();
      pvPars_[4] = vertex.yError();
      pvPars_[5] = vertex.zError();
      pvQual_[0] = vertex.ndof();
      pvQual_[1] = vertex.ndof()>0. ? vertex.chi2()/vertex.ndof() : 0.;
      pvQual_[2] = 0.;
      pvQual_[3] = 0.;
      pvQual_[4] = 0.;      
      for ( reco::Vertex::trackRef_iterator it=vertex.tracks_begin();
	    it!=vertex.tracks_end(); ++it ) {
	double wgt = vertex.trackWeight(*it);
	pvQual_[4] += wgt;
	h_wgts_->Fill(wgt);
	pvQual_[2] += (**it).pt();
	pvQual_[3] += wgt*(**it).pt();

	if ( pvNTk_<NTKMAX ) {
	  tkQualI_[0][pvNTk_] = (**it).hitPattern().numberOfValidPixelHits();
	  tkQualI_[1][pvNTk_] = (**it).hitPattern().numberOfValidStripHits();
	  tkQualI_[2][pvNTk_] = (**it).numberOfLostHits();
	  tkQualF_[0][pvNTk_] = (**it).ndof();
	  tkQualF_[1][pvNTk_] = (**it).normalizedChi2();
	  tkQualF_[2][pvNTk_] = wgt;
	  tkPar_[0][pvNTk_] = (**it).dxy(vertex.position());
	  tkPar_[1][pvNTk_] = (**it).dz(vertex.position());
	  tkPar_[2][pvNTk_] = (**it).pt();
	  tkPar_[3][pvNTk_] = (**it).eta();
	  tkPar_[4][pvNTk_] = (**it).phi();
	  tkErr_[0][pvNTk_] = (**it).dxyError();
	  tkErr_[1][pvNTk_] = (**it).dzError();
	  tkErr_[2][pvNTk_] = (**it).ptError();
	  tkErr_[3][pvNTk_] = (**it).etaError();
	  tkErr_[4][pvNTk_] = (**it).phiError();
	  ++pvNTk_;
	}
      }
      h_ndf_->Fill(vertex.ndof());
      h_chi2_->Fill(vertex.ndof()>0?vertex.chi2()/vertex.ndof():0.);
      h_vx_->Fill(vertex.x());
      h_vy_->Fill(vertex.y());
      h_vz_->Fill(vertex.z());
      h_vxy_->Fill(vertex.x(),vertex.y());
      h_vzx_->Fill(vertex.z(),vertex.x());
      h_vzy_->Fill(vertex.z(),vertex.y());
      h_evx_->Fill(vertex.xError());
      h_evy_->Fill(vertex.yError());
      h_evz_->Fill(vertex.zError());
      h_dvx_->Fill(vertex.x()-beamspot->x0());
      h_dvy_->Fill(vertex.y()-beamspot->y0());
    }

    if ( iv==1 ) {
      dzpvpv_ = vertex.z() - (*vertices)[0].z();
      h_dzpvpv_->Fill(dzpvpv_);
    }
  }
  h_mult_->Fill(npv_);

  tree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
PVAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

  h_mult_ = fs->make<TH1F>("pvmult","PV multiplicity",10,-0.5,9.5);

  h_tkmult_ = fs->make<TH1F>("pvtkmult","PV track multiplicity",100,-0.5,99.5);
  h_wgts_ = fs->make<TH1F>("pvwgts","PV track weights",100,0.,1.);
  h_ndf_ = fs->make<TH1F>("pvndf","PV ndf",100,0.,200.);
  h_chi2_ = fs->make<TH1F>("pvchi2","PV norm chi2",100,0.,20.);

  h_vx_ = fs->make<TH1F>("pvx","PV x position",100,-1.,1.);
  h_vy_ = fs->make<TH1F>("pvy","PV y position",100,-1.,1.);
  h_vz_ = fs->make<TH1F>("pvz","PV z position",100,-50.,50.);
  h_vxy_ = fs->make<TH2F>("pvxy","PV x-y position",50,-1.,1.,50,-1.,1.);
  h_vzx_ = fs->make<TH2F>("pvzx","PV z-x position",50,-50.,50.,50,-1.,1.);
  h_vzy_ = fs->make<TH2F>("pvzy","PV z-y position",50,-50.,50.,50,-1.,1.);

  h_evx_ = fs->make<TH1F>("pvex","PV x error",100,0.,0.2);
  h_evy_ = fs->make<TH1F>("pvey","PV y error",100,0.,0.2);
  h_evz_ = fs->make<TH1F>("pvez","PV z error",100,0.,0.5);

  h_dvx_ = fs->make<TH1F>("pvdvx","dx PV-BS",100,-1.,1.);
  h_dvy_ = fs->make<TH1F>("pvdvy","dy PV-BS",100,-1.,1.);

  h_dzpvpv_ = fs->make<TH1F>("dzpvpv","dz PV1-PV2",100,-30.,30.);

  tree_ = fs->make<TTree>("pv","Primary Vertices");
  tree_->Branch("beamspot",bsPar_,"bsx/F:bsy/F:bsz/F");
  tree_->Branch("npv",&npv_,"npv/s");
  tree_->Branch("dzpvpv",&dzpvpv_,"dzpvpv/F");
  tree_->Branch("pvmult",&pvTkMult_,"tkmult/s");
  tree_->Branch("pvpos",pvPars_,"x/F:y/F:z/F:ex/F:ey/F:ez/F");
  tree_->Branch("pvqual",pvQual_,"ndf/F:nchi2/F:sumpt/F:wsumpt/F:sumw/F");
  tree_->Branch("pvntk",&pvNTk_,"ntk/I");
  tree_->Branch("tknpix",&tkQualI_[0][0],"npix[ntk]/s");
  tree_->Branch("tknstrip",&tkQualI_[1][0],"nstrip[ntk]/s");
  tree_->Branch("tknstrip",&tkQualI_[1][0],"nstrip[ntk]/s");
  tree_->Branch("tklost",&tkQualI_[2][0],"lost[ntk]/s");
  tree_->Branch("tkndf",&tkQualF_[0][0],"ndf[ntk]/i");
  tree_->Branch("tknchi2",&tkQualF_[1][0],"nchi2[ntk]/F");
  tree_->Branch("tkwgt",&tkQualF_[2][0],"wgt[ntk]/F");
  tree_->Branch("tkd0",&tkPar_[0][0],"d0[ntk]/F");
  tree_->Branch("tkdz",&tkPar_[1][0],"dz[ntk]/F");
  tree_->Branch("tkpt",&tkPar_[2][0],"pt[ntk]/F");
  tree_->Branch("tketa",&tkPar_[3][0],"eta[ntk]/F");
  tree_->Branch("tkphi",&tkPar_[4][0],"phi[ntk]/F");
  tree_->Branch("tksigd0",&tkErr_[0][0],"sigd0[ntk]/F");
  tree_->Branch("tksigdz",&tkErr_[1][0],"sigdz[ntk]/F");
  tree_->Branch("tksigpt",&tkErr_[2][0],"sigpt[ntk]/F");
  tree_->Branch("tksigeta",&tkErr_[3][0],"sigeta[ntk]/F");
  tree_->Branch("tksigphi",&tkErr_[4][0],"sigphi[ntk]/F");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PVAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PVAnalyzer);
