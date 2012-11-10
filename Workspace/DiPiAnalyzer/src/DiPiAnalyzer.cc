// -*- C++ -*-
//
// Package:    DiPiAnalyzer
// Class:      DiPiAnalyzer
// 
/**\class DiPiAnalyzer DiPiAnalyzer.cc Workspace/DiPiAnalyzer/src/DiPiAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wolfgang Adam,40 4-A28,,
//         Created:  Fri Nov  2 16:59:58 CET 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <vector>
#include <algorithm>
#include <iostream>
//
// class declaration
//

using namespace std;

class DiPiAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DiPiAnalyzer(const edm::ParameterSet&);
      ~DiPiAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  typedef reco::Candidate::PolarLorentzVector PolarLorentzVector;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  double selectPV (const reco::Vertex* selectedPV, const reco::Vertex* selectedPV2,
		   const reco::VertexCollection& pvs, const reco::Vertex& sv, 
		   InvariantMassFromVertex::LorentzVector p4Vtx) const;
  int sigTipSelection (const vector<const reco::Track*>& selectedTracks,
		       size_t iTkVtx1, size_t iTkVtx2, const reco::Vertex& pv) const;

private:
  edm::InputTag pvTag_;
  edm::InputTag tkTag_;

  float minRawMass_;
  float maxRawMass_;
  float maxVtxChi2_;
  float minVtxMass_;
  float maxVtxMass_;
  float sigTipLimit_;

  KalmanVertexFitter vtxFitter_;
  VertexDistanceXY vdistXY_;
  InvariantMassFromVertex massCalculator_;

  TTree* tree_;
  unsigned int npv_;
  float xyzPv_[3];
  unsigned int ntk_;
  float dzPv_;
  float dzPvPv2_;
  int nSigTip_;
  float chi2Vtx_;
  float flightDist_[2];
  float xyzVtx_[3];
  TLorentzVector pRaw_;
  float massVtx_[2];
  TLorentzVector pVtx_;
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
DiPiAnalyzer::DiPiAnalyzer(const edm::ParameterSet& iConfig) :
  pvTag_(iConfig.getParameter<edm::InputTag>("pvSrc")),
  tkTag_(iConfig.getParameter<edm::InputTag>("trackSrc")),
  minRawMass_(iConfig.getParameter<double>("minRawMass")),
  maxRawMass_(iConfig.getParameter<double>("maxRawMass")),
  maxVtxChi2_(iConfig.getParameter<double>("maxVtxChi2")),
  minVtxMass_(iConfig.getParameter<double>("minVtxMass")),
  maxVtxMass_(iConfig.getParameter<double>("maxVtxMass")),
  sigTipLimit_(iConfig.getParameter<double>("sigTipLimit")),
  vtxFitter_(true) {
   //now do what ever initialization is needed

}


DiPiAnalyzer::~DiPiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DiPiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   const vector<double> muMasses(2,0.105);

   Handle<reco::VertexCollection> vertexHandle;
   iEvent.getByLabel(pvTag_,vertexHandle);
   if ( vertexHandle->empty() )  return;
   npv_ = vertexHandle->size();

   Handle<reco::TrackCollection> trackHandle;
   iEvent.getByLabel(tkTag_,trackHandle);

   vector<const reco::Track*> selectedTracks;
   selectedTracks.reserve(trackHandle->size());
   for ( size_t i=0; i<trackHandle->size(); ++i ) {
     const reco::Track& track = (*trackHandle)[i];
     double pt = track.pt();
     double eta = fabs(track.eta());
     if ( eta<1.2 ) {
       if ( pt<4.5 )  continue;
     }
     else if ( eta<1.4 ) {
       if ( pt<3.5 )  continue;
     }
     else if ( eta<1.6 ) {
       if ( pt<3.0 )  continue;
     }
     else {
       continue;
     }
     if ( fabs(track.dxy())>3.0 || fabs(track.dz())>15.0 )  continue;
     if ( track.normalizedChi2()>1.8 )  continue;
     if ( track.found()<=10 )  continue;
     if ( track.hitPattern().numberOfValidPixelHits()<=1 )  continue;
     selectedTracks.push_back(&track);
   }
//    cout << "#selected tracks = " << selectedTracks.size() << endl;
   if ( selectedTracks.size()<2 )  return;
   ntk_ = selectedTracks.size();

  edm::ESHandle<TransientTrackBuilder> ttBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttBuilder);

  vector<reco::TransientTrack> tTracks;
  tTracks.reserve(selectedTracks.size());
  for ( size_t i=0; i<selectedTracks.size(); ++i )
    tTracks.push_back(ttBuilder->build(*selectedTracks[i]));

  PolarLorentzVector p1;
  PolarLorentzVector p2;
  PolarLorentzVector pDiPi;
  vector<reco::TransientTrack> vtxTTs(2);
  for ( size_t i1=0; i1<tTracks.size()-1; ++i1 ) {
    vtxTTs[0] = tTracks[i1];
    const reco::Track& tk1 = vtxTTs[0].track();
    p1.SetCoordinates(tk1.pt(),tk1.eta(),tk1.phi(),muMasses[0]);
    for ( size_t i2=i1+1; i2<tTracks.size(); ++i2 ) {
      vtxTTs[1] = tTracks[i2];
      const reco::Track& tk2 = vtxTTs[1].track();
      p2.SetCoordinates(tk2.pt(),tk2.eta(),tk2.phi(),muMasses[1]);
      pDiPi = p1 + p2;
      pRaw_.SetPtEtaPhiM(pDiPi.Pt(),pDiPi.Eta(),pDiPi.Phi(),pDiPi.M());
      if ( pRaw_.M()<minRawMass_ || pRaw_.M()>maxRawMass_ )  continue;
//       massRaw_ = pDiPi.M();

      CachingVertex<5> dipiCacheVertex = vtxFitter_.vertex(vtxTTs);
      if ( !dipiCacheVertex.isValid() )  continue;
      chi2Vtx_ = dipiCacheVertex.totalChiSquared();
      if ( chi2Vtx_>maxVtxChi2_ )  continue;

//       CachingVertex<5> dipiCacheVertex = dipiVertex;
      Measurement1D massWErr = massCalculator_.invariantMass(dipiCacheVertex,muMasses);
      massVtx_[0] = massWErr.value();
      massVtx_[1] = massWErr.error();
      if ( massVtx_[0]<minVtxMass_ || massVtx_[0]>maxVtxMass_ )  continue;

      InvariantMassFromVertex::LorentzVector p4Vtx = massCalculator_.p4(dipiCacheVertex,muMasses);
      pVtx_.SetXYZM(p4Vtx.X(),p4Vtx.Y(),p4Vtx.Z(),p4Vtx.M());

      TransientVertex dipiTVertex = dipiCacheVertex;
      reco::Vertex dipiVertex = dipiTVertex;
      xyzVtx_[0] = dipiVertex.x(); xyzVtx_[1] = dipiVertex.y(); xyzVtx_[2] = dipiVertex.z();

      TVector3 vtx;
      vtx.SetXYZ(dipiVertex.position().x(),dipiVertex.position().y(),0);

      const reco::Vertex* pv(0);
      const reco::Vertex* pv2(0);
      dzPv_ = selectPV(pv,pv2,*vertexHandle,dipiVertex,p4Vtx);
      if ( pv==0 )  continue;
      dzPvPv2_ = pv2 ? pv2->z()-pv->z() : 0.;
      xyzPv_[0] = pv->x(); xyzPv_[1] = pv->y(); xyzPv_[2] = pv->z();

      nSigTip_ = sigTipSelection(selectedTracks,i1,i2,*pv);

      TVector3 pvtx;
      pvtx.SetXYZ(pv->x(),pv->y(),0.);

      TVector3 pperp(pDiPi.px(), pDiPi.py(), 0);
      AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

      TVector3 vdiff = vtx - pvtx;
      double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
      Measurement1D distXY = vdistXY_.distance(dipiVertex,*pv);
      //double ctauPV = distXY.value()*cosAlpha*3.09688/pperp.Perp();
      flightDist_[0] = distXY.value()*cosAlpha * pDiPi.M()/pperp.Perp();
      GlobalError v1e = (dipiVertex).error();
      GlobalError v2e = pv->error();
      AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
      //double ctauErrPV = sqrt(vXYe.similarity(vpperp))*3.09688/(pperp.Perp2());
      flightDist_[1] = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*pDiPi.M()/(pperp.Perp2());

      tree_->Fill();
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
DiPiAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("dipi","dipi",1);
  tree_->Branch("npv",&npv_,"npv/i");
  tree_->Branch("ntk",&ntk_,"ntk/i");
//   tree_->Branch("massRaw",&massRaw_,"mRaw/F");
  tree_->Branch("chi2Vtx",&chi2Vtx_,"chi2Vtx/F");
  tree_->Branch("flightDist",flightDist_,"ct/F:ect/F");
  tree_->Branch("xyzPv",xyzPv_,"xPv/F:yPv/F:zPv/F");
  tree_->Branch("dzPv",&dzPv_,"dzPv/F");
  tree_->Branch("dzPvPv2",&dzPvPv2_,"dzPvPv2/F");
  tree_->Branch("nSigTip",&nSigTip_,"nSigTip/I");
  tree_->Branch("xyzVtx",xyzVtx_,"xVtx/F:yVtx/F:zVtx/F");
  tree_->Branch("pRaw","TLorentzVector",&pRaw_);
  tree_->Branch("pVtx","TLorentzVector",&pVtx_);
  tree_->Branch("massVtx",massVtx_,"mVtx/F:emVtx/F");
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiPiAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
DiPiAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DiPiAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DiPiAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DiPiAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiPiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

double
DiPiAnalyzer::selectPV (const reco::Vertex* selectedPV, const reco::Vertex* selectedPV2,
			const reco::VertexCollection& pvs, const reco::Vertex& sv, 
			InvariantMassFromVertex::LorentzVector p4Vtx) const
{
  selectedPV = selectedPV2 = 0;
  size_t iDzMin1(pvs.size());
  double dzMin1(1.e30);
  size_t iDzMin2(pvs.size());
  double dzMin2(1.e30);

  double pnx = p4Vtx.X()/p4Vtx.Pt();
  double pny = p4Vtx.Y()/p4Vtx.Pt();
  double pnz = p4Vtx.Z()/p4Vtx.Pt();
  for ( size_t i=0; i<pvs.size(); ++i ) {
//     double dxy = fabs(pnx*(pvs[i].y()-sv.y())-pny*(pvs[i].x()-sv.x()));
    double dz = (sv.z()-pvs[i].z()) - (pnx*(sv.x()-pvs[i].x())+pny*(sv.y()-pvs[i].y()))*pnz;
    if ( iDzMin1>=pvs.size() || fabs(dz)<fabs(dzMin1) ) {
      iDzMin2 = iDzMin1;
      dzMin2 = dzMin1;
      iDzMin1 = i;
      dzMin1 = dz;
    }
    else if ( iDzMin2>=pvs.size() || fabs(dz)<fabs(dzMin2) ) {
      iDzMin2 = i;
      dzMin2 = dz;
    }
  }
  if ( iDzMin1<pvs.size() ) {
    selectedPV = &pvs[iDzMin1];
    if ( iDzMin2<pvs.size() )  selectedPV2 = &pvs[iDzMin2];
  }
  std::cout << "PV sorting; z(SV) = " << sv.z() << std::endl;
  std::cout << "  PVPt2: " << 0 << " " 
	    << (sv.z()-p4Vtx.Z()*(pnx*(pvs[0].x()-sv.x())+pny*(pvs[0].y()-sv.y()))/p4Vtx.Pt())
	    << " " << pvs[0].z() << std::endl;
  if ( iDzMin1<pvs.size() ) 
    std::cout << "  PV1:    " << iDzMin1 << " " << dzMin1 << " " << pvs[iDzMin1].z() << std::endl;
  if ( iDzMin2<pvs.size() ) 
    std::cout << "  PV2:    " << iDzMin2 << " " << dzMin2 << " " << pvs[iDzMin2].z() << std::endl;
  return dzMin1;
}

int
DiPiAnalyzer::sigTipSelection (const vector<const reco::Track*>& selectedTracks,
			       size_t iTkVtx1, size_t iTkVtx2, const reco::Vertex& pv) const
{
  int result(0);
  for ( size_t i=0; i<selectedTracks.size(); ++i ) {
    if ( i==iTkVtx1 || i==iTkVtx2 )  continue;
    const reco::Track& track = *selectedTracks[i];
    if ( fabs(track.dz(pv.position()))>0.5 )  continue;
    double sigTip = fabs(track.dxy(pv.position()))/track.d0Error();
    if ( sigTip>sigTipLimit_ )  ++result;
  }
  return result;
}


//define this as a plug-in
DEFINE_FWK_MODULE(DiPiAnalyzer);
