#include "Workspace/BeamSpotFitPV/interface/BeamSpotFitPV.h"
#include "Workspace/BeamSpotFitPV/interface/FcnBeamSpotFitPV.h"
#include "Minuit2/FCNBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TFitterMinuit.h"
#include <iostream>

BeamSpotFitPV::BeamSpotFitPV(const edm::ParameterSet& iConfig) :
  minNrVertices_(iConfig.getParameter<unsigned int>("minNrVerticesForFit")),
  maxNrVertices_(iConfig.getParameter<unsigned int>("vertexCacheSize")),
  beamspotTag_(iConfig.getParameter<edm::InputTag>("beamspotTag")),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")),
  minVtxNdf_(iConfig.getParameter<double>("minVertexNdf")),
  maxVtxNormChi2_(iConfig.getParameter<double>("maxVertexNormChi2")),
  minVtxTracks_(iConfig.getParameter<unsigned int>("minVertexNTracks")),
  minVtxWgt_(iConfig.getParameter<double>("minVertexMeanWeight")),
  maxVtxR_(iConfig.getParameter<double>("maxVertexR")),
  maxVtxZ_(iConfig.getParameter<double>("maxVertexZ")),
  errorScale_(iConfig.getParameter<double>("errorScale")),
  sigmaCut_(iConfig.getParameter<double>("nSigmaCut")) {
  if ( minNrVertices_>maxNrVertices_ )  maxNrVertices_ = minNrVertices_;
}


BeamSpotFitPV::~BeamSpotFitPV()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BeamSpotFitPV::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel(beamspotTag_,bsHandle);
  if ( !bsHandle.isValid() )  return;

  edm::Handle< std::vector<reco::Vertex> > vertexHandle;
  iEvent.getByLabel(vertexTag_,vertexHandle);
  if ( !vertexHandle.isValid() )  return;
   
  if ( vertexHandle->size()<1 )  return;
  const reco::Vertex& pv = (*vertexHandle)[0];
  if ( pv.isFake() || pv.tracksSize()==0 )  return;
  if ( pv.ndof()<minVtxNdf_ || (pv.ndof()+3.)/pv.tracksSize()<2*minVtxWgt_ )  return;
  if ( fabs(pv.z()-(*bsHandle).z0())>maxVtxZ_ ||
       ((pv.x()-(*bsHandle).x0())*(pv.x()-(*bsHandle).x0())+
	(pv.y()-(*bsHandle).y0())*(pv.y()-(*bsHandle).y0()))>maxVtxR_*maxVtxR_ )  return;
  //
  // check, if vertex cache is full
  //
  if ( pvStore_.size()>=maxNrVertices_ )  fitBeamspot();
  //
  // store vertex
  //
  if ( pvStore_.empty() ) {
    firstEvents_.push_back(iEvent.id());
    lastEvents_.push_back(iEvent.id());
//     runs_.push_back(std::make_pair<iEvent.run(),iEvent.run()>);
//     lumiBlocks_.push_back(std::make_pair<iEvent.luminosityBlock(),iEvent.luminosityBlock()>);
  }
  else {
    if ( iEvent.id()<firstEvents_.back() )  firstEvents_.back() = iEvent.id();
    if ( iEvent.id()>lastEvents_.back() )  lastEvents_.back() = iEvent.id();
  }
  BeamSpotFitPVData pvData;

  pvData.position[0] = pv.x();
  pvData.position[1] = pv.y();
  pvData.position[2] = pv.z();
  pvData.posError[0] = pv.xError();
  pvData.posError[1] = pv.yError();
  pvData.posError[2] = pv.zError();
  pvData.posCorr[0] = pv.covariance(0,1)/pv.xError()/pv.yError();
  pvData.posCorr[1] = pv.covariance(0,2)/pv.xError()/pv.zError();
  pvData.posCorr[2] = pv.covariance(1,2)/pv.yError()/pv.zError();
  pvStore_.push_back(pvData);

}

// ------------ method called once each job just before starting event loop  ------------
void 
BeamSpotFitPV::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BeamSpotFitPV::endJob() {
  if ( pvStore_.size()>minNrVertices_ )  fitBeamspot();
}

void
BeamSpotFitPV::beginRun (edm::Run const& run, edm::EventSetup const& setup)
{
  std::cout << "beginRun for run " << run.run() << std::endl;
}

void
BeamSpotFitPV::endRun (edm::Run const& run, edm::EventSetup const& setup)
{
  std::cout << "endRun for run " << run.run() << std::endl;
}

void
BeamSpotFitPV::beginLuminosityBlock (edm::LuminosityBlock const& ls, edm::EventSetup const& setup)
{
  std::cout << "beginLuminosityBlock for run / LS " << ls.run() 
	    << " / " << ls.luminosityBlock() << std::endl;
}

void
BeamSpotFitPV::endLuminosityBlock (edm::LuminosityBlock const& ls, edm::EventSetup const& setup)
{
  std::cout << "endLuminosityBlock for run / LS " << ls.run() 
	    << " / " << ls.luminosityBlock() 
	    << " pvStore size = " << pvStore_.size() << std::endl;
  if ( pvStore_.size()>minNrVertices_ )  fitBeamspot();
}


void
BeamSpotFitPV::fitBeamspot ()
{
  FcnBeamSpotFitPV* fcn = new FcnBeamSpotFitPV(pvStore_);
  TFitterMinuit* minuitx = new TFitterMinuit();
  minuitx->SetMinuitFCN(fcn); 
  minuitx->SetParameter(0,"x",0.,0.02,-10.,10.);
  minuitx->SetParameter(1,"y",0.,0.02,-10.,10.);
  minuitx->SetParameter(2,"z",0.,0.20,-30.,30.);
  minuitx->SetParameter(3,"ex",0.015,0.01,0.,10.);
  minuitx->SetParameter(4,"corrxy",0.,0.02,-1.,1.);
  minuitx->SetParameter(5,"ey",0.015,0.01,0.,10.);
  minuitx->SetParameter(6,"dxdz",0.,0.0002,-0.1,0.1);
  minuitx->SetParameter(7,"dydz",0.,0.0002,-0.1,0.1);
  minuitx->SetParameter(8,"ez",1.,0.1,0.,30.);
  minuitx->SetParameter(9,"scale",0.9,0.1,0.5,2.);
  minuitx->FixParameter(4);
  minuitx->FixParameter(6);
  minuitx->FixParameter(7);
  minuitx->FixParameter(9);
  minuitx->SetMaxIterations(100);
  minuitx->SetPrintLevel(3);
  minuitx->CreateMinimizer();
  minuitx->Minimize();
  const double sigmaCut3D_ = 5.;
  fcn->setLimits(minuitx->GetParameter(0)-sigmaCut3D_*minuitx->GetParameter(3),
		 minuitx->GetParameter(0)+sigmaCut3D_*minuitx->GetParameter(3),
		 minuitx->GetParameter(1)-sigmaCut3D_*minuitx->GetParameter(5),
		 minuitx->GetParameter(1)+sigmaCut3D_*minuitx->GetParameter(5),
		 minuitx->GetParameter(2)-sigmaCut3D_*minuitx->GetParameter(8),
		 minuitx->GetParameter(2)+sigmaCut3D_*minuitx->GetParameter(8));
  minuitx->Minimize();
  minuitx->ReleaseParameter(4);
  minuitx->ReleaseParameter(6);
  minuitx->ReleaseParameter(7);
  minuitx->Minimize();
//   minuitx->ReleaseParameter(9);
//   minuitx->Minimize();

  std::cout << "Fitted beamspot for " << fcn->nrOfVerticesUsed()
	    << " from run / LS " 
	    << firstEvents_.back().run() << " / "
	    << firstEvents_.back().luminosityBlock()
	    << " to run / LS "
	    << lastEvents_.back().run() << " / "
	    << lastEvents_.back().luminosityBlock() << std::endl;

  pvStore_.clear();

  delete minuitx;
//   delete fcn;
}

//define this as a plug-in
DEFINE_FWK_MODULE(BeamSpotFitPV);
