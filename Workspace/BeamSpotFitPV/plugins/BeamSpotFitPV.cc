#include "Workspace/BeamSpotFitPV/plugins/BeamSpotFitPV.h"
#include "Workspace/BeamSpotFitPV/interface/FcnBeamSpotFitPV.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TFitterMinuit.h"
#include "Minuit2/FCNBase.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include <iostream>

//
// constructor - setting parameters
//
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


//
// selection / caching of vertices
//
void
BeamSpotFitPV::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // offline beamspot
  //
  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel(beamspotTag_,bsHandle);
  if ( !bsHandle.isValid() )  return;
  //
  // primary vertices
  //
  edm::Handle< std::vector<reco::Vertex> > vertexHandle;
  iEvent.getByLabel(vertexTag_,vertexHandle);
  if ( !vertexHandle.isValid() )  return;
  //
  // vertex selection
  //
  if ( vertexHandle->size()<1 )  return;
  const reco::Vertex& pv = (*vertexHandle)[0];
  if ( pv.isFake() || pv.tracksSize()==0 )  return;
  if ( pv.ndof()<minVtxNdf_ || (pv.ndof()+3.)/pv.tracksSize()<2*minVtxWgt_ )  return;
  if ( fabs(pv.z()-(*bsHandle).z0())>maxVtxZ_ ||
       ((pv.x()-(*bsHandle).x0())*(pv.x()-(*bsHandle).x0())+
	(pv.y()-(*bsHandle).y0())*(pv.y()-(*bsHandle).y0()))>maxVtxR_*maxVtxR_ )  return;
  //
  // check, if vertex cache is full. if so, fit and reset
  //
  if ( pvStore_.size()>=maxNrVertices_ )  fitBeamspot();
  // keep track of first and last event
  if ( pvStore_.empty() ) {
    firstEvent_ = lastEvent_ = iEvent.id();
  }
  else {
    if ( iEvent.id()<firstEvent_ )  firstEvent_ = iEvent.id();
    if ( iEvent.id()>lastEvent_ )  lastEvent_ = iEvent.id();
  }
  //
  // store vertex data
  //
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
  //
  // if there are any remaining vertices: fit
  //
  if ( pvStore_.size()>minNrVertices_ )  fitBeamspot();
  //
  // sort list of luminosity blocks (just for safety)
  //
  std::sort(luminosityBins_.begin(),luminosityBins_.end());
  //
  // create histograms and graphs
  //
  edm::Service<TFileService> fs;
  //
  // histogram of number of primary vertices / luminosity block
  //
  unsigned int nls = luminosityBins_.size();
  TH1* h_count = fs->make<TH1F>("pvcounts","Nr. of selected primary vertices",
				nls,0.,nls);
  // labels and bin contents
  char label[128];
  unsigned int ibin(0);
  TAxis* axis = h_count->GetXaxis();
  for ( std::vector<LSBin>::const_iterator i=luminosityBins_.begin();
	i!=luminosityBins_.end(); ++i ) {
    sprintf(label,"%d/%d",i->run,i->luminosityBlock);
    ++ibin;
    h_count->SetBinContent(ibin,i->pvCount);
    axis->SetBinLabel(ibin,label);
  }
  //
  // definition of graphs of all fitted parameters
  //
  std::vector<TGraphErrors*> graphs;
  for ( unsigned int i=0; i<NFITPAR; ++i ) 
    graphs.push_back(fs->make<TGraphErrors>());
  graphs[0]->SetName("x");
  graphs[1]->SetName("y");
  graphs[2]->SetName("z");
  graphs[3]->SetName("ex");
  graphs[4]->SetName("corrxy");
  graphs[5]->SetName("ey");
  graphs[6]->SetName("dxdz");
  graphs[7]->SetName("dydz");
  graphs[8]->SetName("ez");
  graphs[9]->SetName("scale");

  //
  // filling of graphs (x-axis corresponds to the bins
  // in the histogram above)
  //
  LSBin resbin;
  unsigned int np(0);
  for ( unsigned int i=0; i<fitResults_.size(); ++i ) {
    const FitResult& fitResult = fitResults_[i];
    //
    // find bin range by looking up event ids in the list of luminosity blocks
    // 
    resbin.run = fitResult.firstEvent.run();
    resbin.luminosityBlock = fitResult.firstEvent.luminosityBlock();
    std::vector<LSBin>::const_iterator ifirst = 
      std::find(luminosityBins_.begin(),luminosityBins_.end(),resbin);
    resbin.run = fitResult.lastEvent.run();
    resbin.luminosityBlock = fitResult.lastEvent.luminosityBlock();
    std::vector<LSBin>::const_iterator ilast = 
      std::find(luminosityBins_.begin(),luminosityBins_.end(),resbin);
    if ( ifirst==luminosityBins_.end() || ilast==luminosityBins_.end() ) {
      std::cout << "Did not find luminosity bin for result!!" << std::endl;
      continue;
    }
    unsigned int ibfirst = ifirst - luminosityBins_.begin() + 1;
    unsigned int iblast = ilast - luminosityBins_.begin() + 1;
    //
    // store values
    //
    for ( unsigned int j=0; j<NFITPAR; ++j ) {
      graphs[j]->SetPoint(np,(ibfirst+iblast)/2.,fitResult.values[j]);
      graphs[j]->SetPointError(np,(iblast-ibfirst)/2.,fitResult.errors[j]);
    }
    ++np;
  }
  //
  // mandatory "Write" for TGraphErrors
  //
  for ( unsigned int i=0; i<NFITPAR; ++i )  graphs[i]->Write();

}

void
BeamSpotFitPV::beginRun (edm::Run const& run, edm::EventSetup const& setup)
{
}

void
BeamSpotFitPV::endRun (edm::Run const& run, edm::EventSetup const& setup)
{
  //
  // clear cache at end of run
  //
  pvStore_.clear();
}

void
BeamSpotFitPV::beginLuminosityBlock (edm::LuminosityBlock const& ls, edm::EventSetup const& setup)
{
  //
  // store cache size at the start of the luminosity block
  //
  pvCountAtLS_ = pvStore_.size();
}

void
BeamSpotFitPV::endLuminosityBlock (edm::LuminosityBlock const& ls, edm::EventSetup const& setup)
{
  //
  // store id of the luminosity block if any vertex was selected
  //
  if ( pvStore_.size()>pvCountAtLS_ ) {
    LSBin newBin;
    newBin.run = ls.run();
    newBin.luminosityBlock = ls.luminosityBlock();
    newBin.pvCount = pvStore_.size() - pvCountAtLS_;
    luminosityBins_.push_back(newBin);
  }
  //
  // fit results of the luminosity block, if minimum nr. of vertices was reached
  //
  if ( pvStore_.size()>minNrVertices_ )  fitBeamspot();
}


void
BeamSpotFitPV::fitBeamspot ()
{
  //
  // LL function and fitter
  //
  FcnBeamSpotFitPV* fcn = new FcnBeamSpotFitPV(pvStore_);
  TFitterMinuit* minuitx = new TFitterMinuit();
  minuitx->SetMinuitFCN(fcn); 
  //
  // fit parameters: positions, widths, x-y correlations, tilts in xz and yz
  //
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
  //
  // first iteration without correlations
  //
  minuitx->FixParameter(4);
  minuitx->FixParameter(6);
  minuitx->FixParameter(7);
  minuitx->FixParameter(9);
  minuitx->SetMaxIterations(100);
  minuitx->SetPrintLevel(3);
  minuitx->CreateMinimizer();
  minuitx->Minimize();
  //
  // refit with harder selection on vertices
  //
  fcn->setLimits(minuitx->GetParameter(0)-sigmaCut_*minuitx->GetParameter(3),
		 minuitx->GetParameter(0)+sigmaCut_*minuitx->GetParameter(3),
		 minuitx->GetParameter(1)-sigmaCut_*minuitx->GetParameter(5),
		 minuitx->GetParameter(1)+sigmaCut_*minuitx->GetParameter(5),
		 minuitx->GetParameter(2)-sigmaCut_*minuitx->GetParameter(8),
		 minuitx->GetParameter(2)+sigmaCut_*minuitx->GetParameter(8));
  minuitx->Minimize();
  //
  // refit with correlations
  //
  minuitx->ReleaseParameter(4);
  minuitx->ReleaseParameter(6);
  minuitx->ReleaseParameter(7);
  minuitx->Minimize();
  // refit with floating scale factor
//   minuitx->ReleaseParameter(9);
//   minuitx->Minimize();

  FitResult result;
  result.firstEvent = firstEvent_;
  result.lastEvent = lastEvent_;
  for ( unsigned int i=0; i<NFITPAR; ++i ) {
    result.values[i] = minuitx->GetParameter(i);
    result.errors[i] = minuitx->GetParError(i);
  }
  fitResults_.push_back(result);

  std::cout << "Fitted beamspot for " << fcn->nrOfVerticesUsed()
	    << " from run / LS " 
	    << firstEvent_.run() << " / "
	    << firstEvent_.luminosityBlock()
	    << " to run / LS "
	    << lastEvent_.run() << " / "
	    << lastEvent_.luminosityBlock() << std::endl;

  pvStore_.clear();

  delete minuitx;
//   delete fcn;

}

//define this as a plug-in
DEFINE_FWK_MODULE(BeamSpotFitPV);
