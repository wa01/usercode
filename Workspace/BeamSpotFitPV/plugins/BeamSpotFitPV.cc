#include "Workspace/BeamSpotFitPV/plugins/BeamSpotFitPV.h"
#include "Workspace/BeamSpotFitPV/interface/FcnBeamSpotFitPV.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TFitterMinuit.h"
#include "Minuit2/FCNBase.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TGraphErrors.h"
#include <iostream>
#include <map>
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
  sigmaCut_(iConfig.getParameter<double>("nSigmaCut")),
  produceHistograms_(iConfig.getParameter<bool>("histograms")),
  tFileService_(0),
  previousLuminosityBlock_(0) {
  if ( minNrVertices_>maxNrVertices_ )  maxNrVertices_ = minNrVertices_;
//   dynamicMinVtxNdf_ = minVtxNdf_;
  dynamicQualityCut_ = 1.e30;
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
//   if ( pv.isFake() || pv.tracksSize()==0 )  return;
//   if ( pv.ndof()<minVtxNdf_ || (pv.ndof()+3.)/pv.tracksSize()<2*minVtxWgt_ )  return;
//   if ( fabs(pv.z()-(*bsHandle).z0())>maxVtxZ_ ||
//        ((pv.x()-(*bsHandle).x0())*(pv.x()-(*bsHandle).x0())+
// 	(pv.y()-(*bsHandle).y0())*(pv.y()-(*bsHandle).y0()))>maxVtxR_*maxVtxR_ )  return;
  if ( !acceptVertex(pv,*bsHandle) )  return;
  //
  // check, if vertex cache is full. if so, fit and reset
  //
  if ( pvCache_.size()>=maxNrVertices_ ) {
    edm::LogInfo("BeamSpotFitPV") << "Cache exceeded max. size - forcing fit and reset";
    compressCache();
    // recheck vertex with new ndf cut
//     if ( pv.ndof()<dynamicMinVtxNdf_ )  return;
    if ( pvQuality(pv)>dynamicQualityCut_ ) return;
  }
  if ( pvCache_.size()>=maxNrVertices_ )  fitBeamspot();
  // keep track of first and last event
  if ( pvCache_.empty() ) {
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
//   pvData.ndf = pv.ndof();
  pvCache_.push_back(pvData);

}

// ------------ method called once each job just before starting event loop  ------------
void 
BeamSpotFitPV::beginJob()
{
  tFileService_ = new edm::Service<TFileService>();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BeamSpotFitPV::endJob() {

  delete tFileService_;
}

void
BeamSpotFitPV::beginRun (edm::Run const& run, edm::EventSetup const& setup)
{
  previousLuminosityBlock_ = 0;
}

void
BeamSpotFitPV::endRun (edm::Run const& run, edm::EventSetup const& setup)
{
  //
  // if there are any remaining vertices: fit
  //
  if ( pvCache_.size()>minNrVertices_ )  fitBeamspot();
  resetCache();
  //
  // write results
  //
  saveResults(run.run());
}

void
BeamSpotFitPV::beginLuminosityBlock (edm::LuminosityBlock const& ls, edm::EventSetup const& setup)
{
//   //
//   // check contiguity
//   //
//   if ( ls.luminosityBlock() < previousLuminosityBlock_ ) {
//     edm::LogWarning("BeamSpotFitPV")  << "luminosity blocks out of sequence";
//     pvCache_.clear();
//   }
  previousLuminosityBlock_ = ls.luminosityBlock();
  //
  // store cache size at the start of the luminosity block
  //
  pvCountAtLS_ = pvCache_.size();
}

void
BeamSpotFitPV::endLuminosityBlock (edm::LuminosityBlock const& ls, edm::EventSetup const& setup)
{
  //
  // store id of the luminosity block if any vertex was selected
  //
  if ( pvCache_.size()>pvCountAtLS_ ) {
    LSBin newBin;
    newBin.run = ls.run();
    newBin.luminosityBlock = ls.luminosityBlock();
    newBin.pvCount = pvCache_.size() - pvCountAtLS_;
    std::vector<LSBin>::iterator ibin = 
      find(luminosityBins_.begin(),luminosityBins_.end(),newBin);
    if ( ibin!=luminosityBins_.end() ) {
      edm::LogWarning("BeamSpotFitPV")  << "Found identical luminosity bin";
      ibin->pvCount += newBin.pvCount;
    }
    else {
      luminosityBins_.push_back(newBin);
    }
  }
  //
  // fit results of the luminosity block, if minimum nr. of vertices was reached
  //
  if ( pvCache_.size()>minNrVertices_ )  fitBeamspot();
  else {
    edm::LogInfo("BeamSpotFitPV") << "Insufficient nr. of vertices at endLuminosityBlock: " 
				  << " ( run " << ls.run()
				  << " ls " << ls.luminosityBlock() << " ): "
				  << pvCache_.size();
  }
}


bool
BeamSpotFitPV::fitBeamspot ()
{
  edm::LogInfo("BeamSpotFitPV")  << "Calling beamspot fit for " << pvCache_.size()
				 << " vertices from run / LS " 
				 << firstEvent_.run() << " / "
				 << firstEvent_.luminosityBlock()
				 << " to run / LS "
				 << lastEvent_.run() << " / "
				 << lastEvent_.luminosityBlock();
  //
  // LL function and fitter
  //
  FcnBeamSpotFitPV* fcn = new FcnBeamSpotFitPV(pvCache_);
  TFitterMinuit minuitx;
  minuitx.SetMinuitFCN(fcn); 
  //
  // fit parameters: positions, widths, x-y correlations, tilts in xz and yz
  //
  minuitx.SetParameter(0,"x",0.,0.02,-10.,10.);
  minuitx.SetParameter(1,"y",0.,0.02,-10.,10.);
  minuitx.SetParameter(2,"z",0.,0.20,-30.,30.);
  minuitx.SetParameter(3,"ex",0.015,0.01,0.,10.);
  minuitx.SetParameter(4,"corrxy",0.,0.02,-1.,1.);
  minuitx.SetParameter(5,"ey",0.015,0.01,0.,10.);
  minuitx.SetParameter(6,"dxdz",0.,0.0002,-0.1,0.1);
  minuitx.SetParameter(7,"dydz",0.,0.0002,-0.1,0.1);
  minuitx.SetParameter(8,"ez",1.,0.1,0.,30.);
  minuitx.SetParameter(9,"scale",0.9,0.1,0.5,2.);
  //
  // first iteration without correlations
  //
  minuitx.FixParameter(4);
  minuitx.FixParameter(6);
  minuitx.FixParameter(7);
  minuitx.FixParameter(9);
  minuitx.SetMaxIterations(100);
//   minuitx.SetPrintLevel(3);
  minuitx.SetPrintLevel(0);

  int ierr(0);
  minuitx.CreateMinimizer();
  ierr = minuitx.Minimize();
  if ( ierr ) {
    edm::LogWarning("BeamSpotFitPV") << "minimization failed at 1st iteration";
    resetCache();
    return false;
  }
  //
  // refit with harder selection on vertices
  //
  fcn->setLimits(minuitx.GetParameter(0)-sigmaCut_*minuitx.GetParameter(3),
		 minuitx.GetParameter(0)+sigmaCut_*minuitx.GetParameter(3),
		 minuitx.GetParameter(1)-sigmaCut_*minuitx.GetParameter(5),
		 minuitx.GetParameter(1)+sigmaCut_*minuitx.GetParameter(5),
		 minuitx.GetParameter(2)-sigmaCut_*minuitx.GetParameter(8),
		 minuitx.GetParameter(2)+sigmaCut_*minuitx.GetParameter(8));
  ierr = minuitx.Minimize();
  if ( ierr ) {
    edm::LogWarning("BeamSpotFitPV") << "minimization failed at 2nd iteration";
    resetCache();
    return false;
  }
  //
  // refit with correlations
  //
  minuitx.ReleaseParameter(4);
  minuitx.ReleaseParameter(6);
  minuitx.ReleaseParameter(7);
  ierr = minuitx.Minimize();
  if ( ierr ) {
    edm::LogWarning("BeamSpotFitPV") << "minimization failed at 3rd iteration";
    resetCache();
    return false;
  }
  // refit with floating scale factor
//   minuitx.ReleaseParameter(9);
//   ierr = minuitx.Minimize();
//   if ( ierr ) {
//     edm::LogWarning("BeamSpotFitPV") << "minimization failed at 4th iteration";
//     return false;
//   }
  minuitx.PrintResults(0,0);

  FitResult result;
  result.firstEvent = firstEvent_;
  result.lastEvent = lastEvent_;
  for ( unsigned int i=0; i<NFITPAR; ++i ) {
    result.values[i] = minuitx.GetParameter(i);
    result.errors[i] = minuitx.GetParError(i);
  }
  fitResults_.push_back(result);
  edm::LogInfo("BeamSpotFitPV")  << "Fitted beamspot for " << fcn->nrOfVerticesUsed()
				 << " vertices from run / LS " 
				 << firstEvent_.run() << " / "
				 << firstEvent_.luminosityBlock()
				 << " to run / LS "
				 << lastEvent_.run() << " / "
				 << lastEvent_.luminosityBlock();
  //
  // observed and estimated vertex distributions / run
  //
  if ( produceHistograms_ ) {
    unsigned int run = firstEvent_.run();
    std::vector<unsigned int>::iterator irun = 
      find(processedRuns_.begin(),processedRuns_.end(),run);
    if ( irun==processedRuns_.end() ) {
      processedRuns_.push_back(run);
      pdfHistograms_.push_back(HistogramSet());
      char title[128];
      sprintf(title,"Run%d",run);
      runDirectories_.push_back((**tFileService_).mkdir(title));
      irun = processedRuns_.end() - 1;
    }
    unsigned int index = irun - processedRuns_.begin();
    HistogramSet& histograms = pdfHistograms_[index];
    TFileDirectory& runDir = runDirectories_[index];
    if ( histograms.chi2==0 ) {
      histograms.chi2 = runDir.make<TH1F>("chi2","PV chi2 probability",200,0,1);
      double xmin = minuitx.GetParameter(0) - 10*minuitx.GetParameter(3);
      double xmax = minuitx.GetParameter(0) + 10*minuitx.GetParameter(3);
      double ymin = minuitx.GetParameter(1) - 10*minuitx.GetParameter(5);
      double ymax = minuitx.GetParameter(1) + 10*minuitx.GetParameter(5);
      double zmin = minuitx.GetParameter(2) - 10*minuitx.GetParameter(8);
      double zmax = minuitx.GetParameter(2) + 10*minuitx.GetParameter(8);
      const int nbin(40);
      histograms.observed = runDir.make<TH3F>("PVobs","Measured PV positions",
					      nbin,xmin,xmax,nbin,ymin,ymax,nbin,zmin,zmax);
      histograms.estimated = runDir.make<TH3F>("PVest","Estimated PV positions",
					       nbin,xmin,xmax,nbin,ymin,ymax,nbin,zmin,zmax);
    }
    edm::LogInfo("BeamSpotFitPV") << "Filling pdf histograms";
    fcn->fillPdfs(histograms.observed,histograms.estimated,histograms.chi2,
		  result.values);
  }
//   pvCache_.clear();
  resetCache();

  return true;
}


void
BeamSpotFitPV::saveResults (unsigned int run)
{
  if ( luminosityBins_.empty() )  return;

  std::vector<unsigned int>::iterator irun = 
    find(processedRuns_.begin(),processedRuns_.end(),run);
  if ( irun==processedRuns_.end() ) {
    processedRuns_.push_back(run);
    pdfHistograms_.push_back(HistogramSet());
    char title[128];
    sprintf(title,"Run%d",run);
    runDirectories_.push_back((**tFileService_).mkdir(title));
    irun = processedRuns_.end() - 1;
  }
  unsigned int index = irun - processedRuns_.begin();
  TFileDirectory& runDir = runDirectories_[index];
  //
  // sort list of luminosity blocks (just for safety)
  //
  std::sort(luminosityBins_.begin(),luminosityBins_.end());
  //
  // histogram of number of primary vertices / luminosity block
  // 
  unsigned int nls = luminosityBins_.size();
  TH1* h_count = 
    runDir.make<TH1F>("pvcounts","Nr. of selected primary vertices",nls,0.,nls);
  //
  // labels and bin contents
  //
  char title[128];
  unsigned int ibin(0);
  TAxis* axis = h_count->GetXaxis();
  for ( std::vector<LSBin>::const_iterator i=luminosityBins_.begin();
	i!=luminosityBins_.end(); ++i ) {
    if ( i->run!=run ) {
      edm::LogWarning("BeamSpotFitPV")  << "Inconsistent run number in BeamSpotFitPV::saveResults";
    }
    sprintf(title,"%d",i->luminosityBlock);
    ++ibin;
    h_count->SetBinContent(ibin,i->pvCount);
    axis->SetBinLabel(ibin,title);
  }
  sprintf(title,"Luminosity block, run %d",run);
  axis->SetTitle(title);
  //
  // definition of graphs of all fitted parameters
  //
  std::vector<TGraphErrors*> graphs;
  for ( unsigned int i=0; i<NFITPAR; ++i )  
    graphs.push_back(runDir.make<TGraphErrors>());
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
      edm::LogWarning("BeamSpotFitPV")  << "Did not find luminosity bin for result!!";
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
  //
  // clear arrays
  //
  pvCountAtLS_ = 0;
  previousLuminosityBlock_ = 0;
  luminosityBins_.clear();
  fitResults_.clear();
}

bool
BeamSpotFitPV::acceptVertex (const reco::Vertex& pv,
			     const reco::BeamSpot& bs) const
{
  if ( pv.isFake() || pv.tracksSize()==0 )  return false;
  if ( pv.ndof()<minVtxNdf_ || (pv.ndof()+3.)/pv.tracksSize()<2*minVtxWgt_ )   return false;
  if ( fabs(pv.z()-bs.z0())>maxVtxZ_ ||
       ((pv.x()-bs.x0())*(pv.x()-bs.x0())+
	(pv.y()-bs.y0())*(pv.y()-bs.y0()))>maxVtxR_*maxVtxR_ )  return false;
  if ( pvQuality(pv)>dynamicQualityCut_ )  return false;
  return true;
}

void
BeamSpotFitPV::resetCache ()
{
//   dynamicMinVtxNdf_ = minVtxNdf_;
  dynamicQualityCut_ = 1.e30;
  pvCache_.clear();
}

void
BeamSpotFitPV::compressCache ()
{
  //
  // fill vertex ndfs
  //
//   std::cout << "pvQualities_ before " << pvQualities_.size();
  pvQualities_.resize(pvCache_.size());
//   std::cout << " and after " << pvQualities_.size() << " resize" << std::endl;
  for ( unsigned int i=0; i<pvCache_.size(); ++i )
//     pvQualities_[i] = pvCache_[i].ndf;
    pvQualities_[i] = pvQuality(pvCache_[i]);
  sort(pvQualities_.begin(),pvQualities_.end());
  //
  // set new cut to median
  //
//   std::cout << "ndfs " << pvQualities_.front() << " "
// 	    << pvQualities_[pvQualities_.size()/2] << " "
// 	    << pvQualities_.back() << std::endl;
//   dynamicMinVtxNdf_ = pvNdfs_[pvNdfs_.size()/2];
//   std::cout << "Setting dynamicMinVtxNdf_ to " << dynamicMinVtxNdf_ << std::endl;
  dynamicQualityCut_ = pvQualities_[pvQualities_.size()/2];
  //
  // remove all vertices failing the cut from the cache
  //   (to be moved to a more efficient memory management!)
  //
  unsigned int iwrite(0);
  for ( unsigned int i=0; i<pvCache_.size(); ++i ) {
//     if ( pvCache_[i].ndf<dynamicMinVtxNdf_ )  continue;
    if ( pvQuality(pvCache_[i])>dynamicQualityCut_ )  continue;
    if ( i!=iwrite )  pvCache_[iwrite] = pvCache_[i];
    ++iwrite;
  }
  pvCache_.resize(iwrite);
  edm::LogInfo("BeamSpotFitPV") << "Setting dynamic quality cut_ to " << dynamicQualityCut_ 
				<< " , new cache size = " << pvCache_.size();

}


double
BeamSpotFitPV::pvQuality (const reco::Vertex& pv) const
{
  //
  // determinant of the transverse part of the PV covariance matrix
  //
  return 
    pv.covariance(0,0)*pv.covariance(1,1)-
    pv.covariance(0,1)*pv.covariance(0,1);
}

double
BeamSpotFitPV::pvQuality (const BeamSpotFitPVData& pv) const
{
  //
  // determinant of the transverse part of the PV covariance matrix
  //
  double ex = pv.posError[0];
  double ey = pv.posError[1];
  return ex*ex*ey*ey*(1-pv.posCorr[0]*pv.posCorr[0]);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BeamSpotFitPV);
