#ifndef BeamSpotFitPV_
#define BeamSpotFitPV_


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
//
// class declaration
//
#include "Workspace/BeamSpotFitPV/interface/BeamSpotFitPVData.h"
#include <algorithm>
#include <vector>

class BeamSpotFitPV : public edm::EDAnalyzer {
public:
  explicit BeamSpotFitPV(const edm::ParameterSet&);
  ~BeamSpotFitPV();
  
  

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void fitBeamspot ();
      // ----------member data ---------------------------
private:

  unsigned int minNrVertices_;
  unsigned int maxNrVertices_;
  edm::InputTag beamspotTag_;
  edm::InputTag vertexTag_;
  double minVtxNdf_;
  double maxVtxNormChi2_;
  unsigned minVtxTracks_;
  double minVtxWgt_;
  double maxVtxR_;
  double maxVtxZ_;
  double errorScale_;
  double sigmaCut_;

  std::vector<BeamSpotFitPVData> pvStore_;

  std::vector<edm::EventID> firstEvents_;
  std::vector<edm::EventID> lastEvents_;
//   std::vector< std::pair<unsigned int, unsigned int> > runs_;
//   std::vector< std::pair<unsigned int, unsigned int> > lumiBlocks_;
};

#endif
