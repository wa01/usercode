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
#include <set>

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

  edm::EventID firstEvent_;
  edm::EventID lastEvent_;

  class LSBin {
  public:
    LSBin () : run(0), luminosityBlock(0), pvCount(0) {}
    bool operator < (const LSBin& other) const {
      return run<other.run ||
	(run==other.run && luminosityBlock<other.luminosityBlock);
    }
    unsigned int run;
    unsigned int luminosityBlock;
    unsigned int pvCount;
  };

  unsigned int pvCountAtLS_;
  std::set<LSBin> luminosityBins_;

  static const unsigned int NFITPAR = 10;
  class FitResult {
  public:
    FitResult () : values(NFITPAR,0), errors(NFITPAR,0) {}
    std::vector<float> values;
    std::vector<float> errors;
    edm::EventID firstEvent;
    edm::EventID lastEvent;
  };
  std::vector<FitResult> fitResults;
};

#endif
