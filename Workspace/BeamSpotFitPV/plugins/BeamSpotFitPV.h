#ifndef BeamSpotFitPV_
#define BeamSpotFitPV_
/** \class BeamSpotFitPV
 *  EDAnalyzer fitting unfolded beam spot parameters from primary vertices.
 *  Fits are done at the end of a luminosity block (once a sufficient number
 *  of vertices has been accumulated) or when the vertex cache exceeds
 *  a certain size. The target function for the unbinned LL fit is
 *  defined in FcnBeamSpotFitPV.
 *   \author WA, 9/3/2010
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
//
// class declaration
//
#include "Workspace/BeamSpotFitPV/interface/BeamSpotFitPVData.h"
#include <algorithm>
#include <vector>

namespace reco {
  class Vertex;
  class BeamSpot;
}


class BeamSpotFitPV : public edm::EDAnalyzer {
public:
  explicit BeamSpotFitPV(const edm::ParameterSet&);
  ~BeamSpotFitPV() {}
  

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // the actual fit
  void fitBeamspot ();
  // saving results for one run
  void saveResults (unsigned int run);
  // PV quality cuts
  bool acceptVertex (const reco::Vertex& vertex,
		     const reco::BeamSpot& beamspot) const;
  // remove lowest-quality vertices from cache
  void compressCache();
  // clear cache and reset 
  void resetCache();
  // 
  double errorSquare (const reco::Vertex& vertex) const;
  double errorSquare (const BeamSpotFitPVData& vertex) const;
      // ----------member data ---------------------------

private:

  unsigned int minNrVertices_; //< min. nr. of vertices needed for a fit
  unsigned int maxNrVertices_; //< max. nr. of vertices that will be accumulated
  edm::InputTag beamspotTag_;  //< tag for the beamspot
  edm::InputTag vertexTag_;    //< tag for the primary vertices
  double minVtxNdf_;           //< min. ndof for vertices
  double maxVtxNormChi2_;      //< max. normalized chi2 for vertices
  unsigned minVtxTracks_;      //< min  nr. of tracks entering the vertex fit
  double minVtxWgt_;           //< min. average track weight in the vertex
  double maxVtxR_;             //< max. radial distance to beamspot
  double maxVtxZ_;             //< max. longitudinal distance to beamspot
  double errorScale_;          //< error scaling to be applied to the vertex
  double sigmaCut_;            //< vertex selection at 2nd iteration of the fit (nsigma from BS)

  double dynamicMinVtxNdf_;
  double dynamicMaxError2_;

  edm::Service<TFileService>* tFileService_;

  std::vector<BeamSpotFitPVData> pvStore_; //< cache for PV data

  edm::EventID firstEvent_;    //< event id for first PV in cache
  edm::EventID lastEvent_;     //< event id for last PV in cache

  std::vector<float> pvNdfs_;

  //
  // helper class keeping the identification of luminosity blocks with at least one PV
  //
  class LSBin {
  public:
    LSBin () : run(0), luminosityBlock(0), pvCount(0) {}
    bool operator < (const LSBin& other) const {
      return run<other.run ||
	(run==other.run && luminosityBlock<other.luminosityBlock);
    }
    bool operator == (const LSBin& other) const {
      return 
	run==other.run && luminosityBlock==other.luminosityBlock;
    }
    unsigned int run;
    unsigned int luminosityBlock;
    unsigned int pvCount;
  };

  unsigned int pvCountAtLS_;   //< cache size at the start of the current luminosity block  
  unsigned int previousLuminosityBlock_; //< for check of contiguity
  std::vector<LSBin> luminosityBins_;    //< list of all luminosity blocks with at least one PV

  std::vector<unsigned int> processedRuns_;

  static const unsigned int NFITPAR = 10;
  //
  // class storing values and uncertainties for one fit
  //
  class FitResult {
  public:
    FitResult () : values(NFITPAR,0), errors(NFITPAR,0) {}
    // fitted values
    std::vector<float> values;
    // fitted uncertainties
    std::vector<float> errors;
    // first and last event corresponding to the fit
    edm::EventID firstEvent;
    edm::EventID lastEvent;
  };
  std::vector<FitResult> fitResults_; //< list of fit results
};

#endif
