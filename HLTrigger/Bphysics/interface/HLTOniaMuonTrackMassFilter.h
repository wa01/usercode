#ifndef HLTOniaMuonTrackMassFilter_h_
#define HLTOniaMuonTrackMassFilter_h_
// -*- C++ -*-
//
// Package:    HLTOniaMuonTrackMassFilter
// Class:      HLTOniaMuonTrackMassFilter
// 
/**\class HLTOniaMuonTrackMassFilter HLTOniaMuonTrackMassFilter.cc HLTrigger/HLTOniaMuonTrackMassFilter/src/HLTOniaMuonTrackMassFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wolfgang ADAM
//         Created:  Wed Oct 14 17:57:21 CEST 2009
// $Id: HLTOniaMuonTrackMassFilter.h,v 1.1.1.1 2009/10/26 11:42:53 adamwo Exp $
//
//

// user include files
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

// system include files
#include <memory>
#include <vector>



//
// class declaration
//

class HLTOniaMuonTrackMassFilter : public HLTFilter {
public:
  explicit HLTOniaMuonTrackMassFilter(const edm::ParameterSet&);
  ~HLTOniaMuonTrackMassFilter();

private:
//   virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
//   virtual void endJob() ;
      
private:
  edm::InputTag beamspotTag_;
  edm::InputTag muonTag_;
  edm::InputTag trackTag_;
  edm::InputTag prevCandTag_;
  bool saveTag_;
  std::vector<double> minMasses_;
  std::vector<double> maxMasses_;
  bool checkCharge_;
  double minMuonPt_;
  double minMuonP_;
  double maxMuonEta_;
  double maxMuonDxy_;
  double maxMuonDz_;
  double minTrackPt_;
  double minTrackP_;
  double maxTrackEta_;
  double maxTrackDxy_;
  double maxTrackDz_;
  int minTrackHits_;
  double maxTrackNormChi2_;
  double maxDzMuonTrack_;
};

#endif
