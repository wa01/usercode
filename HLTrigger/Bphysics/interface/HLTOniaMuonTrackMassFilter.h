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
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "HLTrigger/HLTcore/interface/HLTFilter.h"


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
  edm::InputTag muonTag_;
  edm::InputTag trackTag_;
  bool saveTag_;
  double minMass_;
  double maxMass_;
  bool checkCharge_;
  double minTrackPt_;
  double minMuonPt_;
};

#endif
