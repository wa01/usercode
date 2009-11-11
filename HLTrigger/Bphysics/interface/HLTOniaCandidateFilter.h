#ifndef HLTOniaCandidateFilter_h_
#define HLTOniaCandidateFilter_h_
// -*- C++ -*-
//
// Package:    HLTOniaCandidateFilter
// Class:      HLTOniaCandidateFilter
// 
/**\class HLTOniaCandidateFilter HLTOniaCandidateFilter.cc HLTrigger/HLTOniaCandidateFilter/src/HLTOniaCandidateFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wolfgang ADAM
//         Created:  Wed Oct 14 17:57:21 CEST 2009
// $Id: HLTOniaCandidateFilter.h,v 1.1.1.1 2009/10/26 11:42:53 adamwo Exp $
//
//


// system include files
#include <memory>

// user include files
#include "HLTrigger/HLTcore/interface/HLTFilter.h"


//
// class declaration
//

class HLTOniaCandidateFilter : public HLTFilter {
public:
  explicit HLTOniaCandidateFilter(const edm::ParameterSet&);
  ~HLTOniaCandidateFilter();

private:
//   virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
//   virtual void endJob() ;
      
private:
  edm::InputTag candTag_;
  bool saveTag_;
  double minMass_;
  double maxMass_;
  bool checkCharge_;
};

#endif
