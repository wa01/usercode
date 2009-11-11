#ifndef QuarkoniaTrackSelector_h_
#define QuarkoniaTrackSelector_h_
// -*- C++ -*-
//
// Package:    QuarkoniaTrackSelector
// Class:      QuarkoniaTrackSelector
// 
/**\class QuarkoniaTrackSelector QuarkoniaTrackSelector.cc HLTrigger/QuarkoniaTrackSelector/src/QuarkoniaTrackSelector.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wolfgang ADAM
//         Created:  Wed Oct 14 17:57:21 CEST 2009
// $Id: QuarkoniaTrackSelector.h,v 1.1.1.1 2009/10/26 11:42:53 adamwo Exp $
//
//

// user include files
#include "FWCore/Framework/interface/EDProducer.h"

// system include files
#include <memory>
#include <vector>



//
// class declaration
//

class QuarkoniaTrackSelector : public edm::EDProducer {
public:
  explicit QuarkoniaTrackSelector(const edm::ParameterSet&);
  ~QuarkoniaTrackSelector();

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
      
private:
  edm::InputTag muonTag_;
  edm::InputTag trackTag_;
  std::vector<double> minMasses_;
  std::vector<double> maxMasses_;
  bool checkCharge_;
};

#endif
