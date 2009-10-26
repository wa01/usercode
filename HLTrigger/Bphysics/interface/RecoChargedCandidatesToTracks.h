#ifndef RecoChargedCandidatesToTracks_h_
#define RecoChargedCandidatesToTracks_h_
// -*- C++ -*-
//
// Package:    RecoChargedCandidatesToTracks
// Class:      RecoChargedCandidatesToTracks
// 
/**\class RecoChargedCandidatesToTracks 
 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wolfgang ADAM
//         Created:  Thu Oct 15 17:29:11 CEST 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class decleration
//

class RecoChargedCandidatesToTracks : public edm::EDProducer {
public:
  explicit RecoChargedCandidatesToTracks(const edm::ParameterSet&);
  ~RecoChargedCandidatesToTracks();
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
private:
  edm::InputTag candTag_;
};

#endif

