#ifndef HLTCompositeMuonTrackProducer_h_
#define HLTCompositeMuonTrackProducer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



class HLTCompositeMuonTrackProducer : public edm::EDProducer {
public:
  explicit HLTCompositeMuonTrackProducer(const edm::ParameterSet&);
  ~HLTCompositeMuonTrackProducer();

private:
//   virtual void beginJob() ;
  virtual void produce (edm::Event&, const edm::EventSetup&);
//   virtual void endJob() ;
      
private:
  edm::InputTag muonTag_;
  edm::InputTag trackTag_;
  double minMass_;
  double maxMass_;
  bool checkCharge_;
};

#endif
