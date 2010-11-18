// -*- C++ -*-
//
// Package:    MuonSimHitCountFilter
// Class:      MuonSimHitCountFilter
// 
/**\class MuonSimHitCountFilter MuonSimHitCountFilter.cc Workspace/MuonSimHitCountFilter/src/MuonSimHitCountFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wolfgang Adam,40 4-A28,+41227671661,
//         Created:  Tue Nov 16 16:09:49 CET 2010
// $Id: MuonSimHitCountFilter.cc,v 1.1 2010/11/18 16:25:53 adamwo Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <set>
#include <map>
#include <iostream>
//
// class declaration
//

struct SimHitClass {
  SimHitClass () : trackId_(-1), subdetId_(0) {}
  SimHitClass (int trackId, DetId detId) :
    trackId_(trackId),subdetId_(0) {
    if ( detId.det()==DetId::Muon ) {
      if ( detId.subdetId()==MuonSubdetId::DT ) {
	wireId_ = DTWireId(detId);
      }
      else if ( detId.subdetId()==MuonSubdetId::CSC ) {
	cscId_ = CSCDetId(detId);
      }
    }
  }
  bool operator< (const SimHitClass& other) const {
    if ( trackId_!=other.trackId_ )  return trackId_<other.trackId_;
    if ( subdetId_!=other.subdetId_ )  return subdetId_<other.subdetId_;
    if ( subdetId_==MuonSubdetId::DT ) {
      if ( wireId_.superLayer()!=other.wireId_.superLayer() )
	return wireId_.superLayer()<other.wireId_.superLayer();
      if ( wireId_.layer()!=other.wireId_.layer() )
	return wireId_.layer()<other.wireId_.layer();
      return wireId_.station()<other.wireId_.station();
    }
    else if ( subdetId_==MuonSubdetId::CSC ) {
      if ( cscId_.endcap()!=other.cscId_.endcap() )
	return cscId_.endcap()<other.cscId_.endcap();
      if ( cscId_.layer()!=other.cscId_.layer() )
	return cscId_.layer()<other.cscId_.layer();
      return cscId_.chamber()<other.cscId_.chamber();
    }
    return false;
  }
  int trackId_;
  int subdetId_;
  CSCDetId cscId_;
  DTWireId wireId_;
};

class MuonSimHitCountFilter : public edm::EDFilter {
public:
  explicit MuonSimHitCountFilter(const edm::ParameterSet&);
  ~MuonSimHitCountFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
private:
  edm::InputTag simTrackTag_;
  std::vector<edm::InputTag> simHitTags_;
  std::vector<int> minHits_;
  int minSumHits_;
  std::set<int> particleTypes_;
  std::set<int> processTypes_;

  bool wrongInputs_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonSimHitCountFilter::MuonSimHitCountFilter(const edm::ParameterSet& iConfig) :
  simTrackTag_(iConfig.getParameter<edm::InputTag>("simTracks")),
  simHitTags_(iConfig.getParameter< std::vector<edm::InputTag> >("simHits")),
  minHits_(iConfig.getParameter< std::vector<int> >("minHits")),
  minSumHits_(iConfig.getParameter<int>("minSumHits")),
  wrongInputs_(false) {
   //now do what ever initialization is needed
  if ( simHitTags_.size()!=minHits_.size() ) {
    edm::LogError("MuonSimHitCountFilter") << "Inconsistent inputs to MuonSimHitCountFilter";
    wrongInputs_ = true;
  }

  std::vector<int> types;
  types = iConfig.getParameter< std::vector<int> >("particleTypes");
  for ( unsigned int i=0; i<types.size(); ++i )  particleTypes_.insert(types[i]);
  types = iConfig.getParameter< std::vector<int> >("processTypes");
  for ( unsigned int i=0; i<types.size(); ++i )  processTypes_.insert(types[i]);
}


MuonSimHitCountFilter::~MuonSimHitCountFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MuonSimHitCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if ( wrongInputs_ )  return false;

   using namespace edm;

   //
   // SimTrack collection 
   //
   std::map<unsigned int, int> hitCounts;
   std::map<unsigned int, int> allHitCounts;
   std::vector<unsigned int> selectedTrackIds;
   Handle< std::vector<SimTrack> > simTrackHandle;
   iEvent.getByLabel(simTrackTag_,simTrackHandle);
   for ( unsigned int i=0; i<simTrackHandle->size(); ++i ) {
     const SimTrack& track = (*simTrackHandle)[i];
     // only primary tracks of selected particle types
     if ( track.vertIndex()==0 && 
	  (particleTypes_.empty() || 
	   particleTypes_.find(track.type())!=particleTypes_.end()) ) {
       selectedTrackIds.push_back(track.trackId());
     }
   }   

   //
   // reset sum of hit counts / selected track ID
   //
   for ( unsigned int i=0; i<selectedTrackIds.size(); ++i ) {
     allHitCounts[selectedTrackIds[i]] = 0;
   }

   for ( unsigned int i=0; i<simHitTags_.size(); ++i ) {

     //
     // reset hit counts / selected track ID
     //
     for ( unsigned int j=0; j<selectedTrackIds.size(); ++j ) {
       hitCounts[selectedTrackIds[j]] = 0;
     }

     //
     // PSimHits
     //
     Handle< std::vector<PSimHit> > simHitHandle;
     iEvent.getByLabel(simHitTags_[i],simHitHandle);
     for ( unsigned j=0; j<simHitHandle->size(); ++j ) {
       const PSimHit& simHit = (*simHitHandle)[j];
       DetId detId(simHit.detUnitId());
       std::cout << " detId " << detId.det() << " " << detId.subdetId();
       if ( detId.det()!=DetId::Muon )  continue;
       if ( detId.subdetId()==MuonSubdetId::DT ) {
	 DTWireId wireId(detId);
	 std::cout << " layer/superlayer/station " 
		   << wireId.layer() << " "
		   << wireId.superLayer() << " "
		   << wireId.station() << std::endl;
	 
       }
       else if ( detId.subdetId()==MuonSubdetId::CSC ) {
	 CSCDetId cscId(detId);
	 std::cout << " ec/layer/chamber " 
		   << cscId.endcap() << " "
		   << cscId.layer() << " "
		   << cscId.chamber() << std::endl;
       }
       else {
	 std::cout << std::endl;
	 continue;
       }
       // only particle types in list
       if ( !particleTypes_.empty() && 
	    particleTypes_.find(simHit.particleType())==particleTypes_.end() )
	 continue;
       // only process types in list
       if ( !processTypes_.empty() && 
	    processTypes_.find(simHit.processType())==processTypes_.end() )
	 continue;
       // only for track IDs in map
       std::map<unsigned int, int>::iterator ihit = hitCounts.find(simHit.trackId());
       if ( ihit!=hitCounts.end() )  ++ihit->second;
     }
//      std::cout << "Muon hits for " << simHitTags_[i];
     //
     // check hit counts for this sub-detector
     // (accept event, if at least one track has the min. nr. of hits / subdet)
     //
     for ( std::map<unsigned int, int>::iterator itrk=hitCounts.begin();
	   itrk!=hitCounts.end(); ++itrk ) {
       allHitCounts[(*itrk).first] += (*itrk).second;
       if ( (*itrk).second>minHits_[i] )  return true;
//         std::cout << " (" << (*itrk).first << "," << (*itrk).second << ")";
     }
      std::cout << std::endl;
     
   }

   //
   // check global hit counts (accept event, if at least one track 
   // has a total nr. of hits >= minimum)
   //
//    std::cout << "All muon hits ";
   for ( std::map<unsigned int, int>::iterator itrk=allHitCounts.begin();
	 itrk!=allHitCounts.end(); ++itrk ) {
     if ( (*itrk).second>minSumHits_ )  return true;
//      std::cout << " (" << (*itrk).first << "," << (*itrk).second << ")";
   }
//    std::cout << std::endl;

   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonSimHitCountFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonSimHitCountFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonSimHitCountFilter);
