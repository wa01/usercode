// -*- C++ -*-
//
// Package:    MuonSimHitCountFilter
// Class:      MuonSimHitCountFilter
// 
/**\class MuonSimHitCountFilter MuonSimHitCountFilter.cc Workspace/MuonSimHitCountFilter/src/MuonSimHitCountFilter.cc

 Description: requires min. nr. of PSimHits / DT or CSC "chamber" for each prompt muon SimTrack

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wolfgang Adam,40 4-A28,+41227671661,
//         Created:  Tue Nov 16 16:09:49 CET 2010
// $Id: MuonSimHitCountFilter.cc,v 1.3 2010/11/18 18:24:52 adamwo Exp $
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

//
// helper class defining a common id for a set of hits to be counted
//  (track id and "chamber" definition in DT or CSC)
//
struct SimHitClass {
  SimHitClass () : trackId_(-1), subdetId_(0), rawId_(0) {}
  SimHitClass (int trackId, DetId detId) :
    trackId_(trackId),subdetId_(detId.subdetId()),rawId_(0) {
    if ( detId.det()==DetId::Muon ) {
      if ( detId.subdetId()==MuonSubdetId::DT ) {
	// for DTs:  sum / chamber (remove other info from DetId)
	DTWireId tmpId(detId);
	rawId_ = DTWireId(tmpId.chamberId(),0,0,0).rawId();
      }
      else if ( detId.subdetId()==MuonSubdetId::CSC ) {
	// for CSCs: sum / (endcap,station,ring,chamber)
	//   (remove other info from DetId)
	CSCDetId tmpId(detId);
	rawId_ = CSCDetId(tmpId.endcap(),tmpId.station(),tmpId.ring(),tmpId.chamber());
      }
    }
    else {
      subdetId_ = 0;
    }
  }
  // comparison
  bool operator< (const SimHitClass& other) const {
    if ( trackId_!=other.trackId_ )  return trackId_<other.trackId_;
    if ( subdetId_!=other.subdetId_ )  return subdetId_<other.subdetId_;
    return rawId_<other.rawId_;
  }
  int trackId_;
  int subdetId_;
  uint32_t rawId_;
};

class MuonSimHitCountFilter : public edm::EDFilter {
public:
  explicit MuonSimHitCountFilter(const edm::ParameterSet&);
  ~MuonSimHitCountFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool acceptEvent(unsigned int trackId,
		   std::map<unsigned int, bool>& selectedTrackIds,
		   unsigned int& nrAcceptedTrackIds) const;
  
  // ----------member data ---------------------------
private:
  edm::InputTag simTrackTag_;
  std::vector<edm::InputTag> simHitTags_;
  std::vector<int> minHitsChamber_;
  std::vector<int> minHitsSubDet_;
  int minHitsTotal_;
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
  minHitsChamber_(iConfig.getParameter< std::vector<int> >("minHitsChamber")),
  minHitsSubDet_(iConfig.getParameter< std::vector<int> >("minHitsSubDet")),
  minHitsTotal_(iConfig.getParameter<int>("minHitsTotal")),
  wrongInputs_(false) {
   //now do what ever initialization is needed
  if ( (!minHitsChamber_.empty() && simHitTags_.size()!=minHitsChamber_.size()) ||
       (!minHitsSubDet_.empty() && simHitTags_.size()!=minHitsSubDet_.size()) ) {
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

//
bool
MuonSimHitCountFilter::acceptEvent(unsigned int trackId,
				   std::map<unsigned int, bool>& selectedTrackIds,
				   unsigned int& nrAcceptedTrackIds) const
{
  std::map<unsigned int, bool>::iterator iter = selectedTrackIds.find(trackId);
  if ( iter!=selectedTrackIds.end() && !iter->second ) {
    iter->second = true;
    ++nrAcceptedTrackIds;
  }
  if ( nrAcceptedTrackIds==selectedTrackIds.size() ) {
    std::cout << "event accepted" << std::endl;
    return true;
  }
  return false;
}

// ------------ method called on each new Event  ------------
bool
MuonSimHitCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if ( wrongInputs_ )  return false;

   using namespace edm;

   //
   // SimTrack collection 
   //

   // hit counts / track / chamber
   typedef SimHitClass SimClassType;
   std::map<SimClassType, int> chamberHitCounts;
   // selected tracks and hit counts / track
   std::map<unsigned int, bool> selectedTrackIds;
   unsigned int nrAcceptedTrackIds;
   std::map<unsigned int, int> trackHitCounts;
   std::map<unsigned int, int> allTrackHitCounts;

   Handle< std::vector<SimTrack> > simTrackHandle;
   iEvent.getByLabel(simTrackTag_,simTrackHandle);
   for ( unsigned int i=0; i<simTrackHandle->size(); ++i ) {
     const SimTrack& track = (*simTrackHandle)[i];
     // only primary tracks of selected particle types
     if ( track.vertIndex()==0 && 
	  (particleTypes_.empty() || 
	   particleTypes_.find(track.type())!=particleTypes_.end()) ) {
       selectedTrackIds[track.trackId()] = false;
     }
   }
   nrAcceptedTrackIds = 0;

   //
   // reset sum of hit counts
   //
   for ( std::map<unsigned int, bool>::iterator i=selectedTrackIds.begin(); 
	 i!=selectedTrackIds.end(); ++i ) {
     allTrackHitCounts[i->first] = 0;
   }

   //
   // loop over PSimHit collections (== sub dets)
   //
   for ( unsigned int i=0; i<simHitTags_.size(); ++i ) {

     //
     // reset hit counts
     //
     chamberHitCounts.clear();
     for ( std::map<unsigned int, bool>::iterator j=selectedTrackIds.begin(); 
	   j!=selectedTrackIds.end(); ++j ) {
       trackHitCounts[j->first] = 0;
     }

     //
     // PSimHits
     //
     Handle< std::vector<PSimHit> > simHitHandle;
     iEvent.getByLabel(simHitTags_[i],simHitHandle);
     std::cout << "Muon hits for " << simHitTags_[i] << std::endl;
     for ( unsigned j=0; j<simHitHandle->size(); ++j ) {
       const PSimHit& simHit = (*simHitHandle)[j];
       DetId detId(simHit.detUnitId());
       if ( detId.det()!=DetId::Muon )  continue;
       if ( detId.subdetId()!=MuonSubdetId::DT &&
	    detId.subdetId()!=MuonSubdetId::CSC )  continue;
       // only particle types in list
       if ( !particleTypes_.empty() && 
	    particleTypes_.find(simHit.particleType())==particleTypes_.end() )
	 continue;
       // only process types in list
       if ( !processTypes_.empty() && 
	    processTypes_.find(simHit.processType())==processTypes_.end() )
	 continue;
       // only for unaccepted track IDs in map
       unsigned int trackId = simHit.trackId();
       std::map<unsigned int, bool>::iterator itk = selectedTrackIds.find(trackId);
//        std::map<SimClassType, int>::iterator ihit = hitCounts.find(simHit.trackId());
//        if ( ihit!=hitCounts.end() )  ++ihit->second;
       if ( itk!=selectedTrackIds.end() && !itk->second ) {
	 SimHitClass simHitClass(simHit.trackId(),simHit.detUnitId());
	 ++chamberHitCounts[simHitClass];
	 if ( !minHitsChamber_.empty() && 
	      chamberHitCounts[simHitClass]>=minHitsChamber_[i] && 
	      acceptEvent(trackId,selectedTrackIds,nrAcceptedTrackIds) ) {
	   std::cout << "event accepted (chamber)" << std::endl;
// 	   return true;
	 }
       }
     }
     //
     // check hit counts for this sub-detector
     // (accept event, if at least one track has the min. nr. of hits / subdet)
     //
     std::cout << "Muon hits";;
     for ( std::map<SimClassType, int>::iterator itrk=chamberHitCounts.begin();
	   itrk!=chamberHitCounts.end(); ++itrk ) {
       unsigned int trackId = (*itrk).first.trackId_;
       int nHits = (*itrk).second;
       trackHitCounts[trackId] += nHits;
       if ( !minHitsSubDet_.empty() && nHits>=minHitsSubDet_[i] && 
	    acceptEvent(trackId,selectedTrackIds,nrAcceptedTrackIds) ) {
	 std::cout << "event accepted (SubDet)" << std::endl;
// 	 return true;
       }
       std::cout << " (" << (*itrk).first.trackId_ << "/" << (*itrk).first.subdetId_
		 << "/" << (*itrk).first.rawId_ << "," << nHits << ")";
     }
     std::cout << std::endl;
     std::cout << "Muon hits / track";;
     for ( std::map<unsigned int, int>::iterator itrk=trackHitCounts.begin();
	   itrk!=trackHitCounts.end(); ++itrk ) {
       allTrackHitCounts[(*itrk).first] += (*itrk).second;
       std::cout << " (" << (*itrk).first << "," << (*itrk).second << ")";
     }
     std::cout << std::endl;

   }

   //
   // check global hit counts (accept event, if at least one track 
   // has a total nr. of hits >= minimum)
   //
   std::cout << "All muon hits / track";
   for ( std::map<unsigned int, int>::iterator itrk=allTrackHitCounts.begin();
	 itrk!=allTrackHitCounts.end(); ++itrk ) {
     if ( minHitsTotal_>=0 && (*itrk).second>=minHitsTotal_ &&
	  acceptEvent((*itrk).first,selectedTrackIds,nrAcceptedTrackIds) ) {
       std::cout << "event accepted (total)" << std::endl;
//        return true;
     }
     std::cout << " (" << (*itrk).first << "," << (*itrk).second << ")";
   }
   std::cout << std::endl;

   return nrAcceptedTrackIds==selectedTrackIds.size();
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
