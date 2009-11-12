#include "HLTrigger/Bphysics/interface/HLTOniaMuonTrackMassFilter.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <memory>
#include <iostream>
#include <sstream>


HLTOniaMuonTrackMassFilter::HLTOniaMuonTrackMassFilter(const edm::ParameterSet& iConfig) :
  beamspotTag_(iConfig.getParameter<edm::InputTag>("beamspot")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonCandidates")),
  trackTag_(iConfig.getParameter<edm::InputTag>("trackCandidates")),
  prevCandTag_(iConfig.getParameter<edm::InputTag>("previousCandidates")),
  saveTag_(iConfig.getParameter<bool>("saveTag")),
  minMasses_(iConfig.getParameter< std::vector<double> >("MinMasses")),
  maxMasses_(iConfig.getParameter< std::vector<double> >("MaxMasses")),
  checkCharge_(iConfig.getParameter<bool>("checkCharge")),
  minTrackPt_(iConfig.getParameter<double>("MinTrackPt")),
  minTrackP_(iConfig.getParameter<double>("MinTrackP")),
  maxTrackEta_(iConfig.getParameter<double>("MaxTrackEta")),
  maxTrackDxy_(iConfig.getParameter<double>("MaxTrackDxy")),
  maxTrackDz_(iConfig.getParameter<double>("MaxTrackDz")),
  minTrackHits_(iConfig.getParameter<int>("MinTrackHits")),
  maxTrackNormChi2_(iConfig.getParameter<double>("MaxTrackNormChi2")),
  maxDzMuonTrack_(iConfig.getParameter<double>("MaxDzMuonTrack"))
{
  //register your products
  produces<trigger::TriggerFilterObjectWithRefs>();
  //
  // verify mass windows
  //
  bool massesValid = minMasses_.size()==maxMasses_.size();
  if ( massesValid ) {
    for ( unsigned int i=0; i<minMasses_.size(); ++i ) {
      if ( minMasses_[i]<0 || maxMasses_[i]<0 || 
	   minMasses_[i]>maxMasses_[i] )  massesValid = false;
    }
  }
  if ( !massesValid ) {
    edm::LogError("HLTOniaMuonTrackMassFilter") << "Inconsistency in definition of mass windows, "
						<< "no event will pass the filter";
    minMasses_.clear();
    maxMasses_.clear();
  }

  std::ostringstream stream;
  stream << "instantiated with parameters\n";
  stream << "  beamspot = " << beamspotTag_ << "\n";
  stream << "  muonCandidates = " << muonTag_ << "\n";
  stream << "  trackCandidates = " << trackTag_ << "\n";
  stream << "  previousCandidates = " << prevCandTag_ << "\n";
  stream << "  saveTag = " << saveTag_ << "\n";
  stream << "  mass windows =";
  for ( size_t i=0; i<minMasses_.size(); ++i )  
    stream << " (" << minMasses_[i] << "," << maxMasses_[i] << ")";
  stream << "\n";
  stream << "  checkCharge = " << checkCharge_ << "\n";
  stream << "  MinTrackPt = " << minTrackPt_ << "\n";
  stream << "  MinTrackP = " << minTrackP_ << "\n";
  stream << "  MaxTrackEta = " << maxTrackEta_ << "\n";
  stream << "  MaxTrackDxy = " << maxTrackDxy_ << "\n";
  stream << "  MaxTrackDz = " << maxTrackDz_ << "\n";
  stream << "  MinTrackHits = " << minTrackHits_ << "\n";
  stream << "  MaxTrackNormChi2 = " << maxTrackNormChi2_ << "\n";
  stream << "  MaxDzMuonTrack = " << maxDzMuonTrack_ << "\n";
  LogDebug("HLTOniaMuonTrackMassFilter") << stream.str();

}

bool
HLTOniaMuonTrackMassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 // The filter object
  std::auto_ptr<trigger::TriggerFilterObjectWithRefs>
    filterproduct (new trigger::TriggerFilterObjectWithRefs(path(),module()));
  if ( saveTag_ ) {
    filterproduct->addCollectionTag(muonTag_);
    filterproduct->addCollectionTag(trackTag_);
  }
  //
  // Beamspot
  //
  edm::Handle<reco::BeamSpot> beamspotHandle;
  iEvent.getByLabel(beamspotTag_,beamspotHandle);
  reco::BeamSpot::Point beamspot(beamspotHandle->position());
  //
  // Muons
  //
  edm::Handle<reco::RecoChargedCandidateCollection> muonHandle;
  iEvent.getByLabel(muonTag_,muonHandle);
  //
  // Tracks
  //
  edm::Handle<reco::RecoChargedCandidateCollection> trackHandle;
  iEvent.getByLabel(trackTag_,trackHandle);
  //
  // Muons from previous filter
  //
  edm::Handle<trigger::TriggerFilterObjectWithRefs> prevCandHandle;
  iEvent.getByLabel(prevCandTag_,prevCandHandle);
  std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
  prevCandHandle->getObjects(trigger::TriggerMuon,prevMuonRefs);
  std::vector<reco::RecoChargedCandidateRef> prevTrackRefs;
  prevCandHandle->getObjects(trigger::TriggerTrack,prevTrackRefs);
  bool checkPrevTracks = prevTrackRefs.size()==prevMuonRefs.size();
//   LogDebug("HLTOniaMuonTrackMassFilter") << "#previous track refs = " << prevTrackRefs.size();
  //
  // access to muons and selection according to configuration
  //   if the previous candidates are taken from a muon+track
  //   quarkonia filter we rely on the fact that only the muons 
  //   are flagged as TriggerMuon 
  //
  std::vector<reco::RecoChargedCandidateRef> selectedMuonRefs;
  selectedMuonRefs.reserve(muonHandle->size());
  std::vector<size_t> prevMuonIndices;
  std::ostringstream stream1;
  for ( unsigned int i=0; i<muonHandle->size(); ++i ) {
    // Ref
    reco::RecoChargedCandidateRef muonRef(muonHandle,i);
//     LogDebug("HLTOniaMuonTrackMassFilter") << "Checking muon with q / pt / p / eta = "
// 					   << muonRef->charge() << " " << muonRef->pt() << " "
// 					   << muonRef->p() << " " << muonRef->eta();
    stream1 << "Checking muon with q / pt / p / eta = "
	    << muonRef->charge() << " " << muonRef->pt() << " "
	    << muonRef->p() << " " << muonRef->eta() << "\n";
    // passed previous filter?
    if ( find(prevMuonRefs.begin(),prevMuonRefs.end(),muonRef)==
	 prevMuonRefs.end() )  continue;
    prevMuonIndices.push_back(find(prevMuonRefs.begin(),prevMuonRefs.end(),muonRef)-
			      prevMuonRefs.begin());
//     LogDebug("HLTOniaMuonTrackMassFilter") << "located in previous candidates";
    // keep muon
    stream1 << "... accepted as #" << selectedMuonRefs.size() << "\n";
    selectedMuonRefs.push_back(muonRef);
  }
  LogDebug("HLTOniaMuonTrackMassFilter") << stream1.str();
  //
  // access to tracks and selection according to configuration
  //
  std::vector<reco::RecoChargedCandidateRef> selectedTrackRefs;
  selectedTrackRefs.reserve(trackHandle->size());
  std::ostringstream stream2;
  for ( unsigned int i=0; i<trackHandle->size(); ++i ) {
    // validity of REF
    reco::RecoChargedCandidateRef trackRef(trackHandle,i);
    const reco::RecoChargedCandidate& trackCand = *trackRef;
//     LogDebug("HLTOniaMuonTrackMassFilter") << "Checking track with q / pt / p / eta = "
// 					   << trackCand.charge() << " " << trackCand.pt() << " "
// 					   << trackCand.p() << " " << trackCand.eta();
    stream2 << "Checking track with q / pt / p / eta = "
	    << trackCand.charge() << " " << trackCand.pt() << " "
	    << trackCand.p() << " " << trackCand.eta() << "\n";
    // cuts on the momentum
    if ( trackCand.pt()<minTrackPt_ || trackCand.p()<minTrackP_ ||
	 fabs(trackCand.eta())>maxTrackEta_ )  continue;
    if ( trackCand.track().isNull() )  continue;
    // cuts on track quality
    const reco::Track& track = *trackCand.track();
//     LogDebug("HLTOniaMuonTrackMassFilter") << "Checking track with dxy / dz / #hits / chi2 = "
// 					   << track.dxy(beamspot) << " "
// 					   << track.dz(beamspot) << " "
// 					   << track.numberOfValidHits() << " "
// 					   << track.normalizedChi2();
    stream2 << "... with dxy / dz / #hits / chi2 = "
	    << track.dxy(beamspot) << " "
	    << track.dz(beamspot) << " "
	    << track.numberOfValidHits() << " "
	    << track.normalizedChi2();
    if ( fabs(track.dxy(beamspot))>maxTrackDxy_ ||
	 fabs(track.dz(beamspot))>maxTrackDz_ ||
	 track.numberOfValidHits()<minTrackHits_ ||
	 track.normalizedChi2()>maxTrackNormChi2_ )  continue;
    // keep track
    stream2 << "... accepted as #" << selectedTrackRefs.size() << "\n";
    selectedTrackRefs.push_back(trackRef);
  }
  LogDebug("HLTOniaMuonTrackMassFilter") << stream2.str();
  //
  // combinations
  //
  unsigned int nDz(0);
  unsigned int nQ(0);
  unsigned int nComb(0);
  reco::Particle::LorentzVector p4Muon;
  reco::Particle::LorentzVector p4JPsi;
  std::ostringstream stream3;
  for ( unsigned int im=0; im<selectedMuonRefs.size(); ++im ) {
    const reco::RecoChargedCandidate& muon = *selectedMuonRefs[im];
    int qMuon = muon.charge();
    p4Muon = muon.p4();
    for ( unsigned int it=0; it<selectedTrackRefs.size(); ++it ) {
      const reco::RecoChargedCandidate& track = *selectedTrackRefs[it];
//       LogDebug("HLTOniaMuonTrackMassFilter") << "combination with dz / q / mass = "
// 					     << muon.track()->dz(beamspot)-track.track()->dz(beamspot) << " "
// 					     << track.charge()+qMuon << " "
// 					     << (p4Muon+track.p4()).mass();
      stream3 << "Combination " << im << " / " << it << " with dz / q / mass = "
	      << muon.track()->dz(beamspot)-track.track()->dz(beamspot) << " "
	      << track.charge()+qMuon << " "
	      << (p4Muon+track.p4()).mass() << "\n";
      if ( fabs(muon.track()->dz(beamspot)-track.track()->dz(beamspot))>
	   maxDzMuonTrack_ )  continue;
      ++nDz;
      if ( checkCharge_ && track.charge()!=-qMuon )  continue;
      ++nQ;
      double mass = (p4Muon+track.p4()).mass();
      for ( unsigned int j=0; j<minMasses_.size(); ++j ) {
	if ( mass>minMasses_[j] && mass<maxMasses_[j] ) {
	  ++nComb;
	  filterproduct->addObject(trigger::TriggerMuon,selectedMuonRefs[im]);
	  filterproduct->addObject(trigger::TriggerTrack,selectedTrackRefs[it]);
	  stream3 << "... accepted\n";
      if ( checkPrevTracks ) {
	const reco::RecoChargedCandidateRef& prevTrack =
	  prevTrackRefs[prevMuonIndices[im]];
// 	LogDebug("HLTOniaMuonTrackMassFilter") << reco::deltaR(track.eta(),track.phi(),
// 							       prevTrack->eta(),
// 							       prevTrack->phi());
	stream3 << " deltaR for index " << prevMuonIndices[im] << " " 
		<< reco::deltaR(track.eta(),track.phi(),prevTrack->eta(),prevTrack->phi());
	pairMatched(prevMuonRefs,prevTrackRefs,selectedMuonRefs[im],selectedTrackRefs[it],stream3);
	stream3 << "\n";
	stream3 << " seedRef isNull " << track.track()->seedRef().isNull() << "\n";
      }
	  break;
	}
      }
    }
  }
  LogDebug("HLTOniaMuonTrackMassFilter") << stream3.str();


  if ( edm::isDebugEnabled() ) {
    std::ostringstream stream;
    stream << "Total number of combinations = " 
	   << selectedMuonRefs.size()*selectedTrackRefs.size() << " , after dz " << nDz
	   << " , after charge " << nQ << " , after mass " << nComb << std::endl;
    stream << "Found " << nComb << " jpsi candidates with # / mass / q / pt / eta" << std::endl;
    std::vector<reco::RecoChargedCandidateRef> muRefs;
    std::vector<reco::RecoChargedCandidateRef> tkRefs;
    filterproduct->getObjects(trigger::TriggerMuon,muRefs);
    filterproduct->getObjects(trigger::TriggerTrack,tkRefs);
    reco::Particle::LorentzVector p4Mu;
    reco::Particle::LorentzVector p4Tk;
    reco::Particle::LorentzVector p4JPsi;
    if ( muRefs.size()==tkRefs.size() ) {
      for ( unsigned int i=0; i<muRefs.size(); ++i ) {
	p4Mu = muRefs[i]->p4();
	p4Tk = tkRefs[i]->p4();
	p4JPsi = p4Mu + p4Tk;
	stream << "  " << i << " " 
	       << p4JPsi.M() << " "
	       << muRefs[i]->charge()+tkRefs[i]->charge() << " "
	       << p4JPsi.P() << " "
	       << p4JPsi.Eta() << "\n";
      }
      LogDebug("HLTOniaMuonTrackMassFilter") << stream.str();
    }
    else {
      LogDebug("HLTOniaMuonTrackMassFilter") << "different sizes for muon and track containers!!!";
    }
  }

  iEvent.put(filterproduct);

  return nComb>0;
}

bool
HLTOniaMuonTrackMassFilter::pairMatched (std::vector<reco::RecoChargedCandidateRef>& prevMuonRefs,
					 std::vector<reco::RecoChargedCandidateRef>& prevTrackRefs,
					 const reco::RecoChargedCandidateRef& muonRef,
					 const reco::RecoChargedCandidateRef& trackRef,
					 std::ostringstream& stream) const
{
  //
  // check only if references to tracks are available
  //
  if ( prevTrackRefs.empty() )  return true;
  //
  // validity
  //
  if ( prevMuonRefs.size()!=prevTrackRefs.size() )  return false;
  edm::RefToBase<TrajectorySeed> seedRef = trackRef->track()->seedRef();
  if ( seedRef.isNull() )  return false;
  //
  // comparison by hits of TrajectorySeed of the new track
  // with the previous candidate
  //
  TrajectorySeed::range seedHits = seedRef->recHits();
  trackingRecHit_iterator prevTrackHitEnd;
  trackingRecHit_iterator iprev;
  TrajectorySeed::const_iterator iseed;
  for ( size_t i=0; i<prevMuonRefs.size(); ++i ) {
    // identity of muon
    if ( prevMuonRefs[i]==muonRef )  continue;
    // validity of Ref to previous track candidate
    reco::TrackRef prevTrackRef = prevTrackRefs[i]->track();
    if ( prevTrackRef.isNull() )  continue;
    // same #hits
    if ( seedRef->nHits()!=prevTrackRef->recHitsSize() )  continue;
    // hit-by-hit comparison based on the sharesInput method
    iseed = seedHits.first;
    iprev = prevTrackRef->recHitsBegin();
    prevTrackHitEnd = prevTrackRef->recHitsEnd();
    bool identical(true);
    for ( ; iseed!=seedHits.second&&iprev!=prevTrackHitEnd; ++iseed,++iprev ) {
      if ( (*iseed).isValid()!=(**iprev).isValid() ||
	   !(*iseed).sharesInput(&**iprev,TrackingRecHit::all) ) {
	// terminate loop over hits on first mismatch
	identical = false;
	break;
      }
    }
    // if seed and previous track candidates are identical : return success
    if ( identical )  return true;
  }
  // no match found
  return false;
}

				       



//define this as a plug-in
DEFINE_FWK_MODULE(HLTOniaMuonTrackMassFilter);
