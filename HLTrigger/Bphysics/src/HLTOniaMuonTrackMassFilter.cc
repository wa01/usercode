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
  std::vector<reco::RecoChargedCandidateRef> prevCandRefs;
  prevCandHandle->getObjects(trigger::TriggerMuon,prevCandRefs);
  std::vector<reco::RecoChargedCandidateRef> prevTrackRefs;
  prevCandHandle->getObjects(trigger::TriggerTrack,prevTrackRefs);
  LogDebug("HLTOniaMuonTrackMassFilter") << "#previous track refs = " << prevTrackRefs.size();
  //
  // access to muons and selection according to configuration
  //   if the previous candidates are taken from a muon+track
  //   quarkonia filter we rely on the fact that only the muons 
  //   are flagged as TriggerMuon 
  //
  std::vector<reco::RecoChargedCandidateRef> selectedMuonRefs;
  selectedMuonRefs.reserve(muonHandle->size());
  for ( unsigned int i=0; i<muonHandle->size(); ++i ) {
    // Ref
    reco::RecoChargedCandidateRef muonRef(muonHandle,i);
    LogDebug("HLTOniaMuonTrackMassFilter") << "Checking muon with q / pt / p / eta = "
					   << muonRef->charge() << " " << muonRef->pt() << " "
					   << muonRef->p() << " " << muonRef->eta() << " ";
    // passed previous filter?
    if ( find(prevCandRefs.begin(),prevCandRefs.end(),muonRef)==
	 prevCandRefs.end() )  continue;
    LogDebug("HLTOniaMuonTrackMassFilter") << "located in previous candidates";
    // keep muon
    selectedMuonRefs.push_back(muonRef);
  }
  //
  // access to tracks and selection according to configuration
  //
  std::vector<reco::RecoChargedCandidateRef> selectedTrackRefs;
  selectedTrackRefs.reserve(trackHandle->size());
  for ( unsigned int i=0; i<trackHandle->size(); ++i ) {
    // validity of REF
    reco::RecoChargedCandidateRef trackRef(trackHandle,i);
    const reco::RecoChargedCandidate& trackCand = *trackRef;
    LogDebug("HLTOniaMuonTrackMassFilter") << "Checking track with q / pt / p / eta = "
					   << trackCand.charge() << " " << trackCand.pt() << " "
					   << trackCand.p() << " " << trackCand.eta();
    // cuts on the momentum
    if ( trackCand.pt()<minTrackPt_ || trackCand.p()<minTrackP_ ||
	 fabs(trackCand.eta())>maxTrackEta_ )  continue;
    if ( trackCand.track().isNull() )  continue;
    // cuts on track quality
    const reco::Track& track = *trackCand.track();
    LogDebug("HLTOniaMuonTrackMassFilter") << "Checking track with dxy / dz / #hits / chi2 = "
					   << track.dxy(beamspot) << " "
					   << track.dz(beamspot) << " "
					   << track.numberOfValidHits() << " "
					   << track.normalizedChi2();
    if ( fabs(track.dxy(beamspot))>maxTrackDxy_ ||
	 fabs(track.dz(beamspot))>maxTrackDz_ ||
	 track.numberOfValidHits()<minTrackHits_ ||
	 track.normalizedChi2()>maxTrackNormChi2_ )  continue;
    // keep track
    selectedTrackRefs.push_back(trackRef);
  }
  //
  // combinations
  //
  unsigned int nDz(0);
  unsigned int nQ(0);
  unsigned int nComb(0);
  reco::Particle::LorentzVector p4Muon;
  reco::Particle::LorentzVector p4JPsi;
  for ( unsigned int im=0; im<selectedMuonRefs.size(); ++im ) {
    const reco::RecoChargedCandidate& muon = *selectedMuonRefs[im];
    int qMuon = muon.charge();
    p4Muon = muon.p4();
    for ( unsigned int it=0; it<selectedTrackRefs.size(); ++it ) {
      const reco::RecoChargedCandidate& track = *selectedTrackRefs[it];
      LogDebug("HLTOniaMuonTrackMassFilter") << "combination with dz / q / mass = "
					     << muon.track()->dz(beamspot)-track.track()->dz(beamspot) << " "
					     << track.charge()+qMuon << " "
					     << (p4Muon+track.p4()).mass();
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
	  break;
	}
      }
    }
  }


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


//define this as a plug-in
DEFINE_FWK_MODULE(HLTOniaMuonTrackMassFilter);
