#include "HLTrigger/Bphysics/interface/HLTOniaMuonTrackMassFilter.h"

// user include files
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

// system include files
#include <memory>
#include <iostream>
#include <sstream>

//
// constructors and destructor
//
HLTOniaMuonTrackMassFilter::HLTOniaMuonTrackMassFilter(const edm::ParameterSet& iConfig) :
  beamspotTag_(iConfig.getParameter<edm::InputTag>("beamspot")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonCandidates")),
  trackTag_(iConfig.getParameter<edm::InputTag>("trackCandidates")),
  prevCandTag_(iConfig.getParameter<edm::InputTag>("previousCandidates")),
  saveTag_(iConfig.getParameter<bool>("saveTag")),
  minMasses_(iConfig.getParameter< std::vector<double> >("MinMasses")),
  maxMasses_(iConfig.getParameter< std::vector<double> >("MaxMasses")),
  checkCharge_(iConfig.getParameter<bool>("checkCharge")),
  minMuonPt_(iConfig.getParameter<double>("MinMuonPt")),
  minMuonP_(iConfig.getParameter<double>("MinMuonP")),
  maxMuonEta_(iConfig.getParameter<double>("MaxMuonEta")),
  maxMuonDxy_(iConfig.getParameter<double>("MaxMuonDxy")),
  maxMuonDz_(iConfig.getParameter<double>("MaxMuonDz")),
  minTrackPt_(iConfig.getParameter<double>("MinTrackPt")),
  minTrackP_(iConfig.getParameter<double>("MinTrackP")),
  maxTrackEta_(iConfig.getParameter<double>("MaxTrackEta")),
  maxTrackDxy_(iConfig.getParameter<double>("MaxTrackDxy")),
  maxTrackDz_(iConfig.getParameter<double>("MaxTrackDz")),
  minTrackHits_(iConfig.getParameter<int>("MinTrackHits")),
  maxTrackNormChi2_(iConfig.getParameter<double>("MaxTrackNormChi2")),
  maxDzMuonTrack_(iConfig.getParameter<double>("MaxDzMuonTrack"))
{
  LogDebug("HLTTrackFilter") << "instantiated with parameters\n"
			     << "  muonTag  = " << muonTag_ << "\n"
			     << "  trackTag = " << trackTag_ << "\n"
			     << "  saveTag = " << saveTag_ << "\n"
// 			     << "  MinMass  = " << minMass_ << "\n"
// 			     << "  MaxMass  = " << maxMass_ << "\n"
			     << "  checkCharge  = " << checkCharge_;
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
}


HLTOniaMuonTrackMassFilter::~HLTOniaMuonTrackMassFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
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
  edm::Handle<trigger::TriggerFilterObjectWithRefs> muonHandle;
  iEvent.getByLabel(muonTag_,muonHandle);
  //
  // Tracks
  //
  edm::Handle<trigger::TriggerFilterObjectWithRefs> trackHandle;
  iEvent.getByLabel(trackTag_,trackHandle);
  //
  // Muons from previous filter
  //
  edm::Handle<trigger::TriggerFilterObjectWithRefs> prevCandHandle;
  iEvent.getByLabel(prevCandTag_,prevCandHandle);
  std::vector<reco::RecoChargedCandidateRef> prevCandRefs;
  prevCandHandle->getObjects(trigger::TriggerMuon,prevCandRefs);
  //
  // access to muons and selection according to configuration
  //   if the previous candidates are taken from a muon+track
  //   quarkonia filter we rely on the fact that only the muons 
  //   are flagged as TriggerMuon 
  //
  std::vector<reco::RecoChargedCandidateRef> muonRefs;
  muonHandle->getObjects(trigger::TriggerMuon,muonRefs);
  //
  std::vector<reco::RecoChargedCandidateRef> selectedMuonRefs;
  selectedMuonRefs.reserve(muonRefs.size());
  for ( unsigned int i=0; i<muonRefs.size(); ++i ) {
    // validity of the Ref
    if ( muonRefs[i].isNull() )  continue;
    // passed previous filter?
    if ( find(prevCandRefs.begin(),prevCandRefs.end(),muonRefs[i])==prevCandRefs.end() )
      continue;
    const reco::RecoChargedCandidate& muon = *muonRefs[i];
    // cuts on the momentum
    if ( muon.pt()<minMuonPt_ || muon.p()<minMuonP_ ||
	 fabs(muon.eta())>maxMuonEta_ )  continue;
    // cuts on the track associated to the muon
    if ( muon.track().isNull() )  continue;
    const reco::Track& track = *muon.track();
    if ( fabs(track.dxy(beamspot))>maxMuonDxy_ ||
	 fabs(track.dz(beamspot))>maxMuonDz_ )  continue;
    // keep muon
    selectedMuonRefs.push_back(muonRefs[i]);
  }
  //
  // access to tracks and selection according to configuration
  //
  std::vector<reco::RecoChargedCandidateRef> trackRefs;
  trackHandle->getObjects(trigger::TriggerTrack,trackRefs);
  //
  std::vector<reco::RecoChargedCandidateRef> selectedTrackRefs;
  selectedTrackRefs.reserve(trackRefs.size());
  for ( unsigned int i=0; i<trackRefs.size(); ++i ) {
    // validity of REF
    if ( trackRefs[i].isNull() )  continue;
    const reco::RecoChargedCandidate& trackCand = *trackRefs[i];
    // cuts on the momentum
    if ( trackCand.pt()<minTrackPt_ || trackCand.p()<minTrackP_ ||
	 fabs(trackCand.eta())>maxTrackEta_ )  continue;
    if ( trackCand.track().isNull() )  continue;
    // cuts on track quality
    const reco::Track& track = *trackCand.track();
    if ( fabs(track.dxy(beamspot))>maxTrackDxy_ ||
	 fabs(track.dz(beamspot))>maxTrackDz_ ||
	 track.numberOfValidHits()<minTrackHits_ ||
	 track.normalizedChi2()>maxTrackNormChi2_ )  continue;
    // keep track
    selectedTrackRefs.push_back(trackRefs[i]);
  }

  unsigned int nQ(0);
  unsigned int nSel(0);
  reco::Particle::LorentzVector p4Muon;
  reco::Particle::LorentzVector p4JPsi;
  for ( unsigned int im=0; im<muonRefs.size(); ++im ) {
    int qMuon = muonRefs[im]->charge();
    p4Muon = muonRefs[im]->p4();
    for ( unsigned int it=0; it<trackRefs.size(); ++it ) {
      if ( checkCharge_ && trackRefs[it]->charge()!=-qMuon )  continue;
      ++nQ;
      double mass = (p4Muon+trackRefs[it]->p4()).mass();
      for ( unsigned int j=0; j<minMasses_.size(); ++j ) {
	if ( mass>minMasses_[j] && mass<maxMasses_[j] ) {
	  ++nSel;
	  filterproduct->addObject(trigger::TriggerMuon,muonRefs[im]);
	  filterproduct->addObject(trigger::TriggerTrack,trackRefs[it]);
	  break;
	}
      }
    }
  }


  if ( edm::isDebugEnabled() ) {
    std::ostringstream stream;
    stream << "Total number of combinations = " << muonRefs.size()*trackRefs.size()
	   << " , after charge " << nQ << " , after mass " << nSel << std::endl;
    stream << "Found " << nSel << " jpsi candidates with # / mass / q / pt / eta" << std::endl;
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

  return nSel>0;
}


//define this as a plug-in
DEFINE_FWK_MODULE(HLTOniaMuonTrackMassFilter);
