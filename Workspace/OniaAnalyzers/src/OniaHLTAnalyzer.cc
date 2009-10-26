// -*- C++ -*-
//
// Package:    OniaHLTAnalyzer
// Class:      OniaHLTAnalyzer
// 
/**\class OniaHLTAnalyzer OniaHLTAnalyzer.cc Workspace/OniaHLTAnalyzer/src/OniaHLTAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wolfgang ADAM
//         Created:  Fri Oct 16 11:26:41 CEST 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
//
// class decleration
//

class OniaHLTAnalyzer : public edm::EDAnalyzer {
public:
  struct HistogramSet {
    TH1* h_mult;
    TH1* h_p;
    TH1* h_pt;
    TH1* h_eta;
    TH1* h_phi;
    TH1* h_m;
    TH1* h_q;
  };

public:
  explicit OniaHLTAnalyzer(const edm::ParameterSet&);
  ~OniaHLTAnalyzer();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void analyzeCandidates (const reco::RecoChargedCandidateCollection& candidates,
			  HistogramSet& histos);
  void analyzeCandidates (const std::vector<reco::RecoChargedCandidateRef>& candidates,
			  HistogramSet& histos);
  
  // ----------member data ---------------------------
private:
  edm::InputTag muCandTag_;
  edm::InputTag tkCandTag_;
  edm::InputTag muFilterTag_;
  edm::InputTag tkFilterTag_;
  edm::InputTag jpsiTag_;

  HistogramSet muCandHistos_;
  HistogramSet tkCandHistos_;
  HistogramSet muFilterHistos_;
  HistogramSet tkFilterHistos_;
  HistogramSet jpsiHistos_;
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
OniaHLTAnalyzer::OniaHLTAnalyzer(const edm::ParameterSet& iConfig) :
  muCandTag_(iConfig.getParameter<edm::InputTag>("muonCandidates")),
  tkCandTag_(iConfig.getParameter<edm::InputTag>("trackCandidates")),
  muFilterTag_(iConfig.getParameter<edm::InputTag>("muonFilteredCandidates")),
  tkFilterTag_(iConfig.getParameter<edm::InputTag>("trackFilteredCandidates")),
  jpsiTag_(iConfig.getParameter<edm::InputTag>("jpsiCandidates"))
{
}


OniaHLTAnalyzer::~OniaHLTAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
OniaHLTAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

//   std::cout << "1" << std::endl;
   Handle<reco::RecoChargedCandidateCollection> muCandHandle;
   iEvent.getByLabel(muCandTag_,muCandHandle);
   if ( muCandHandle.isValid() )  analyzeCandidates(*muCandHandle,muCandHistos_);

//   std::cout << "2" << std::endl;
   Handle<reco::RecoChargedCandidateCollection> trackCandHandle;
   iEvent.getByLabel(tkCandTag_,trackCandHandle);
   if ( trackCandHandle.isValid() )  analyzeCandidates(*trackCandHandle,tkCandHistos_);

//   std::cout << "3" << std::endl;
   Handle<trigger::TriggerFilterObjectWithRefs> muFilterHandle;
   iEvent.getByLabel(muFilterTag_,muFilterHandle);
   if ( muFilterHandle.isValid() ) {
//   std::cout << "3a" << std::endl;
     std::vector<reco::RecoChargedCandidateRef> filteredMuons;
     muFilterHandle->getObjects(trigger::TriggerMuon,filteredMuons);
//   std::cout << "3b" << std::endl;
     analyzeCandidates(filteredMuons,muFilterHistos_);
   }

//   std::cout << "4" << std::endl;
   Handle<trigger::TriggerFilterObjectWithRefs> tkFilterHandle;
   iEvent.getByLabel(tkFilterTag_,tkFilterHandle);
   if ( tkFilterHandle.isValid() ) {
//   std::cout << "4a" << std::endl;
     std::vector<reco::RecoChargedCandidateRef> filteredTracks;
     tkFilterHandle->getObjects(trigger::TriggerTrack,filteredTracks);
//   std::cout << "4b" << std::endl;
     analyzeCandidates(filteredTracks,tkFilterHistos_);
   }

//   std::cout << "5" << std::endl;
   Handle<trigger::TriggerFilterObjectWithRefs> jpsiHandle;
   iEvent.getByLabel(jpsiTag_,jpsiHandle);

   if ( jpsiHandle.isValid() ) {
//   std::cout << "5a" << std::endl;

     std::vector<reco::RecoChargedCandidateRef> jpsiMuons;
     jpsiHandle->getObjects(trigger::TriggerMuon,jpsiMuons);
     std::vector<reco::RecoChargedCandidateRef> jpsiTracks;
     jpsiHandle->getObjects(trigger::TriggerTrack,jpsiTracks);
//   std::cout << "5b" << std::endl;

     if ( jpsiMuons.size()==jpsiTracks.size() ) {
       jpsiHistos_.h_mult->Fill(jpsiMuons.size());
       for ( unsigned int i=0; i<jpsiMuons.size(); ++i ) {
	 reco::Particle::LorentzVector jpsiP4 = 
	   jpsiMuons[i]->p4() + jpsiTracks[i]->p4();
	 jpsiHistos_.h_p->Fill(jpsiP4.P());
	 jpsiHistos_.h_pt->Fill(jpsiP4.Pt());
	 jpsiHistos_.h_eta->Fill(jpsiP4.Eta());
	 jpsiHistos_.h_phi->Fill(jpsiP4.Phi());
	 jpsiHistos_.h_m->Fill(jpsiP4.M());
	 jpsiHistos_.h_q->Fill(jpsiMuons[i]->charge()+jpsiTracks[i]->charge());
       }
     }
     else {
       edm::LogError("OniaHLTAnalyzer") << "inconsistent size of jpsi muon and track candidates: "
					<< jpsiMuons.size() << " / " << jpsiTracks.size();
     }
   }
   
}

void
OniaHLTAnalyzer::analyzeCandidates (const reco::RecoChargedCandidateCollection& candidates,
				    HistogramSet& histos) 
{
  histos.h_mult->Fill(candidates.size());
  for ( unsigned int i=0; i<candidates.size(); ++i ) {
//    std::cout << "analyzeCandidates(1): " << i << std::endl;
    const reco::RecoChargedCandidate& cand = candidates[i];
    histos.h_p->Fill(cand.p());
    histos.h_pt->Fill(cand.pt());
    histos.h_eta->Fill(cand.eta());
    histos.h_phi->Fill(cand.phi());
    histos.h_m->Fill(cand.mass());
    histos.h_q->Fill(cand.charge());
  }
}

void
OniaHLTAnalyzer::analyzeCandidates (const std::vector<reco::RecoChargedCandidateRef>& candidates,
				    HistogramSet& histos) 
{
  histos.h_mult->Fill(candidates.size());
  for ( unsigned int i=0; i<candidates.size(); ++i ) {
//    std::cout << "analyzeCandidates(2): " << i << " " << candidates[i].isNonnull()
//	      << " " << candidates[i].isAvailable() 
//	      << " " << candidates[i].get() << std::endl;
    const reco::RecoChargedCandidate& cand = *candidates[i];
//    std::cout << "...1" << std::endl;
    histos.h_p->Fill(cand.p());
    histos.h_pt->Fill(cand.pt());
    histos.h_eta->Fill(cand.eta());
    histos.h_phi->Fill(cand.phi());
    histos.h_m->Fill(cand.mass());
    histos.h_q->Fill(cand.charge());
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
OniaHLTAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

  muCandHistos_.h_mult = fs->make<TH1F>("muonCandMult","muon candidate multiplicity",21,-0.5,20.5);
  muCandHistos_.h_p = fs->make<TH1F>("muonCandP","muon candidate momentum",150,0,75);
  muCandHistos_.h_pt = fs->make<TH1F>("muonCandPt","muon candidate transverse momentum",125,0,25);
  muCandHistos_.h_eta = fs->make<TH1F>("muonCandEta","muon candidate |eta|",100,0,2.5);
  muCandHistos_.h_phi = fs->make<TH1F>("muonCandPhi","muon candidate phi",100,-3.14159265,3.14159265);
  muCandHistos_.h_m = fs->make<TH1F>("muonCandMass","muon candidate mass",100,0.,1.);
  muCandHistos_.h_q = fs->make<TH1F>("muonCandCharge","muon candidate charge",5,-2.5,2.5);

  tkCandHistos_.h_mult = fs->make<TH1F>("trackCandMult","track candidate multiplicity",101,-0.5,100.5);
  tkCandHistos_.h_p = fs->make<TH1F>("trackCandP","track candidate momentum",150,0,75);
  tkCandHistos_.h_pt = fs->make<TH1F>("trackCandPt","track candidate transverse momentum",125,0,25);
  tkCandHistos_.h_eta = fs->make<TH1F>("trackCandEta","track candidate |eta|",100,0,2.5);
  tkCandHistos_.h_phi = fs->make<TH1F>("trackCandPhi","track candidate phi",100,-3.14159265,3.14159265);
  tkCandHistos_.h_m = fs->make<TH1F>("trackCandMass","track candidate mass",100,0.,1.);
  tkCandHistos_.h_q = fs->make<TH1F>("trackCandCharge","track candidate charge",5,-2.5,2.5);

  muFilterHistos_.h_mult = fs->make<TH1F>("muonFiltMult","filtered muon multiplicity",21,-0.5,20.5);
  muFilterHistos_.h_p = fs->make<TH1F>("muonFiltP","filtered muon momentum",150,0,75);
  muFilterHistos_.h_pt = fs->make<TH1F>("muonFiltPt","filtered muon transverse momentum",125,0,25);
  muFilterHistos_.h_eta = fs->make<TH1F>("muonFiltEta","filtered muon |eta|",100,0,2.5);
  muFilterHistos_.h_phi = fs->make<TH1F>("muonFiltPhi","filtered muon phi",100,-3.14159265,3.14159265);
  muFilterHistos_.h_m = fs->make<TH1F>("muonFiltMass","filtered muon mass",100,0.,1.);
  muFilterHistos_.h_q = fs->make<TH1F>("muonFiltCharge","filtered muon charge",5,-2.5,2.5);

  tkFilterHistos_.h_mult = fs->make<TH1F>("trackFiltMult","filtered track multiplicity",101,-0.5,100.5);
  tkFilterHistos_.h_p = fs->make<TH1F>("trackFiltP","filtered track momentum",150,0,75);
  tkFilterHistos_.h_pt = fs->make<TH1F>("trackFiltPt","filtered track transverse momentum",125,0,25);
  tkFilterHistos_.h_eta = fs->make<TH1F>("trackFiltEta","filtered track |eta|",100,0,2.5);
  tkFilterHistos_.h_phi = fs->make<TH1F>("trackFiltPhi","filtered track phi",100,-3.14159265,3.14159265);
  tkFilterHistos_.h_m = fs->make<TH1F>("trackFiltMass","filtered track mass",100,0.,1.);
  tkFilterHistos_.h_q = fs->make<TH1F>("trackFiltCharge","filtered track charge",5,-2.5,2.5);

  jpsiHistos_.h_mult = fs->make<TH1F>("jpsiMult","jpsi multiplicity",21,-0.5,20.5);
  jpsiHistos_.h_p = fs->make<TH1F>("jpsiP","jpsi momentum",150,0,75);
  jpsiHistos_.h_pt = fs->make<TH1F>("jpsiPt","jpsi transverse momentum",125,0,25);
  jpsiHistos_.h_eta = fs->make<TH1F>("jpsiEta","jpsi |eta|",100,0,2.5);
  jpsiHistos_.h_phi = fs->make<TH1F>("jpsiPhi","jpsi phi",100,-3.14159265,3.14159265);
  jpsiHistos_.h_m = fs->make<TH1F>("jpsiMass","jpsi mass",500,0.,10.);
  jpsiHistos_.h_q = fs->make<TH1F>("jpsiCharge","jpsi charge",5,-2.5,2.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
OniaHLTAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(OniaHLTAnalyzer);
