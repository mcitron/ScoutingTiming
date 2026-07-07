// -*- C++ -*-
//
// Package:    ScoutingTiming/ScoutingTimingAnalyzer
// Class:      ScoutingTimingAnalyzer
//
/**\class ScoutingTimingAnalyzer ScoutingTimingAnalyzer.cc ScoutingTiming/ScoutingTimingAnalyzer/plugins/ScoutingTimingAnalyzer.cc

 Description: Compute per-PF-jet ECAL timing from Run3 scouting EB rec hits.

 Implementation:
     For each scouting PF jet, take an ET-weighted mean of the rec hit time
     over EB rec hits with E > 0.5 GeV within dR < 0.4 of the jet axis.
     Cell-level info (eta, phi, E, time, flag bitfield) is also stored so
     flag-based selections can be applied offline.
*/
//
// Original Author:  Matthew Daniel Citron
//         Created:  Sat, 09 Aug 2025 04:38:14 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingEBRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

using reco::TrackCollection;

class ScoutingTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingTimingAnalyzer(const edm::ParameterSet&);
  ~ScoutingTimingAnalyzer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override {}
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  // ----------member data ---------------------------
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> esToken;
  edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>>    pfJetsToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingEBRecHit>> ebRecHitsToken_;
  edm::EDGetTokenT<edm::TriggerResults>               triggerResultsToken;
  edm::EDGetTokenT<edm::TriggerResults>               triggerResultsRerunToken;

  TTree* timeTree;

  // EB rec hit (calo cell) info. Per-flag booleans are deliberately not
  // stored; decode caloCell_flags offline via (flags >> bit) & 0x1.
  std::vector<double>   v_caloCellPhi;
  std::vector<double>   v_caloCellEta;
  std::vector<double>   v_caloCellEnergy;
  std::vector<double>   v_caloCellEcalTime;
  std::vector<uint32_t> v_caloCellFlags;

  // PF jet info
  std::vector<double>   v_pfJetE;
  std::vector<double>   v_pfJetPt;
  std::vector<double>   v_pfJetPhi;
  std::vector<double>   v_pfJetEta;
  std::vector<double>   v_pfJetWeightedTimeCell;
  std::vector<double>   v_pfJetTotalPtCell;
  std::vector<uint32_t> v_pfJetNCell;
  std::vector<double>   v_pfJetChargedHadEnergy;
  std::vector<double>   v_pfJetNeutralHadEnergy;
  std::vector<double>   v_pfJetChargedEmEnergy;
  std::vector<double>   v_pfJetNeutralEmEnergy;

  // Trigger decisions
  bool b_delayedJetPathPass;
  bool b_scoutingJetPathPass;
};

//
// constructors and destructor
//
ScoutingTimingAnalyzer::ScoutingTimingAnalyzer(const edm::ParameterSet& iPSet)
    : esToken(esConsumes()),
      pfJetsToken_(consumes(iPSet.getParameter<edm::InputTag>("pfJetsTag"))),
      ebRecHitsToken_(consumes(iPSet.getParameter<edm::InputTag>("ebRecHitsTag"))),
      triggerResultsToken(consumes(iPSet.getParameter<edm::InputTag>("triggerResultsTag"))),
      triggerResultsRerunToken(consumes(iPSet.getParameter<edm::InputTag>("triggerResultsRerunTag"))) {
  edm::Service<TFileService> fs;
  timeTree = fs->make<TTree>("timeTree", "timeTree");

  // calo cell branches
  timeTree->Branch("caloCell_eta",      &v_caloCellEta);
  timeTree->Branch("caloCell_phi",      &v_caloCellPhi);
  timeTree->Branch("caloCell_e",        &v_caloCellEnergy);
  timeTree->Branch("caloCell_ecalTime", &v_caloCellEcalTime);
  timeTree->Branch("caloCell_flags",    &v_caloCellFlags);

  // PF jet branches
  timeTree->Branch("pfJet_e",                &v_pfJetE);
  timeTree->Branch("pfJet_pt",               &v_pfJetPt);
  timeTree->Branch("pfJet_phi",              &v_pfJetPhi);
  timeTree->Branch("pfJet_eta",              &v_pfJetEta);
  timeTree->Branch("pfJet_weightedTime",     &v_pfJetWeightedTimeCell);
  timeTree->Branch("pfJet_totalPtCell",      &v_pfJetTotalPtCell);
  timeTree->Branch("pfJet_nCell",            &v_pfJetNCell);
  timeTree->Branch("pfJet_chargedHadEnergy", &v_pfJetChargedHadEnergy);
  timeTree->Branch("pfJet_neutralHadEnergy", &v_pfJetNeutralHadEnergy);
  timeTree->Branch("pfJet_chargedEmEnergy",  &v_pfJetChargedEmEnergy);
  timeTree->Branch("pfJet_neutralEmEnergy",  &v_pfJetNeutralEmEnergy);

  // trigger decision branches
  timeTree->Branch("delayedJetPathPass",  &b_delayedJetPathPass);
  timeTree->Branch("scoutingJetPathPass", &b_scoutingJetPathPass);
}

//
// member functions
//

// ------------ method called for each event  ------------
void ScoutingTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  auto const& pG                  = iSetup.getData(esToken);
  auto const& pfJets              = iEvent.get(pfJetsToken_);
  auto const& ebRecHits           = iEvent.get(ebRecHitsToken_);
  auto const& triggerResults      = iEvent.get(triggerResultsToken);
  auto const& triggerResultsRerun = iEvent.get(triggerResultsRerunToken);

  // clear per-event containers
  v_caloCellPhi.clear();
  v_caloCellEta.clear();
  v_caloCellEnergy.clear();
  v_caloCellEcalTime.clear();
  v_caloCellFlags.clear();

  v_pfJetE.clear();
  v_pfJetPt.clear();
  v_pfJetPhi.clear();
  v_pfJetEta.clear();
  v_pfJetWeightedTimeCell.clear();
  v_pfJetTotalPtCell.clear();
  v_pfJetNCell.clear();
  v_pfJetChargedHadEnergy.clear();
  v_pfJetNeutralHadEnergy.clear();
  v_pfJetChargedEmEnergy.clear();
  v_pfJetNeutralEmEnergy.clear();

  // -------- trigger decisions --------
  b_delayedJetPathPass  = false;
  b_scoutingJetPathPass = false;

  {
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(triggerResults);
    for (unsigned int i = 0; i < triggerResults.size(); ++i) {
      const auto& trigname = triggerNames.triggerName(i);
      if (trigname.find("HLT_HT430_DelayedJet40_SingleDelay2nsInclusive_v") != std::string::npos)
        b_delayedJetPathPass = triggerResults.accept(i);
    }
  }

  {
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(triggerResultsRerun);
    for (unsigned int i = 0; i < triggerResultsRerun.size(); ++i) {
      const auto& trigname = triggerNames.triggerName(i);
      if (trigname.find("HLT_HT430_DelayedJet40_SingleDelay2nsInclusive_v") != std::string::npos)
        b_scoutingJetPathPass = triggerResultsRerun.accept(i);
    }
  }

  // -------- EB rec hits (calo cells) --------
  for (auto const& ebRecHit : ebRecHits) {
    if (ebRecHit.energy() < 0.5) continue;
    GlobalPoint pCell = pG.getPosition(ebRecHit.detId());

    v_caloCellPhi.push_back(pCell.phi());
    v_caloCellEta.push_back(pCell.eta());
    v_caloCellEnergy.push_back(ebRecHit.energy());
    v_caloCellEcalTime.push_back(ebRecHit.time());
    v_caloCellFlags.push_back(ebRecHit.flags());
  }

  // -------- PF jets: ET-weighted mean cell time within dR < 0.4 --------
  for (auto const& pfJet : pfJets) {
    TLorentzVector pfJetVecTemp;
    pfJetVecTemp.SetPtEtaPhiM(pfJet.pt(), pfJet.eta(), pfJet.phi(), 0);

    double totalPtCell        = 0;
    double weightedTimeCellPt = 0;
    uint32_t nCell            = 0;

    for (auto const& ebRecHit : ebRecHits) {
      if (ebRecHit.energy() < 0.5) continue;
      GlobalPoint pCell = pG.getPosition(ebRecHit.detId());

      TLorentzVector caloCellVecTemp;
      caloCellVecTemp.SetPtEtaPhiM(1, pCell.eta(), pCell.phi(), 0);
      if (caloCellVecTemp.DeltaR(pfJetVecTemp) > 0.4) continue;

      // Deliberately no flag-based cell-level rejection here: definition
      // matches earlier CMS delayed-jet analyses, and tighter selections
      // can be applied offline via caloCell_flags.

      const double cellPt = ebRecHit.energy() * TMath::Sin(pCell.theta());
      weightedTimeCellPt += ebRecHit.time() * cellPt;
      totalPtCell        += cellPt;
      ++nCell;
    }

    if (totalPtCell > 0) weightedTimeCellPt /= totalPtCell;
    else                 weightedTimeCellPt = -200.0;

    v_pfJetE.push_back(pfJet.chargedHadronEnergy() + pfJet.neutralHadronEnergy()
                       + pfJet.electronEnergy() + pfJet.photonEnergy());
    v_pfJetPt.push_back(pfJet.pt());
    v_pfJetPhi.push_back(pfJet.phi());
    v_pfJetEta.push_back(pfJet.eta());
    v_pfJetWeightedTimeCell.push_back(weightedTimeCellPt);
    v_pfJetTotalPtCell.push_back(totalPtCell);
    v_pfJetNCell.push_back(nCell);
    v_pfJetChargedHadEnergy.push_back(pfJet.chargedHadronEnergy());
    v_pfJetNeutralHadEnergy.push_back(pfJet.neutralHadronEnergy());
    v_pfJetChargedEmEnergy.push_back(pfJet.electronEnergy());
    v_pfJetNeutralEmEnergy.push_back(pfJet.photonEnergy());
  }

  timeTree->Fill();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("pfJetsTag");
  desc.add<edm::InputTag>("ebRecHitsTag");
  desc.add<edm::InputTag>("triggerResultsTag");
  desc.add<edm::InputTag>("triggerResultsRerunTag");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingTimingAnalyzer);
