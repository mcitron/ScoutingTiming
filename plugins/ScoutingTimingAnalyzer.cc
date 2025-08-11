// -*- C++ -*-
//
// Package:    ScoutingTiming/ScoutingTimingAnalyzer
// Class:      ScoutingTimingAnalyzer
//
/**\class ScoutingTimingAnalyzer ScoutingTimingAnalyzer.cc ScoutingTiming/ScoutingTimingAnalyzer/plugins/ScoutingTimingAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
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

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "TLorentzVector.h"
#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class ScoutingTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit ScoutingTimingAnalyzer(const edm::ParameterSet&);
        ~ScoutingTimingAnalyzer() override;

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        void beginJob() override;
        bool checkFlag(uint32_t flag,uint32_t flagBits);
        void analyze(const edm::Event&, const edm::EventSetup&) override;
        void endJob() override;

        // ----------member data ---------------------------
        edm::ESGetToken<CaloGeometry, CaloGeometryRecord> esToken;
        edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>> pfJetsToken_;
        edm::EDGetTokenT<std::vector<Run3ScoutingEBRecHit>> ebRecHitsToken_;
        TTree * timeTree;

        std::vector<double> *v_caloCellPhi = new std::vector<double>();
        std::vector<double> *v_caloCellEta = new std::vector<double>();
	std::vector<uint32_t> *v_caloCellFlags= new std::vector<uint32_t>();
	std::vector<bool> *v_caloCellkGood= new std::vector<bool>();
	std::vector<bool> *v_caloCellkPoorReco= new std::vector<bool>();
	std::vector<bool> *v_caloCellkOutOfTime= new std::vector<bool>();
	std::vector<bool> *v_caloCellkFaultyHardware= new std::vector<bool>();
	std::vector<bool> *v_caloCellkNoisy= new std::vector<bool>();
	std::vector<bool> *v_caloCellkPoorCalib= new std::vector<bool>();
	std::vector<bool> *v_caloCellkSaturated= new std::vector<bool>();
	std::vector<bool> *v_caloCellkLeadingEdgeRecovered= new std::vector<bool>();
	std::vector<bool> *v_caloCellkNeighboursRecovered= new std::vector<bool>();
	std::vector<bool> *v_caloCellkTowerRecovered= new std::vector<bool>();
	std::vector<bool> *v_caloCellkDead= new std::vector<bool>();
	std::vector<bool> *v_caloCellkKilled= new std::vector<bool>();
	std::vector<bool> *v_caloCellkTPSaturated= new std::vector<bool>();
	std::vector<bool> *v_caloCellkL1SpikeFlag= new std::vector<bool>();
	std::vector<bool> *v_caloCellkWeird= new std::vector<bool>();
	std::vector<bool> *v_caloCellkDiWeird= new std::vector<bool>();
	std::vector<bool> *v_caloCellkHasSwitchToGain6= new std::vector<bool>();
	std::vector<bool> *v_caloCellkHasSwitchToGain1= new std::vector<bool>();
	std::vector<double> *v_caloCellEnergy = new std::vector<double>();
	std::vector<double> *v_caloCellEcalTime = new std::vector<double>();
	std::vector<double> *v_caloCellEcalTimeError = new std::vector<double>();
	std::vector<double> *v_pfJetE = new std::vector<double>();
	std::vector<double> *v_pfJetPt = new std::vector<double>();
	std::vector<double> *v_pfJetPhi = new std::vector<double>();
	std::vector<double> *v_pfJetEta = new std::vector<double>();
	std::vector<double> *v_pfJetWeightedTimeCell = new std::vector<double>();
	std::vector<double> *v_pfJetTotalPtCell = new std::vector<double>();
	std::vector<uint32_t> *v_pfJetNCell = new std::vector<uint32_t>();
	std::vector<double> *v_pfJetChargedHadEnergy = new std::vector<double>();
	std::vector<double> *v_pfJetNeutralHadEnergy = new std::vector<double>();
	std::vector<double> *v_pfJetChargedEmEnergy = new std::vector<double>();
	std::vector<double> *v_pfJetNeutralEmEnergy = new std::vector<double>();
	// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	//   edm::ESGetToken<SetupData, SetupRecord> setupToken_;
	// #endif
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
ScoutingTimingAnalyzer::ScoutingTimingAnalyzer(const edm::ParameterSet& iPSet)
    : esToken(esConsumes()),
    pfJetsToken_(consumes(iPSet.getParameter<edm::InputTag>("pfJetsTag"))),
    ebRecHitsToken_(consumes(iPSet.getParameter<edm::InputTag>("ebRecHitsTag"))){
	// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	//   setupDataToken_ = esConsumes<SetupData, SetupRecord>();
	// #endif
	//now do what ever initialization is needed
	edm::Service<TFileService> fs;
	timeTree = fs->make<TTree>("timeTree","timeTree");
	timeTree->Branch("caloCell_eta",&v_caloCellEta);
	timeTree->Branch("caloCell_kGood",v_caloCellkGood );
	timeTree->Branch("caloCell_kPoorReco",v_caloCellkPoorReco);
	timeTree->Branch("caloCell_kOutOfTime",v_caloCellkOutOfTime);
	timeTree->Branch("caloCell_kFaultyHardware",v_caloCellkFaultyHardware);
	timeTree->Branch("caloCell_kNoisy",v_caloCellkNoisy);
	timeTree->Branch("caloCell_kPoorCalib",v_caloCellkPoorCalib);
	timeTree->Branch("caloCell_kSaturated",v_caloCellkSaturated);
	timeTree->Branch("caloCell_kLeadingEdgeRecovered",v_caloCellkLeadingEdgeRecovered);
	timeTree->Branch("caloCell_kNeighboursRecovered",v_caloCellkNeighboursRecovered);
	timeTree->Branch("caloCell_kTowerRecovered",v_caloCellkTowerRecovered);
	timeTree->Branch("caloCell_kDead",v_caloCellkDead);
	timeTree->Branch("caloCell_kKilled",v_caloCellkKilled);
	timeTree->Branch("caloCell_kTPSaturated",v_caloCellkTPSaturated);
	timeTree->Branch("caloCell_kL1SpikeFlag",v_caloCellkL1SpikeFlag);
	timeTree->Branch("caloCell_kWeird",v_caloCellkWeird);
	timeTree->Branch("caloCell_kDiWeird",v_caloCellkDiWeird);
	timeTree->Branch("caloCell_kHasSwitchToGain6",v_caloCellkHasSwitchToGain6);
	timeTree->Branch("caloCell_kHasSwitchToGain1",v_caloCellkHasSwitchToGain1);

	timeTree->Branch("caloCell_flags",&v_caloCellFlags);
	// timeTree->Branch("caloCell_detId",&v_caloCelldetId);
	// timeTree->Branch("caloCell_good",&v_caloCellGood);
	// timeTree->Branch("caloCell_goodNoOOT",&v_caloCellGoodNoOOT);
	// timeTree->Branch("caloCell_OOT",&v_caloCellOOT);
	timeTree->Branch("caloCell_e",&v_caloCellEnergy);
	timeTree->Branch("caloCell_phi",&v_caloCellPhi);
	// timeTree->Branch("caloCell_iphi",&v_caloCelliPhi);
	// timeTree->Branch("caloCell_ieta",&v_caloCelliEta);
	timeTree->Branch("caloCell_ecalTime",&v_caloCellEcalTime);
	timeTree->Branch("caloCell_ecalTimeError",&v_caloCellEcalTimeError);

	timeTree->Branch("pfJet_e",&v_pfJetE);
	timeTree->Branch("pfJet_pt",&v_pfJetPt);
	timeTree->Branch("pfJet_phi",&v_pfJetPhi);
	timeTree->Branch("pfJet_eta",&v_pfJetEta);

	timeTree->Branch("pfJet_weightedTime", &v_pfJetWeightedTimeCell);
	timeTree->Branch("pfJet_totalPtCell", &v_pfJetTotalPtCell);
	timeTree->Branch("pfJet_nCell", &v_pfJetNCell);

	timeTree->Branch("pfJet_chargedHadEnergy",&v_pfJetChargedHadEnergy);
	timeTree->Branch("pfJet_neutralHadEnergy",&v_pfJetNeutralHadEnergy);
	timeTree->Branch("pfJet_chargedEmEnergy",&v_pfJetChargedEmEnergy);
	timeTree->Branch("pfJet_neutralEmEnergy", &v_pfJetNeutralEmEnergy);
    }

ScoutingTimingAnalyzer::~ScoutingTimingAnalyzer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    //
    // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
bool ScoutingTimingAnalyzer::checkFlag(uint32_t flag,uint32_t flagBits){ return flagBits & (0x1 << flag); }
void ScoutingTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    auto const& pG = iSetup.getData(esToken);
    auto const& pfJets = iEvent.get(pfJetsToken_);
    auto const& ebRecHits = iEvent.get(ebRecHitsToken_);
    v_caloCellPhi->clear();
    v_caloCellEta->clear();
    v_caloCellFlags->clear();
    v_caloCellEnergy->clear();
    v_caloCellEcalTime->clear();
    v_caloCellEcalTimeError->clear();
    v_caloCellkGood->clear();
    v_caloCellkPoorReco->clear();
    v_caloCellkOutOfTime->clear();
    v_caloCellkFaultyHardware->clear();
    v_caloCellkNoisy->clear();
    v_caloCellkPoorCalib->clear();
    v_caloCellkSaturated->clear();
    v_caloCellkLeadingEdgeRecovered->clear();
    v_caloCellkNeighboursRecovered->clear();
    v_caloCellkTowerRecovered->clear();
    v_caloCellkDead->clear();
    v_caloCellkKilled->clear();
    v_caloCellkTPSaturated->clear();
    v_caloCellkL1SpikeFlag->clear();
    v_caloCellkWeird->clear();
    v_caloCellkDiWeird->clear();
    v_caloCellkHasSwitchToGain6->clear();
    v_caloCellkHasSwitchToGain1->clear();
    v_pfJetE->clear();
    v_pfJetPt->clear();
    v_pfJetPhi->clear();
    v_pfJetEta->clear();
    v_pfJetWeightedTimeCell->clear();
    v_pfJetNCell->clear();
    v_pfJetTotalPtCell->clear();
    v_pfJetChargedHadEnergy->clear();
    v_pfJetNeutralHadEnergy->clear();
    v_pfJetChargedEmEnergy->clear();
    v_pfJetNeutralEmEnergy->clear();
    for (auto const& ebRecHit : ebRecHits) {
	if (ebRecHit.energy() < 0.5) continue;
	GlobalPoint pCell=pG.getPosition(ebRecHit.detId());
	TLorentzVector caloCellVecTemp;
	caloCellVecTemp.SetPtEtaPhiM(1,pCell.eta(),pCell.phi(),0);
	v_caloCellPhi->push_back(pCell.phi());
	v_caloCellEta->push_back(pCell.eta());
	v_caloCellEnergy->push_back(ebRecHit.energy());
	v_caloCellEcalTime->push_back(ebRecHit.time());
	uint32_t flags = ebRecHit.flags();
	v_caloCellFlags->push_back(flags);
	v_caloCellkGood->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kGood,flags));
	v_caloCellkPoorReco->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kPoorReco,flags));
	v_caloCellkOutOfTime->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kOutOfTime,flags));
	v_caloCellkFaultyHardware->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kFaultyHardware,flags));
	v_caloCellkNoisy->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kNoisy,flags));
	v_caloCellkPoorCalib->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kPoorCalib,flags));
	v_caloCellkSaturated->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kSaturated,flags));
	v_caloCellkLeadingEdgeRecovered->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kLeadingEdgeRecovered,flags));
	v_caloCellkNeighboursRecovered->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kNeighboursRecovered,flags));
	v_caloCellkTowerRecovered->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kTowerRecovered,flags));
	v_caloCellkDead->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kDead,flags));
	v_caloCellkKilled->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kKilled,flags));
	v_caloCellkTPSaturated->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kTPSaturated,flags));
	v_caloCellkL1SpikeFlag->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kL1SpikeFlag,flags));
	v_caloCellkWeird->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kWeird,flags));
	v_caloCellkDiWeird->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kDiWeird,flags));
	v_caloCellkHasSwitchToGain6->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kHasSwitchToGain6,flags));
	v_caloCellkHasSwitchToGain1->push_back(ScoutingTimingAnalyzer::checkFlag(EcalRecHit::kHasSwitchToGain1,flags));
	// v_caloCellEcalTimeError->push_back(ebRecHit.timeError());
    }
    for (auto const& pfJet : pfJets) {
	TLorentzVector pfJetVecTemp;
	pfJetVecTemp.SetPtEtaPhiM(pfJet.pt(),pfJet.eta(),pfJet.phi(),0);
	double totalPtCell = 0;
	double  weightedTimeCellPt =0;
	uint32_t nCell = 0;

	for (auto const& ebRecHit : ebRecHits) {
	    if (ebRecHit.energy() < 0.5) continue;
	    GlobalPoint pCell=pG.getPosition(ebRecHit.detId());
	    TLorentzVector caloCellVecTemp;
	    caloCellVecTemp.SetPtEtaPhiM(1,pCell.eta(),pCell.phi(),0);
	    if (caloCellVecTemp.DeltaR(pfJetVecTemp) > 0.4) continue;
	    // if (ebRecHit.checkFlag(EcalRecHit::kSaturated) || ebRecHit.checkFlag(EcalRecHit::kLeadingEdgeRecovered) || ebRecHit.checkFlag(EcalRecHit::kPoorReco)) {
	    //      continue;
	    //  }
	    //  if (ebRecHit.checkFlag(EcalRecHit::kWeird)){
	    //      continue;
	    //  }
	    //  if (ebRecHit.checkFlag(EcalRecHit::kDiWeird)){
	    //      continue;
	    //  }
	    // if (ebRecHit.timeError() < 0 || ebRecHit.timeError() > 100) {
	    //     continue;
	    // }
	    weightedTimeCellPt += ebRecHit.time()*ebRecHit.energy()*TMath::Sin(pCell.theta());	
	    totalPtCell += ebRecHit.energy()*TMath::Sin(pCell.theta());
	    nCell ++;
	}
	if (totalPtCell>0) weightedTimeCellPt/=totalPtCell;
	else{weightedTimeCellPt = -200.0;}



	v_pfJetE->push_back(pfJet.chargedHadronEnergy()+pfJet.neutralHadronEnergy()+pfJet.electronEnergy()+pfJet.photonEnergy());
	v_pfJetPt->push_back(pfJet.pt());
	v_pfJetPhi->push_back(pfJet.phi());
	v_pfJetEta->push_back(pfJet.eta());
	v_pfJetWeightedTimeCell->push_back(weightedTimeCellPt);
	v_pfJetTotalPtCell->push_back(totalPtCell);
	v_pfJetNCell->push_back(nCell);
	v_pfJetChargedHadEnergy->push_back(pfJet.chargedHadronEnergy());
	v_pfJetNeutralHadEnergy->push_back(pfJet.neutralHadronEnergy());
	v_pfJetChargedEmEnergy->push_back(pfJet.electronEnergy());
	v_pfJetNeutralEmEnergy->push_back(pfJet.photonEnergy());
	// for (auto vertex = vertexCollection->begin(); vertex != vertexCollection->end(); vertex++){
	//   double ptPVTracks = 0.;	
	//   int nTracksPVTemp = 0;
	//   for(auto pvTrack=vertex->tracks_begin(); pvTrack!=vertex->tracks_end(); pvTrack++){
	//       FreeTrajectoryState ftspv = trajectoryStateTransform::initialFreeState (**pvTrack, magneticFieldPV_); 
	//       TrajectoryStateOnSurface outerpv = stateOnTrackerPV(ftspv); 
	//       if(!outerpv.isValid()) continue; 
	//       GlobalPoint outerpvPos = outerpv.globalPosition();
	//
	//       TLorentzVector pvTrackVecTemp;
	//       pvTrackVecTemp.SetPtEtaPhiM((*pvTrack)->pt(),outerpvPos.eta(),outerpvPos.phi(),0);
	//       //If pv track associated with jet add pt to ptPVTracks
	//       double dR = pvTrackVecTemp.DeltaR(pfJetVecTemp); 
	//       if ((*pvTrack)->pt() > 0.5){
	// 	  if (dR < dRClosestPVTrack0p5) dRClosestPVTrack0p5 = dR;
	//       }
	//   }
	// }
    }
    timeTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingTimingAnalyzer::beginJob() {
    // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingTimingAnalyzer::endJob() {
    // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingTimingAnalyzer);
