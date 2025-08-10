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
        void analyze(const edm::Event&, const edm::EventSetup&) override;
        void endJob() override;

        // ----------member data ---------------------------
        edm::ESGetToken<CaloGeometry, CaloGeometryRecord> esToken;
        edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>> pfJetsToken_;
        edm::EDGetTokenT<std::vector<Run3ScoutingEBRecHit>> ebRecHitsToken_;
        TTree * timeTree;

        std::vector<double> *v_caloCellPhi = new std::vector<double>();
        std::vector<double> *v_caloCellEta = new std::vector<double>();
        std::vector<double> *v_caloCellEnergy = new std::vector<double>();
        std::vector<double> *v_caloCellEcalTime = new std::vector<double>();
        std::vector<double> *v_caloCellEcalTimeError = new std::vector<double>();
        std::vector<double> *v_pfJetE = new std::vector<double>();
        std::vector<double> *v_pfJetPt = new std::vector<double>();
        std::vector<double> *v_pfJetPhi = new std::vector<double>();
        std::vector<double> *v_pfJetEta = new std::vector<double>();
        std::vector<double> *v_pfJetWeightedTimeCell = new std::vector<double>();
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
void ScoutingTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    auto const& pG = iSetup.getData(esToken);
    auto const& pfJets = iEvent.get(pfJetsToken_);
    auto const& ebRecHits = iEvent.get(ebRecHitsToken_);
    v_caloCellPhi->clear();
    v_caloCellEta->clear();
    v_caloCellEnergy->clear();
    v_caloCellEcalTime->clear();
    v_caloCellEcalTimeError->clear();
    v_pfJetE->clear();
    v_pfJetPt->clear();
    v_pfJetPhi->clear();
    v_pfJetEta->clear();
    v_pfJetWeightedTimeCell->clear();
    v_pfJetChargedHadEnergy->clear();
    v_pfJetNeutralHadEnergy->clear();
    v_pfJetChargedEmEnergy->clear();
    v_pfJetNeutralEmEnergy->clear();
    v_caloCellPhi->clear();
    v_caloCellEta->clear();
    v_caloCellEnergy->clear();
    v_caloCellEcalTime->clear();
    v_caloCellEcalTimeError->clear();
    for (auto const& ebRecHit : ebRecHits) {
        if (ebRecHit.energy() < 0.5) continue;
        GlobalPoint pCell=pG.getPosition(ebRecHit.detId());
        TLorentzVector caloCellVecTemp;
        caloCellVecTemp.SetPtEtaPhiM(1,pCell.eta(),pCell.phi(),0);
        v_caloCellPhi->push_back(pCell.phi());
        v_caloCellEta->push_back(pCell.eta());
        v_caloCellEnergy->push_back(ebRecHit.energy());
        v_caloCellEcalTime->push_back(ebRecHit.time());
        // v_caloCellEcalTimeError->push_back(ebRecHit.timeError());
    }
    for (auto const& pfJet : pfJets) {
        TLorentzVector pfJetVecTemp;
        pfJetVecTemp.SetPtEtaPhiM(pfJet.pt(),pfJet.eta(),pfJet.phi(),0);
        double totalPtCell = 0;
        double  weightedTimeCellPt =0;
        int nCell = 0;

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
        else{weightedTimeCellPt = -200.0;totalPtCell=-1.0;}



        v_pfJetE->push_back(pfJet.chargedHadronEnergy()+pfJet.neutralHadronEnergy()+pfJet.electronEnergy()+pfJet.photonEnergy());
        v_pfJetPt->push_back(pfJet.pt());
        v_pfJetPhi->push_back(pfJet.phi());
        v_pfJetEta->push_back(pfJet.eta());
        v_pfJetWeightedTimeCell->push_back(weightedTimeCellPt);
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
