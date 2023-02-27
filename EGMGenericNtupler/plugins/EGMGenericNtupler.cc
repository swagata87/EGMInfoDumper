// -*- C++ -*-
// EGMGenericNtupler
// Author:  Swagata Mukherjee

// system include files
#include <memory>

// user include files
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "TTree.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

class EGMGenericNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit EGMGenericNtupler(const edm::ParameterSet&);
  ~EGMGenericNtupler();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  //  edm::Service<TFileService> fs;
  //TTree   *tree = fs->make<TTree>("EventTree", "EventData");

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------

  //  edm::EDGetTokenT<edm::View<reco::GsfElectron> > eleToken_;
  //edm::EDGetTokenT<edm::View<reco::Photon> > phoToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > eleToken_;
  edm::EDGetTokenT<edm::View<pat::Photon> > phoToken_;
};

EGMGenericNtupler::EGMGenericNtupler(const edm::ParameterSet& iConfig)
  :
  eleToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  phoToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("Photons")))
{}

EGMGenericNtupler::~EGMGenericNtupler()
{
}

void
EGMGenericNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  std::cout << "\n \n --New Event-- \n" ;

  using namespace edm;
  
  ////photon
  for(const auto& pho : iEvent.get(phoToken_) ) {
    std::cout << "\n new photon " << std::endl;
    
    //Run2 ID
    std::cout << "PhoID Cutbased Fall17V2 loose pass? " << pho.photonID("cutBasedPhotonID-Fall17-94X-V2-loose") << std::endl;
    //Run3 ID
    std::cout << "PhoID Cutbased Winter22 loose pass? " << pho.photonID("cutBasedPhotonID-RunIIIWinter22-122X-V1-loose") << std::endl;       
  }

  /////electron
  for(const auto& ele : iEvent.get(eleToken_) ) {
    std::cout << "\n new electron .... " << std::endl;

    //Run2 ID
    std::cout << "EleID MVA Fall17V2 iso wp90 pass? " << ele.electronID("mvaEleID-Fall17-iso-V2-wp90") << std::endl;
    //Run3 ID
    std::cout << "EleID MVA Winter22 iso wp90 pass? " << ele.electronID("mvaEleID-RunIIIWinter22-iso-V1-wp90") << std::endl;

  }
  
  //  tree->Fill();
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
EGMGenericNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
EGMGenericNtupler::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EGMGenericNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGMGenericNtupler);
