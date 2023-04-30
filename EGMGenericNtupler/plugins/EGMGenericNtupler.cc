// // -*- C++ -*-
// // EGMGenericNtupler
// // Author:  Swagata Mukherjee

// // system include files
// #include <memory>

// // user include files
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
  
private:
   virtual void beginJob() override;
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endJob() override;

//   // ----------member data ---------------------------
  edm::EDGetToken electronsToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
  
};

EGMGenericNtupler::EGMGenericNtupler(const edm::ParameterSet& iConfig)
  :
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap")))
 {
   electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"));
 }

 EGMGenericNtupler::~EGMGenericNtupler()
 {}

void
EGMGenericNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "\n \n --New Event-- \n" ;
  using namespace edm;
  edm::Handle<edm::ValueMap<float> > mvaValues;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsToken_, electrons);
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);
    std::cout << (*mvaValues)[el] << std::endl;
  }
  
  
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
