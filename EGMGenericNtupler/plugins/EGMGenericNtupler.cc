// -*- C++ -*-
// EGMGenericNtupler
// Author:  Swagata Mukherjee

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "TTree.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"

class EGMGenericNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit EGMGenericNtupler(const edm::ParameterSet&);
  ~EGMGenericNtupler();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  edm::Service<TFileService> fs;
  TTree   *tree = fs->make<TTree>("EventTree", "EventData");
  //std::vector<float>  rho;
  //std::vector<double> pho_Pt;

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  double getEAcorrectedVar(EffectiveAreas thisEA, double thisAbsEta, double thisRho, double thisVar);

  // ----------member data ---------------------------
  EffectiveAreas _effectiveAreas_charged;
  EffectiveAreas _effectiveAreas_neutral;
  EffectiveAreas _effectiveAreas_photon;
  edm::EDGetTokenT<double> rhoLabel_;
  edm::EDGetTokenT<edm::View<pat::Photon> > phoToken_;
};

EGMGenericNtupler::EGMGenericNtupler(const edm::ParameterSet& iConfig)
  :
  _effectiveAreas_charged((iConfig.getParameter<edm::FileInPath>("effAreas_charged")).fullPath()),
  _effectiveAreas_neutral((iConfig.getParameter<edm::FileInPath>("effAreas_neutral")).fullPath()),
  _effectiveAreas_photon((iConfig.getParameter<edm::FileInPath>("effAreas_photon")).fullPath()),
  rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"))),
  phoToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("Photons")))
{
  //tree->Branch("rho_", &rho);
  //  tree->Branch("pho_Pt_",&pho_Pt);
}

EGMGenericNtupler::~EGMGenericNtupler()
{
}

// ------------ method called for each event  ------------
void
EGMGenericNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
    // rho.clear();
  //  pho_Pt.clear();

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  //  rho.push_back(*(rhoHandle.product()));
  double rho = rhoHandle.isValid() ? *rhoHandle : 0;
  std::cout << "rho " << rho << std::endl;

  ////photon
  for(const auto& pho : iEvent.get(phoToken_) ) {

    double absEta = std::abs(pho.superCluster()->eta());
    std::cout << "photon pT/ absEta " << pho.pt() << " " << absEta << std::endl;

    double chargedHadronIsoWithEA = getEAcorrectedVar(_effectiveAreas_charged, absEta, rho, pho.chargedHadronIso());
    double neutralHadronIsoWithEA = getEAcorrectedVar(_effectiveAreas_neutral, absEta, rho, pho.neutralHadronIso());
    double photonIsoWithEA        = getEAcorrectedVar(_effectiveAreas_photon,  absEta, rho, pho.photonIso());

    std::cout << "charged/neutral/photon isolation (PU corrected) " << 
      chargedHadronIsoWithEA << " " << neutralHadronIsoWithEA << " " << photonIsoWithEA << std::endl;
    std::cout << "H/E and SigmaIEtaIEta " << pho.hcalOverEcalBc() << " " << pho.full5x5_sigmaIetaIeta() << std::endl;
    std::cout << "\n"; 
  }

  tree->Fill();
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

//function for PU correction
double EGMGenericNtupler::getEAcorrectedVar(EffectiveAreas thisEA, double thisAbsEta, double thisRho, double thisVar) {
  double valEA = thisEA.getEffectiveArea(thisAbsEta);
  double thisVar_EAcorr = std::max(0.0, thisVar - thisRho*valEA);
  return thisVar_EAcorr;
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
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGMGenericNtupler);
