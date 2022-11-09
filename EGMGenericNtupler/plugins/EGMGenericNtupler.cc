// -*- C++ -*-
// EGMGenericNtupler
// Author:  Swagata Mukherjee

// system include files
#include <memory>

// user include files
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
//#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
//#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
//#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
//#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "TTree.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"

class EGMGenericNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit EGMGenericNtupler(const edm::ParameterSet&);
  ~EGMGenericNtupler();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  edm::Service<TFileService> fs;
  TTree   *tree = fs->make<TTree>("EventTree", "EventData");

  //
  //  std::vector<float>  puTrue;
  //std::vector<float>  rho;
  //  std::vector<int>    gen_status;
  // std::vector<int>    gen_pdgid;
  //std::vector<float>  gen_pt;

  /// photon
  //  std::vector<float>  ptRecoPho_by_ptGenPho;
  //std::vector<float>  dR_recoPho_genPho;
  //std::vector<int>  pho_genmatched;
  //std::vector<std::vector<int> > pho_IDbits;
  //std::vector<int> pho_IDLoose;
  //std::vector<int> pho_IDMedium;
  //std::vector<int> pho_IDTight;
  //std::vector<double> pho_Pt;
  //  std::vector<double> pho_HoE;
  //std::vector<double> pho_SCrawE;
  //std::vector<double> pho_ChargedHadronWorstVtxIso;
  //std::vector<double> pho_ChargedHadronIso;
  //std::vector<double> pho_PhotonIso;
  //std::vector<double> pho_EtaWidth;
  //std::vector<double> pho_PhiWidth;
  //std::vector<double> pho_Phi;
  //std::vector<double> pho_ScEn;
  //std::vector<double> pho_ScEta;
  //std::vector<double> pho_SigmaIetaIeta;

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //  virtual  void setbit(UShort_t& x, UShort_t bit);
  
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
  //eleToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  phoToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("Photons")))
  //ootphoToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("ootPhotons"))),
  //puCollection_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupCollection"))),
  //genParticlesCollection_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc")))

{



  //  tree->Branch("puTrue_", &puTrue);
  //tree->Branch("rho_", &rho);
  //tree->Branch("gen_status_",&gen_status);
  //tree->Branch("gen_pdgid_",&gen_pdgid);
  // tree->Branch("gen_pt_",&gen_pt);

  ////photon
  /*
  tree->Branch("ptRecoPho_by_ptGenPho_",&ptRecoPho_by_ptGenPho);
  tree->Branch("dR_recoPho_genPho_",&dR_recoPho_genPho);
  tree->Branch("pho_genmatched_",&pho_genmatched);
  tree->Branch("pho_IDbits_", &pho_IDbits);
  tree->Branch("pho_IDLoose_", &pho_IDLoose);
  tree->Branch("pho_IDMedium_", &pho_IDMedium);
  tree->Branch("pho_IDTight_", &pho_IDTight);
  tree->Branch("pho_Pt_",&pho_Pt);
  tree->Branch("pho_R9_",&pho_R9);
  tree->Branch("pho_HoE_",&pho_HoE);
  tree->Branch("pho_SCrawE_",&pho_SCrawE);
  tree->Branch("pho_ChargedHadronWorstVtxIso_",&pho_ChargedHadronWorstVtxIso);
  tree->Branch("pho_ChargedHadronIso_",&pho_ChargedHadronIso);
  tree->Branch("pho_PhotonIso_",&pho_PhotonIso);
  tree->Branch("pho_EtaWidth_",&pho_EtaWidth);
  tree->Branch("pho_PhiWidth_",&pho_PhiWidth);
  tree->Branch("pho_Phi_",&pho_Phi);
  tree->Branch("pho_ScEn_",&pho_ScEn);
  tree->Branch("pho_ScEta_",&pho_ScEta);
  tree->Branch("pho_SigmaIetaIeta_",&pho_SigmaIetaIeta);
  */
}

EGMGenericNtupler::~EGMGenericNtupler()
{
}

//
// member functions
//

// ------------ method called for each event  ------------
void
EGMGenericNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  
  //
  //  puTrue.clear();
  // rho.clear();
  // gen_status.clear();
  // gen_pdgid.clear();
  // gen_pt.clear();

  //
  // photon
  //
  /*
  pho_genmatched.clear();
  dR_recoPho_genPho.clear();
  ptRecoPho_by_ptGenPho.clear();
  pho_IDbits.clear();
  pho_IDLoose.clear();
  pho_IDMedium.clear();
  pho_IDTight.clear();
  pho_Pt.clear();
  pho_R9.clear();
  pho_HoE.clear();
  pho_SCrawE.clear();
  pho_ChargedHadronWorstVtxIso.clear();
  pho_ChargedHadronIso.clear();
  pho_PhotonIso.clear();
  pho_EtaWidth.clear();
  pho_PhiWidth.clear();
  pho_Phi.clear();
  pho_ScEn.clear();
  pho_ScEta.clear();
  pho_SigmaIetaIeta.clear();
  */
  // edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  
  ////
  //  if ( ! iEvent.isRealData()) {
  //edm::Handle<std::vector<PileupSummaryInfo> > genPileupHandle;
  //iEvent.getByToken(puCollection_, genPileupHandle);
  //if (genPileupHandle.isValid()) {
  //  for (std::vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
  //	puTrue.push_back(pu->getTrueNumInteractions());
  //  }
  //}
    
  //iEvent.getByToken(genParticlesCollection_, genParticlesHandle);
    
  //   if (genParticlesHandle.isValid()) {
  //  for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
  //	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
  //	gen_status.push_back(p->status());
  //	gen_pdgid.push_back(p->pdgId());
  //	gen_pt.push_back(p->pt());
  //  }  
  //}
  //} // check real data
  
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  //  rho.push_back(*(rhoHandle.product()));
  std::cout << "rho = " << *(rhoHandle.product()) << std::endl;
  ////photon
  for(const auto& pho : iEvent.get(phoToken_) ) {

    // bool pass_tight = pho.photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
    //pho_IDbits.push_back({pho.userInt("cutBasedPhotonID-Fall17-94X-V2-loose"),
    //	  pho.userInt("cutBasedPhotonID-Fall17-94X-V2-medium"),
    //	  pho.userInt("cutBasedPhotonID-Fall17-94X-V2-tight")});
    
    //    int thispho_genmatched=0;
    // double thispho_min_dr=9999.9;
    // double thispho_ptR=9999.9;
    //if ( ! iEvent.isRealData()) {
      
    //if (genParticlesHandle.isValid()) {
    //	for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
    //	  const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	  
    //	  if ( (p->status()==1) && ((abs(p->pdgId())) == 22) )  {
	    
    //	    double this_dr=reco::deltaR(pho,*p);
    //	    if (this_dr<thispho_min_dr) {
    //	      thispho_min_dr=this_dr;
    //	      thispho_ptR=pho.pt()/p->pt();
    //	    }
    //	  }
    //	}
    //}
    //}
    
    //    if ( (thispho_min_dr<0.1) && (thispho_ptR>0.7) && (thispho_ptR<1.3) )  thispho_genmatched=1;
    
    // pho_genmatched.push_back(thispho_genmatched);
    //dR_recoPho_genPho.push_back(thispho_min_dr);
    //ptRecoPho_by_ptGenPho.push_back(thispho_ptR);
    //pho_Pt.push_back(pho.pt());
    //pho_Phi.push_back(pho.phi());
    //pho_ScEn.push_back(pho.superCluster()->energy());
    //pho_ScEta.push_back(pho.superCluster()->eta());
    //pho_SigmaIetaIeta.push_back(pho.full5x5_sigmaIetaIeta());
    //pho_HoE.push_back(pho.hadTowOverEm());
    //pho_R9.push_back(pho.full5x5_r9());
    //pho_EtaWidth.push_back(pho.superCluster()->etaWidth());
    //pho_PhiWidth.push_back(pho.superCluster()->phiWidth());
    //pho_PhotonIso.push_back(pho.photonIso());  
    //pho_ChargedHadronIso.push_back(pho.chargedHadronIso());  

    //    if  ( (pho.chargedHadronIso() < 0 ) || isnan(pho.chargedHadronIso()) )
    //{
	std::cout << "\n\n ** pho.chargedHadronIso() = " << pho.chargedHadronIso() << " pho.chargedHadronPFPVIso() = " << pho.chargedHadronPFPVIso() << std::endl;
	std::cout << "pho.pt() " << pho.pt() << " pho.photonIso() " << pho.photonIso() << " pho.superCluster()->eta() " << pho.superCluster()->eta() << std::endl;
	std::cout << "tightID " << pho.photonID("cutBasedPhotonID-Fall17-94X-V2-tight") << std::endl;
	std::cout << "looseID " << pho.photonID("cutBasedPhotonID-Fall17-94X-V2-loose") << std::endl;
	std::cout << "\n"; 
	//}
    //    pho_ChargedHadronWorstVtxIso.push_back(pho.chargedHadronWorstVtxIso());
    //  pho_SCrawE.push_back(pho.superCluster()->rawEnergy());    
   

          
  }

  
  tree->Fill();
  
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
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);


}

//define this as a plug-in
DEFINE_FWK_MODULE(EGMGenericNtupler);
