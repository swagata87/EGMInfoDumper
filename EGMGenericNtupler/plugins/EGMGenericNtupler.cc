// -*- C++ -*-
// EGMGenericNtupler
// Author:  Swagata Mukherjee

// system include files
#include <memory>

// user include files
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "DataFormats/TrackReco/interface/TrackExtra.h"
//#include "DataFormats/TrackReco/interface/HitPattern.h"
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
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/plugins/CaloTopologyBuilder.h"
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
  
  edm::Service<TFileService> fs;
  TTree   *tree = fs->make<TTree>("EventTree", "EventData");

  std::vector<std::vector<int> > ele_IDbits;
  std::vector<bool>  ele_IDVeto;
  std::vector<bool>  ele_IDLoose;
  std::vector<bool>  ele_IDMedium;
  std::vector<bool>  ele_IDTight;
  std::vector<bool>  ele_IDMVAiso90;
  std::vector<int>  ele_isEB;
  std::vector<int>  ele_isEE;
  std::vector<double>  ele_SigmaIetaIeta;
  std::vector<double>  ele_SigmaIphiIphi;
  std::vector<float>  ele_R9;
  std::vector<float>  ele_HoE;
  std::vector<float>  ele_HoEfull5x5;
  std::vector<float>  ele_ScEta;
  std::vector<float>  ele_ScPhi;
  std::vector<float>  ele_Eta;
  std::vector<float>  ele_Phi;
  std::vector<float>  ele_ScEn;
  std::vector<float>  ele_Pt;
  std::vector<bool>  ele_isEcalDriven;
  std::vector<float>  ele_fbrem;
  //  std::vector<float>  ele_missingHit;
  //  std::vector<float>  ele_pixelLayersTotallyOffOrBad;
  std::vector<float>  ele_psEorawE;
  std::vector<float>  ele_1oEm1op;
  std::vector<float>  ele_eSCoPout;
  std::vector<float>  ele_eSCoP;
   std::vector<float>  ele_convVtxFitProb;
    std::vector<float>  ele_gsfTrackChi2;
  ///  std::vector<float>  ele_nHit;
  std::vector<float>  ele_deltaEtaSuperClusterTrackAtVtx;
  std::vector<float>  ele_deltaPhiSuperClusterTrackAtVtx;
  std::vector<float>  ele_deltaEtaSeedClusterTrackAtCalo;
  std::vector<float>  ele_PFChIso;
  std::vector<float>  ele_PFPhoIso;
  std::vector<float>  ele_PFNeuIso;
  std::vector<float>  ele_etaWidth;
  std::vector<float>  ele_phiWidth;
  std::vector<float>  ele_closestCtfTrackNLayers;
  std::vector<float>  ele_closestCtfTrackNormChi2;
  std::vector<float>  ele_1mE1x5oE5x5;
  std::vector<float>  ele_gsfTrack_vz;
  std::vector<float>  ele_gsfTrack_outerZ;
  std::vector<float>  ele_gsfTrack_eta;
  std::vector<float>  ptRecoEle_by_ptGenEle;
  std::vector<float>  dR_recoEle_genEle;
  std::vector<int>  ele_genmatched;
  //
  std::vector<float>  puTrue;
  std::vector<float>  rho;
  std::vector<int>    gen_status;
  std::vector<int>    gen_pdgid;
  std::vector<float>  gen_pt;

  /// photon
  std::vector<float>  ptRecoPho_by_ptGenPho;
  std::vector<float>  dR_recoPho_genPho;
  std::vector<int>  pho_genmatched;
  std::vector<std::vector<int> > pho_IDbits;
  std::vector<int> pho_IDLoose;
  std::vector<int> pho_IDMedium;
  std::vector<int> pho_IDTight;
  std::vector<double> pho_Pt;
  std::vector<double> pho_R9;
  std::vector<double> pho_HoE;
  std::vector<double> pho_SCrawE;
  std::vector<double> pho_ChargedHadronWorstVtxIso;
  std::vector<double> pho_ChargedHadronIso;
  std::vector<double> pho_PhotonIso;
  std::vector<double> pho_EtaWidth;
  std::vector<double> pho_PhiWidth;
  std::vector<double> pho_Phi;
  std::vector<double> pho_ScEn;
  std::vector<double> pho_ScEta;
  std::vector<double> pho_SigmaIetaIeta;

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //  virtual  void setbit(UShort_t& x, UShort_t bit);
  
  // ----------member data ---------------------------

  edm::EDGetTokenT<double> rhoLabel_;
  //  edm::EDGetTokenT<edm::View<reco::GsfElectron> > eleToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > eleToken_;
  //edm::EDGetTokenT<edm::View<reco::Photon> > phoToken_;
  edm::EDGetTokenT<edm::View<pat::Photon> > phoToken_;
  edm::EDGetTokenT<edm::View<pat::Photon> > ootphoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> >     genParticlesCollection_;

};

EGMGenericNtupler::EGMGenericNtupler(const edm::ParameterSet& iConfig)
  :
  rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"))),
  eleToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  phoToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("Photons"))),
  ootphoToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("ootPhotons"))),
  puCollection_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupCollection"))),
  genParticlesCollection_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleSrc")))

{

  tree->Branch("ele_IDbits_", &ele_IDbits);
  tree->Branch("ele_IDVeto_", &ele_IDVeto);
  tree->Branch("ele_IDLoose_", &ele_IDLoose);
  tree->Branch("ele_IDMedium_", &ele_IDMedium);
  tree->Branch("ele_IDTight_", &ele_IDTight);
  tree->Branch("ele_IDMVAiso90_", &ele_IDMVAiso90);
  //
  tree->Branch("ele_isEB_",&ele_isEB);
  tree->Branch("ele_isEE_",&ele_isEE);
  tree->Branch("ele_genmatched_",&ele_genmatched);
  tree->Branch("ele_SigmaIetaIeta_",&ele_SigmaIetaIeta);
  tree->Branch("ele_SigmaIphiIphi_",&ele_SigmaIphiIphi);
  tree->Branch("ele_R9_",&ele_R9);
  tree->Branch("ele_HoE_",&ele_HoE);
  tree->Branch("ele_HoEfull5x5_",&ele_HoEfull5x5);
  tree->Branch("ele_ScEta_",&ele_ScEta);
  tree->Branch("ele_ScPhi_",&ele_ScPhi);
  tree->Branch("ele_Eta_",&ele_Eta);
  tree->Branch("ele_Phi_",&ele_Phi);
  tree->Branch("ele_ScEn_",&ele_ScEn);
  tree->Branch("ele_Pt_",&ele_Pt);
  tree->Branch("ele_isEcalDriven_",&ele_isEcalDriven);
  tree->Branch("ele_fbrem_",&ele_fbrem);
  // tree->Branch("ele_missingHit_",&ele_missingHit);
  //tree->Branch("ele_pixelLayersTotallyOffOrBad_",&ele_pixelLayersTotallyOffOrBad);
  tree->Branch("ele_psEorawE_",&ele_psEorawE);
  tree->Branch("ele_1oEm1op_",&ele_1oEm1op);
  tree->Branch("ele_eSCoP_",&ele_eSCoP);
  tree->Branch("ele_eSCoPout_",&ele_eSCoPout);
  tree->Branch("ele_convVtxFitProb_",&ele_convVtxFitProb);
  tree->Branch("ele_gsfTrackChi2_",&ele_gsfTrackChi2);
  // tree->Branch("ele_nHit_",&ele_nHit);
  tree->Branch("ele_gsfTrack_vz_",&ele_gsfTrack_vz);
  tree->Branch("ele_gsfTrack_outerZ_",&ele_gsfTrack_outerZ);
  tree->Branch("ele_gsfTrack_eta_",&ele_gsfTrack_eta);
  tree->Branch("ele_deltaEtaSuperClusterTrackAtVtx_",&ele_deltaEtaSuperClusterTrackAtVtx);
  tree->Branch("ele_deltaPhiSuperClusterTrackAtVtx_",&ele_deltaPhiSuperClusterTrackAtVtx);
  tree->Branch("ele_deltaEtaSeedClusterTrackAtCalo_",&ele_deltaEtaSeedClusterTrackAtCalo);
  tree->Branch("ele_PFChIso_",&ele_PFChIso);
  tree->Branch("ele_PFPhoIso_",&ele_PFPhoIso);
  tree->Branch("ele_PFNeuIso_",&ele_PFNeuIso);
  tree->Branch("ele_etaWidth_",&ele_etaWidth);
  tree->Branch("ele_phiWidth_",&ele_phiWidth);
  tree->Branch("ele_closestCtfTrackNLayers_",&ele_closestCtfTrackNLayers);
  tree->Branch("ele_closestCtfTrackNormChi2_",&ele_closestCtfTrackNormChi2);
  tree->Branch("ele_1mE1x5oE5x5_",&ele_1mE1x5oE5x5);
  tree->Branch("ptRecoEle_by_ptGenEle_",&ptRecoEle_by_ptGenEle);
  tree->Branch("dR_recoEle_genEle_",&dR_recoEle_genEle);

  tree->Branch("puTrue_", &puTrue);
  tree->Branch("rho_", &rho);
  tree->Branch("gen_status_",&gen_status);
  tree->Branch("gen_pdgid_",&gen_pdgid);
  tree->Branch("gen_pt_",&gen_pt);

  ////photon
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
  
  std::cout << "\n \n --New Event-- \n" ;
  using namespace edm;
  
  ele_IDVeto.clear();
  ele_IDLoose.clear();
  ele_IDMedium.clear();
  ele_IDTight.clear();
  ele_IDMVAiso90.clear();
  ele_IDbits.clear();
  ele_isEB.clear();
  ele_isEE.clear();
  ele_genmatched.clear();
  ele_SigmaIetaIeta.clear();
  ele_SigmaIphiIphi.clear();
  ele_R9.clear();
  ele_HoE.clear();
  ele_HoEfull5x5.clear();
  ele_ScEta.clear();
  ele_ScPhi.clear();
  ele_Eta.clear();
  ele_Phi.clear();
  ele_ScEn.clear();
  ele_Pt.clear();
  ele_isEcalDriven.clear();
  ele_fbrem.clear();
  //ele_missingHit.clear();
  //ele_pixelLayersTotallyOffOrBad.clear();
  ele_psEorawE.clear();
  ele_1oEm1op.clear();
  ele_eSCoP.clear();
  ele_eSCoPout.clear();
  ele_convVtxFitProb.clear();
  ele_gsfTrackChi2.clear();
  //ele_nHit.clear();
  ele_deltaEtaSuperClusterTrackAtVtx.clear();
  ele_deltaPhiSuperClusterTrackAtVtx.clear();
  ele_deltaEtaSeedClusterTrackAtCalo.clear();
  ele_PFChIso.clear();
  ele_PFPhoIso.clear();
  ele_PFNeuIso.clear();
  ele_etaWidth.clear();
  ele_phiWidth.clear();
  ele_gsfTrack_vz.clear();
  ele_gsfTrack_outerZ.clear();
  ele_gsfTrack_eta.clear();
  ele_closestCtfTrackNLayers.clear();
  ele_closestCtfTrackNormChi2.clear();
  ele_1mE1x5oE5x5.clear();
  ptRecoEle_by_ptGenEle.clear();
  dR_recoEle_genEle.clear();
 
  //
  puTrue.clear();
  rho.clear();
  gen_status.clear();
  gen_pdgid.clear();
  gen_pt.clear();

  //
  // photon
  //
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

  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
  
  ////
  if ( ! iEvent.isRealData()) {
    edm::Handle<std::vector<PileupSummaryInfo> > genPileupHandle;
    iEvent.getByToken(puCollection_, genPileupHandle);
    if (genPileupHandle.isValid()) {
      for (std::vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
	puTrue.push_back(pu->getTrueNumInteractions());
      }
    }
    
    iEvent.getByToken(genParticlesCollection_, genParticlesHandle);
    
    if (genParticlesHandle.isValid()) {
      for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	gen_status.push_back(p->status());
	gen_pdgid.push_back(p->pdgId());
	gen_pt.push_back(p->pt());
      }  
    }
  } // check real data
  
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  rho.push_back(*(rhoHandle.product()));

  ///oot pho
  //std::cout << "will enter oot now " << std::endl;
  for(const auto& ootpho : iEvent.get(ootphoToken_) ) {
    std::cout << "\n new oot photon " << std::endl;
    std::cout << "var " << ootpho.pt() << " " << ootpho.chargedHadronIso() 
	      << " " << ootpho.photonIso() << " " << ootpho.neutralHadronIso() 
	      << " UF " << ootpho.userFloat("phoChargedIsolation") << std::endl;
  }


  ////photon
  for(const auto& pho : iEvent.get(phoToken_) ) {
    //std::cout << "\n new photon " << std::endl;
    
    //    bool pass_loose = pho.photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
    // bool pass_medium = pho.photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
    // bool pass_tight = pho.photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
    
    //  pho_IDLoose.push_back(pass_loose);
    // pho_IDMedium.push_back(pass_medium);
    // pho_IDTight.push_back(pass_tight);
    
    //pho_IDbits.push_back({pho.userInt("cutBasedPhotonID-Fall17-94X-V2-loose"),
    //	  pho.userInt("cutBasedPhotonID-Fall17-94X-V2-medium"),
    //	  pho.userInt("cutBasedPhotonID-Fall17-94X-V2-tight")});
    
    int thispho_genmatched=0;
    double thispho_min_dr=9999.9;
    double thispho_ptR=9999.9;
    if ( ! iEvent.isRealData()) {
      
      if (genParticlesHandle.isValid()) {
	for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	  const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	  
	  if ( (p->status()==1) && ((abs(p->pdgId())) == 22) )  {
	    
	    double this_dr=reco::deltaR(pho,*p);
	    if (this_dr<thispho_min_dr) {
	      thispho_min_dr=this_dr;
	      thispho_ptR=pho.pt()/p->pt();
	    }
	  }
	}
      }
    }
    
    if ( (thispho_min_dr<0.1) && (thispho_ptR>0.7) && (thispho_ptR<1.3) )  thispho_genmatched=1;
    
    pho_genmatched.push_back(thispho_genmatched);
    dR_recoPho_genPho.push_back(thispho_min_dr);
    ptRecoPho_by_ptGenPho.push_back(thispho_ptR);
    pho_Pt.push_back(pho.pt());
    pho_Phi.push_back(pho.phi());
    pho_ScEn.push_back(pho.superCluster()->energy());
    pho_ScEta.push_back(pho.superCluster()->eta());
    pho_SigmaIetaIeta.push_back(pho.full5x5_sigmaIetaIeta());
    pho_HoE.push_back(pho.hadTowOverEm());
    pho_R9.push_back(pho.full5x5_r9());
    pho_EtaWidth.push_back(pho.superCluster()->etaWidth());
    pho_PhiWidth.push_back(pho.superCluster()->phiWidth());
    pho_PhotonIso.push_back(pho.photonIso());  
    pho_ChargedHadronIso.push_back(pho.chargedHadronIso());  
    pho_ChargedHadronWorstVtxIso.push_back(pho.chargedHadronWorstVtxIso());
    pho_SCrawE.push_back(pho.superCluster()->rawEnergy());    
   
    // older versions:
    // phoPhotonIso.push_back(pho.userFloat("phoPhotonIsolation"));  
    //phoChargedHadronIso.push_back(pho.userFloat("phoChargedIsolation"));  
    //phoChargedHadronWorstVtxIso.push_back(pho.userFloat("phoWorstChargedIsolation"));
          
  }

  for(const auto& ele : iEvent.get(eleToken_) ) {
    // std::cout << "\n ---/// New Electron .... " << std::endl;
    
    //    bool isPassVeto   = ele.electronID("cutBasedElectronID-Fall17-94X-V2-veto");
    // bool isPassLoose   = ele.electronID("cutBasedElectronID-Fall17-94X-V2-loose");
    //  bool isPassMedium   = ele.electronID("cutBasedElectronID-Fall17-94X-V2-medium");
    // bool isPassTight   = ele.electronID("cutBasedElectronID-Fall17-94X-V2-tight");
    // bool isPassMVAiso90 = ele.electronID("mvaEleID-Fall17-iso-V2-wp90");

    // ele_IDVeto.push_back(isPassVeto);
    // ele_IDLoose.push_back(isPassLoose);
    // ele_IDMedium.push_back(isPassMedium);
    // ele_IDTight.push_back(isPassTight);
    // ele_IDMVAiso90.push_back(isPassMVAiso90);

    // ele_IDbits.push_back({          
    //	  ele.userInt("cutBasedElectronID-Fall17-94X-V2-veto"),
    //	    ele.userInt("cutBasedElectronID-Fall17-94X-V2-loose"),
    //	    ele.userInt("cutBasedElectronID-Fall17-94X-V2-medium"),
    ///	    ele.userInt("cutBasedElectronID-Fall17-94X-V2-tight"),
    //	    ele.userInt("heepElectronID-HEEPV70")});
      
      int genmatched=0;
      double min_dr=9999.9;
      double ptR=9999.9;
   
      if ( ! iEvent.isRealData()) {
	
	if (genParticlesHandle.isValid()) {
	  for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	    const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	    if ( (p->status()==1) &&  ((abs(p->pdgId())) == 11) ) {
	      double this_dr=reco::deltaR(ele,*p);
	      if (this_dr<min_dr) {
		min_dr=this_dr;
		ptR=ele.pt()/p->pt();
	      }
	    }
	  }
	}
      }

    
      if ( (min_dr<0.1) && (ptR>0.7) && (ptR<1.3) )  genmatched=1;

      ele_genmatched.push_back(genmatched);
      dR_recoEle_genEle.push_back(min_dr);
      ptRecoEle_by_ptGenEle.push_back(ptR);
      ele_isEB.push_back(ele.isEB()) ;
      ele_isEE.push_back(ele.isEE()) ;
      ele_SigmaIetaIeta.push_back(ele.full5x5_sigmaIetaIeta());
      ele_SigmaIphiIphi.push_back(ele.full5x5_sigmaIphiIphi());
      ele_R9.push_back(ele.full5x5_r9());
      ele_HoE.push_back(ele.hcalOverEcal());
      ele_HoEfull5x5.push_back(ele.full5x5_hcalOverEcal());
      ele_ScEta.push_back(ele.superCluster()->eta());
      ele_ScPhi.push_back(ele.superCluster()->phi());
      ele_Eta.push_back(ele.eta());
      ele_Phi.push_back(ele.phi());
      ele_ScEn.push_back(ele.superCluster()->energy());
      ele_Pt.push_back(ele.pt());
      ele_isEcalDriven.push_back(ele.ecalDriven());
      ele_1mE1x5oE5x5.push_back(1-(ele.full5x5_e1x5()/ele.full5x5_e5x5()));      
      ele_etaWidth.push_back(ele.superCluster()->etaWidth());
      ele_phiWidth.push_back(ele.superCluster()->phiWidth());
      ele_deltaEtaSuperClusterTrackAtVtx.push_back(ele.deltaEtaSuperClusterTrackAtVtx());
      ele_deltaPhiSuperClusterTrackAtVtx.push_back(ele.deltaPhiSuperClusterTrackAtVtx());
      ele_deltaEtaSeedClusterTrackAtCalo.push_back(ele.deltaEtaSeedClusterTrackAtCalo());
      ele_fbrem.push_back(ele.fbrem());
      //      ele_missingHit.push_back(ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      //      ele_pixelLayersTotallyOffOrBad.push_back(ele.gsfTrack()->hitPattern().pixelLayersTotallyOffOrBad());
      //      ele_nHit.push_back(ele.gsfTrack()->hitPattern().trackerLayersWithMeasurement());
      // ele_gsfTrackChi2.push_back(ele.gsfTrack()->normalizedChi2());
      // ele_gsfTrack_vz.push_back(ele.gsfTrack()->vz());
      //    ele_gsfTrack_outerZ.push_back(ele.gsfTrack()->outerZ());
      //  ele_gsfTrack_eta.push_back(ele.gsfTrack()->eta());
      ele_eSCoP.push_back(ele.eSuperClusterOverP());
      ele_eSCoPout.push_back(ele.eEleClusterOverPout());
      ele_1oEm1op.push_back( (1.0/ele.ecalEnergy()) - (1.0/ele.trackMomentumAtVtx().r()) );
      ele_psEorawE.push_back( ele.superCluster()->preshowerEnergy()/ele.superCluster()->rawEnergy() );

      reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
      ele_PFChIso.push_back(pfIso.sumChargedHadronPt);
      ele_PFPhoIso.push_back(pfIso.sumPhotonEt);
      ele_PFNeuIso.push_back(pfIso.sumNeutralHadronEt);

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

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EGMGenericNtupler);
