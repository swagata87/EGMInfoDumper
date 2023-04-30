import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 50

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')   

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                '/store/data/Run2018A/ParkingBPH1/AOD/20Jun2021_UL2018-v1/260004/E9E566F2-3418-5E44-B11B-99E4947B0F6F.root'
                                )
)

process.demo = cms.EDAnalyzer('EGMGenericNtupler',
       electrons = cms.InputTag('gedGsfElectrons'),
       Photons       = cms.InputTag("gedPhotons"),
       mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values"),
)

#Set up Egamma post-reco tool 
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq,makeEgammaPATWithUserData

setupEgammaPostRecoSeq(process,
                       era='2018-UL',
                       isMiniAOD=False,
                       runVID=True
)  

#Run Egamma post reco tool
process.p = cms.Path( process.egammaPostRecoSeq * process.demo )
