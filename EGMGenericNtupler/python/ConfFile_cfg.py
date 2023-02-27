import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_mcRun3_2022_realistic_postEE_v1', '')   

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                '/store/mc/Run3Summer22EEMiniAODv3/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/MINIAODSIM/forPOG_124X_mcRun3_2022_realistic_postEE_v1-v3/2810000/02d4983f-f18f-4ca1-923c-7496924691cc.root'
                                )
)

process.demo = cms.EDAnalyzer('EGMGenericNtupler',
       electrons = cms.InputTag('slimmedElectrons'),
       Photons       = cms.InputTag("slimmedPhotons"),
                              )

#Set up Egamma post-reco tool for Run3
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       era='run3',
                       runEnergyCorrections=False,
                       runVID=True,
)  

#Run Egamma post reco tool
process.p = cms.Path(process.egammaPostRecoSeq * process.demo)
