import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '122X_mcRun3_2021_realistic_v9', '')     # 2022 MC
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                '/store/mc/Run3Winter22MiniAOD/GJet_Pt-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8/MINIAODSIM/FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/2430000/b30da4dc-1771-491e-b4ce-f424c076aab5.root',
                                )
)

process.demo = cms.EDAnalyzer('EGMGenericNtupler',
       rhoLabel  = cms.InputTag("fixedGridRhoFastjetAll"),
       Photons       = cms.InputTag("slimmedPhotons"),
       effAreas_charged = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt"),
       effAreas_neutral = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt"),
       effAreas_photon = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt"),
                              )

process.TFileService = cms.Service("TFileService", fileName = cms.string('genericEGMntuple_mc.root'))
process.p = cms.Path(process.demo)
