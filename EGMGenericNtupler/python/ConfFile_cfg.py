import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#Right now for Run3 samples EgammaPostRecoTools is not needed. But it will be needed later. 
#from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2024_realistic_v4', '')     # 2024 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2023_realistic_v3', '')     # 2023 MC
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2021_realistic_v3', '')     # 2021 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')     # 2018 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018Rereco_v1', '')     # 2018 data

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1),
#    moduleMemorySummary = cms.untracked.bool(True),
#)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

#setupEgammaPostRecoSeq(process,
#                       era='2017-UL',
#                       era='2018-Prompt',
#                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']
#)  

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                '/store/mc/Run3Summer19MiniAOD/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/MINIAODSIM/2021Scenario_106X_mcRun3_2021_realistic_v3-v2/130000/00DF0005-F507-2C4B-BF8B-C46342D7194E.root'
                                #'/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/E7AE77A2-59DD-6C4B-AF91-8263DA41EFD8.root'
#                                '/store/mc/Run3Summer19MiniAOD/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/2023Scenario_106X_mcRun3_2023_realistic_v3-v2/260000/FE3A7D24-F46E-8744-B92B-F6115FD395A0.root'
                                )
)

process.demo = cms.EDAnalyzer('EGMGenericNtupler',
       rhoLabel  = cms.InputTag("fixedGridRhoFastjetAll"),
       electrons = cms.InputTag('slimmedElectrons'),
       pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
       genParticleSrc       = cms.InputTag("prunedGenParticles"),
       Photons       = cms.InputTag("slimmedPhotons"),
                              )
process.TFileService = cms.Service("TFileService", fileName = cms.string('genericEGMntuple.root'))

#process.p = cms.Path(process.egammaPostRecoSeq * process.demo)
process.p = cms.Path(process.demo)
