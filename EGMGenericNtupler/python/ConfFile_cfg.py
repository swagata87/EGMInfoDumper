import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#Right now for Run3 samples EgammaPostRecoTools is not needed. But it will be needed later. 
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

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
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2021_realistic_v3', '')     # 2021 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')     # 2018 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018Rereco_v1', '')     # 2018 data
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v6', '')     # 2017 data 

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1),
#    moduleMemorySummary = cms.untracked.bool(True),
#)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

setupEgammaPostRecoSeq(process,
#                       era='2017-UL',
                       era='2018-Prompt',
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']
)  

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#                                '/store/mc/Run3Summer19MiniAOD/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/MINIAODSIM/2021Scenario_106X_mcRun3_2021_realistic_v3-v2/130000/00DF0005-F507-2C4B-BF8B-C46342D7194E.root'
#                                '/store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/80000/E7AE77A2-59DD-6C4B-AF91-8263DA41EFD8.root'
#                                '/store/mc/Run3Summer19MiniAOD/DYToEE_M-50_NNPDF31_TuneCP5_14TeV-powheg-pythia8/MINIAODSIM/2023Scenario_106X_mcRun3_2023_realistic_v3-v2/260000/FE3A7D24-F46E-8744-B92B-F6115FD395A0.root'
#                                '/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/10000/07184556-FFC3-A94A-8124-9CBD85FCD2C9.root',
#                                '/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/50000/F5215AC7-4F83-3C42-A188-46198D345570.root',
#                                '/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/50000/EDBE4672-2F34-154A-9A46-30E4146A7851.root',
#                                '/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/50000/EDA7C170-3345-9D4B-917E-16EA730D6C2B.root',
#                                '/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/50000/E9BF7ABB-46DF-9849-938F-6A1C9B5A1E8E.root',
#                                '/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/50000/E8108FFB-A1F2-7A47-92E6-2667EA6582A3.root',
#                                '/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/50000/E50D3AA3-B57B-BB4E-83FC-E2D4D5C0457D.root',
#                                '/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/50000/E1434D5D-4400-9346-8204-588CEC1BEDA5.root'
#                                '/store/mc/Run3Winter20DRMiniAOD/QCD_Pt-15to7000_TuneCP5_Flat_14TeV_pythia8/MINIAODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v1/130000/FB300C18-FEEA-7341-9B96-1A835CFC65EB.root'
                                '/store/data/Run2017C/SinglePhoton/MINIAOD/31Mar2018-v1/00000/80157357-DD37-E811-A5F2-0025907D244A.root'
                                )
)

process.demo = cms.EDAnalyzer('EGMGenericNtupler',
       rhoLabel  = cms.InputTag("fixedGridRhoFastjetAll"),
       electrons = cms.InputTag('slimmedElectrons'),
       pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
       genParticleSrc       = cms.InputTag("prunedGenParticles"),
       Photons       = cms.InputTag("slimmedPhotons"),
       ootPhotons       = cms.InputTag("slimmedOOTPhotons"),
                              )

#process.TFileService = cms.Service("TFileService", fileName = cms.string('genericEGMntuple_QCD_2021.root'))
process.TFileService = cms.Service("TFileService", fileName = cms.string('genericEGMntuple_data_2017.root'))

#process.p = cms.Path(process.egammaPostRecoSeq * process.demo)
process.p = cms.Path(process.demo)
