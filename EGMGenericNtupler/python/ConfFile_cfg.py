import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

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

process.GlobalTag = GlobalTag(process.GlobalTag, '122X_mcRun3_2021_realistic_v9', '')     # 2022 MC



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(

#                                '/store/mc/Run3Winter22MiniAOD/GJet_Pt-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8/MINIAODSIM/FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/2430000/ad93b873-e478-4f1d-82ea-2962217328a2.root'
                                '/store/mc/Run3Winter22MiniAOD/GJet_Pt-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8/MINIAODSIM/FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/2430000/b30da4dc-1771-491e-b4ce-f424c076aab5.root',
#                                '/store/mc/Run3Winter22MiniAOD/GJet_Pt-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8/MINIAODSIM/FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/2430000/b52ef7ee-f58e-4a94-95fa-dd980c814966.root',

                                )
)

#process.source.eventsToProcess = cms.untracked.VEventRange('1:500-1:700')

process.demo = cms.EDAnalyzer('EGMGenericNtupler',
       rhoLabel  = cms.InputTag("fixedGridRhoFastjetAll"),
       #electrons = cms.InputTag('slimmedElectrons'),
       #pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
       #genParticleSrc       = cms.InputTag("prunedGenParticles"),
       Photons       = cms.InputTag("slimmedPhotons"),
       effAreas_charged = cms.FileInPath("/afs/cern.ch/user/s/swmukher/work/spoonFeedArrogantUsers/CMSSW_12_4_6/src/EGMInfoDumper/EGMGenericNtupler/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt"),
       effAreas_neutral = cms.FileInPath("/afs/cern.ch/user/s/swmukher/work/spoonFeedArrogantUsers/CMSSW_12_4_6/src/EGMInfoDumper/EGMGenericNtupler/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt"),
       effAreas_photon = cms.FileInPath("/afs/cern.ch/user/s/swmukher/work/spoonFeedArrogantUsers/CMSSW_12_4_6/src/EGMInfoDumper/EGMGenericNtupler/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt"),
       #ootPhotons       = cms.InputTag("slimmedOOTPhotons"),
                              )

process.TFileService = cms.Service("TFileService", fileName = cms.string('genericEGMntuple_mc.root'))

process.p = cms.Path(process.demo)
