import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/caf/user/sdevissc/QCD_RECODEBUG_TEST.root'
    )
)

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOEventContent.outputCommands,
    fileName = cms.untracked.string('file:out_inRECO.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('RECOfromRECO'),
        dataTier = cms.untracked.string('AODSIM')
    )
)


process.demo = cms.EDAnalyzer('SeedEfficiency',
 genParticleTag=cms.InputTag("genParticles"),
 gedgsfElectronTag= cms.InputTag("gedGsfElectrons"),
 GenPTag=cms.InputTag("genParticles"),
 JetsTag=cms.InputTag("ak5PFJets"),
 PVerTag=cms.InputTag("offlinePrimaryVertices"),
 simtracksTag = cms.InputTag("g4SimHits"),
 tracksTag = cms.InputTag("electronGsfTracks"),
 TrackingParticleTag=cms.InputTag("mix","MergedTrackTruth"),
 minHits = cms.uint32(8),
 UseRECO=cms.bool(True),
 Verbose=cms.bool(True)

)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_72_V3', '')


import SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi
from SimTracker.TrackerHitAssociation.clusterTpAssociationProducer_cfi import *
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.useClusterTPAssociation = cms.bool(True)
process.load("SimTracker.TrackerHitAssociation.clusterTpAssociationProducer_cfi")

# Dominik workaround
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.MixingModule.trackingTruthProducerSelection_cfi")
process.trackingParticles.simHitCollections = cms.PSet( )
process.mix.playback = cms.untracked.bool(True)
process.mix.digitizers = cms.PSet(
    mergedtruth = cms.PSet(process.trackingParticles)
)
for a in process.aliases: delattr(process, a)
process.reconstruction_step = cms.Path(process.reconstruction_fromRECO)
process.load("RecoParticleFlow.PFTracking.particleFlowTrack_cff")
process.make_pftracks = cms.Path(process.pfTrackingGlobalReco)


process.p = cms.Path(process.demo)
process.schedule = cms.Schedule(process.reconstruction_step,process.make_pftracks,process.p)

