import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/tmp/adamwo/skimPV-Feb9.root'
    )
)

process.demo = cms.EDAnalyzer('BeamSpotFitPV',
  minNrVerticesForFit = cms.uint32(100),
  vertexCacheSize = cms.uint32(5000),
  beamspotTag = cms.InputTag("offlineBeamSpot"),
  vertexTag = cms.InputTag("offlinePrimaryVertices"),
  minVertexNdf = cms.double(4.),
  maxVertexNormChi2 = cms.double(10.),
  minVertexNTracks = cms.uint32(0),
  minVertexMeanWeight = cms.double(0.5),
  maxVertexR = cms.double(999),
  maxVertexZ = cms.double(999),
  errorScale = cms.double(0.9),
  nSigmaCut = cms.double(3.)
)


process.p = cms.Path(process.demo)
