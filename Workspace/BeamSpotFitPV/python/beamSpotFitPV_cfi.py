import FWCore.ParameterSet.Config as cms

beamSpotFitPV = cms.EDAnalyzer('BeamSpotFitPV',
                               minNrVerticesForFit = cms.uint32(100),
                               vertexCacheSize = cms.uint32(5000),
                               beamspotTag = cms.InputTag("offlineBeamSpot"),
                               vertexTag = cms.InputTag("offlinePrimaryVertices"),
                               minVertexNdf = cms.double(4.),
                               maxVertexNormChi2 = cms.double(10.),
                               minVertexNTracks = cms.uint32(0),
                               minVertexMeanWeight = cms.double(0.5),
                               maxVertexR = cms.double(2),
                               maxVertexZ = cms.double(30),
                               errorScale = cms.double(0.9),
                               nSigmaCut = cms.double(5.),
                               assumeContiguousRuns = cms.bool(True),
                               histograms = cms.bool(False)
                               )

