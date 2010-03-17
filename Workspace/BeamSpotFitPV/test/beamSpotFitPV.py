import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO.limit = 999999

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/tmp/adamwo/skimPV-Feb9.root'
#    #        '/store/data/BeamCommissioning09/ZeroBias/RECO/Feb9ReReco_v2/0026/E6723795-0B16-DF11-8555-00163E0101D9.root'
#    '/store/data/BeamCommissioning09/ZeroBias/RECO/Feb9ReReco_v2/0027/20A22C24-8716-DF11-9FFA-00163E010269.root',
#    '/store/data/BeamCommissioning09/ZeroBias/RECO/Feb9ReReco_v2/0026/ACA52F85-0716-DF11-8D6E-00163E010524.root',
#    '/store/data/BeamCommissioning09/ZeroBias/RECO/Feb9ReReco_v2/0026/7CA4E4AF-1216-DF11-9757-00163E010547.root'
    ) #,   lumisToProcess = cms.untracked.VLuminosityBlockRange('124120:1-124120:max')
                   
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("histo.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )


process.load("Workspace.BeamSpotFitPV.beamSpotFitPV_cfi")
#process.beamSpotFitPV.vertexCacheSize = 500
##process.beamSpotFitPV.minNrVerticesForFit = 100
#process.beamSpotFitPV.vertexTag = cms.InputTag("hltPixelVertices")
#process.beamSpotFitPV.histograms = True

process.p = cms.Path(process.beamSpotFitPV)
