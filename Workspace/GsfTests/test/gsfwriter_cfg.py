import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.GlobalTag.globaltag = 'START44_V6::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories = cms.untracked.vstring('GsfWriter')
process.MessageLogger.cerr.threshold = "INFO"
process.MessageLogger.cerr.INFO.limit = 0
process.MessageLogger.cerr.GsfWriter = cms.untracked.PSet(
    limit = cms.untracked.int32(10000000)
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    '/store/relval/CMSSW_4_4_1/RelValSingleElectronPt35/GEN-SIM-RECO/START44_V6-v1/0057/9212A461-B7F7-E011-A19E-002354EF3BE4.root'
   '/store/relval/CMSSW_4_4_1/RelValSingleElectronPt10/GEN-SIM-RECO/START44_V6-v1/0058/6A8C17E5-B9F7-E011-8487-0026189437F2.root'
    )
)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("histo.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )


#process.anyGsfElectronPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
#    PropagationDirection = cms.string('anyDirection'),
#    ComponentName = cms.string('anyGsfElectronPropagator'),
#    Mass = cms.double(0.000511),
#    ptMin = cms.double(-1.0),
#    MaxDPhi = cms.double(1.6),
#    useRungeKutta = cms.bool(False)
#)

process.demo = cms.EDAnalyzer('GsfWriter',
                              gsfTracks = cms.InputTag("electronGsfTracks"),
                              simTracks = cms.InputTag("g4SimHits"),
                              simVertices = cms.InputTag("g4SimHits")
)


process.p = cms.Path(process.demo)
