import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
# process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
    'file:/tmp/adamwo/Hermine/src/JpsiMM_cfi_py_GEN_SIM.root'
    )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

process.demo = cms.EDFilter('MuonSimHitCountFilter',
                            simTracks = cms.InputTag("g4SimHits"),
                            simHits = cms.VInputTag(cms.InputTag("g4SimHits","MuonCSCHits"),
                                                    cms.InputTag("g4SimHits","MuonDTHits")),
                            minHits = cms.vint32(2,2),
                            minSumHits = cms.int32(2),
                            particleTypes = cms.vint32(-13,13),
                            processTypes = cms.vint32()
                            )


process.p = cms.Path(process.demo)
