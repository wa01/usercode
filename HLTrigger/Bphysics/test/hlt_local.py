import FWCore.ParameterSet.Config as cms

process = cms.Process("MYHLT")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#    'file:/tmp/adamwo/JPsiMuMu/GEN-SIM-RAW/MC_31X_V3_7TeV-v1/16D3A6FF-94A5-DE11-83AB-002481DE4C86.root'
#    'rfio:/castor/cern.ch/user/a/apana/L1Skims/3_1_2/1e31/MinBias_Summer09_1e31_0.root',
#    'rfio:/castor/cern.ch/user/a/apana/L1Skims/3_1_2/1e31/MinBias_Summer09_1e31_1.root'
    'file:/tmp/adamwo/JPsiMuMu/GEN-SIM-RAW/MC_31X_V3-v1/0041/622EDA2D-7E8B-DE11-A882-00304877A648.root'
    ),
#                            lumisToProcess = cms.untracked.VLuminosityBlockRange("1:100890-1:100890")
#                            eventsToProcess = cms.untracked.VEventRange(
#"1:12262561",
#"1:12262770",
#"1:12262804",
#"1:12263592",
#"1:12264699",
#"1:12265155",
#"1:12265835",
#"1:12265849",
#"1:12266072",
#"1:12266214",
#"1:12266242",   
#    )
                            )

process.maxEvents = cms.untracked.PSet(   input = cms.untracked.int32( -1 )   )

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('HLTrigger/Configuration/HLT_8E29_cff')
#process.load('HLTrigger/Configuration/HLT_1E31_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

#process.GlobalTag.globaltag = 'MC_31X_V3::All'
process.GlobalTag.globaltag = 'STARTUP31X_V2::All'

process.MessageLogger.categories.append('HLTTrackFilter')
process.MessageLogger.categories.append('HLTJPsiFilter')
process.MessageLogger.categories.append('HLTMuonL3PreFilter')
process.MessageLogger.debugModules = ['hltJPsiPixelTrackFilter','hltJPsiPixelMassFilter',
                                      'hltJPsiCtfTrackFilter','hltJPsiCtfMassFilter',
                                      'hltJPsiMuL3Filter' ]
process.MessageLogger.cerr.threshold = "DEBUG"
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

#process.Tracer = cms.Service("Tracer")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histo.root")
                                   )

#process.load("L1Trigger.GlobalTriggerAnalyzer.l1GtTrigReport_cfi")
#process.l1GtTrigReport.PrintOutput = 0


process.load("HLTrigger/HLTcore/triggerSummaryAnalyzerRAW_cfi")
process.triggerSummaryAnalyzerRAW.inputTag = cms.InputTag("hltTriggerSummaryRAW","","MYHLT")

process.load("HLTrigger/HLTcore/triggerSummaryAnalyzerAOD_cfi")
process.triggerSummaryAnalyzerAOD.inputTag = cms.InputTag("hltTriggerSummaryAOD","","MYHLT")


process.contentAnalyzer = cms.EDAnalyzer("EventContentAnalyzer")

process.l3analyzer = cms.EDAnalyzer("JPsiHltAnalyzer",
                                    muonCandidates = cms.InputTag("hltL3MuonCandidates"),
                                    trackCandidates = cms.InputTag("hltPixelTrackCands"),
                                    muonFilteredCandidates = cms.InputTag("hltJPsiMuL3Filter"),
                                    trackFilteredCandidates = cms.InputTag("hltJPsiPixelTrackFilter"),
                                    jpsiCandidates = cms.InputTag("hltJPsiPixelMassFilter")
                                    )
process.l35analyzer = cms.EDAnalyzer("JPsiHltAnalyzer",
                                     muonCandidates = cms.InputTag("hltL3MuonCandidates"),
                                     trackCandidates = cms.InputTag("hltJPsiCtfTrackCands"),
                                     muonFilteredCandidates = cms.InputTag("hltJPsiMuL3Filter"),
                                     trackFilteredCandidates = cms.InputTag("hltJPsiCtfTrackFilter"),
                                     jpsiCandidates = cms.InputTag("hltJPsiCtfMassFilter")
                                     )
process.analyzerPath = cms.Path(
#    process.triggerSummaryAnalyzerRAW+
    process.triggerSummaryAnalyzerAOD #+
#    process.contentAnalyzer #+
#                                process.l3analyzer+process.l35analyzer
#                                +process.l1GtTrigReport
    )
process.schedule = process.HLTSchedule
#process.schedule.extend([process.jpsiPath,process.analyzerPath])
#process.schedule.extend([ #process.analyzerPath,
#    process.jpsiPath,process.jpsiAnalyzerPath])
process.schedule.extend([process.analyzerPath])

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('JPsiMuMu_HLT_JPsi_8E29.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-HLTDEBUG'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('HLT_JPsi_8E29'
                                   #,'HLT_Mu*','HLT_DoubleMu*'
                                   )
    )
)

#process.out_step = cms.EndPath(process.output)
#process.schedule.extend([process.out_step])
