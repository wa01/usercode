import FWCore.ParameterSet.Config as cms

process = cms.Process("MYHLT")


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#    'file:/tmp/adamwo/JPsiMuMu/GEN-SIM-RAW/MC_31X_V3_7TeV-v1/16D3A6FF-94A5-DE11-83AB-002481DE4C86.root'
#    'rfio:/castor/cern.ch/user/a/apana/L1Skims/3_1_2/1e31/MinBias_Summer09_1e31_0.root',
#    'rfio:/castor/cern.ch/user/a/apana/L1Skims/3_1_2/1e31/MinBias_Summer09_1e31_1.root'
#    'file:/tmp/adamwo/622EDA2D-7E8B-DE11-A882-00304877A648.root'
#    '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v1/0000/DAD8907F-B1C8-DE11-ACA7-001D09F2546F.root'
 '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/B65667F8-9EC8-DE11-ABA3-0018F3D095EA.root',
 '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/A4805FB4-9BC8-DE11-9562-003048678F9C.root',
 '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/8AA3E877-96C8-DE11-BF14-0018F3D09688.root',
 '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/66934AD9-15C9-DE11-B455-001A92971B36.root',
 '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/5ACC083B-A6C8-DE11-8AF2-001731AF68C9.root',
 '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/44357DC6-99C8-DE11-985F-001A92810AA4.root',
 '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/405B5744-A2C8-DE11-A0E2-0018F3D096A2.root',
 '/store/relval/CMSSW_3_3_2/RelValJpsiMM/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V8-v2/0000/28A228FA-A0C8-DE11-9531-001A92811738.root'
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

#process.GlobalTag.globaltag = 'MC_31X_V9::All'
process.GlobalTag.globaltag = 'STARTUP31X_V8::All'

#process.MessageLogger.categories.append('HLTTrackFilter')
#process.MessageLogger.categories.append('HLTJPsiFilter')
#process.MessageLogger.categories.append('HLTMuonL3PreFilter')
process.MessageLogger.categories.append('QuarkoniaTrackSelector')
process.MessageLogger.categories.append('HLTMuonTrackMassFilter')
process.MessageLogger.debugModules = [
  # 'hltJPsiPixelTrackFilter','hltJPsiPixelMassFilter',
   #                                   'hltJPsiCtfTrackFilter','hltJPsiCtfMassFilter',
   #                                   'hltJPsiMuL3Filter',
#  'muonPlusPixelTrack'
  'hltOniaPixelTrackSelector',
  'hltOniaPixelMassFilter', 'hltOniaCtfMassFilter'
  ]
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

process.l1analyzer = cms.EDAnalyzer("OniaHLTAnalyzer",
                                    muonCandidates = cms.InputTag("hltL1extraParticles"),
                                    trackCandidates = cms.InputTag("none"),
                                    muonFilteredCandidates = cms.InputTag("hltOniaMuL1Filter"),
                                    trackFilteredCandidates = cms.InputTag("none"),
                                    jpsiCandidates = cms.InputTag("none")
                                    )
process.l2analyzer = cms.EDAnalyzer("OniaHLTAnalyzer",
                                    muonCandidates = cms.InputTag("hltL2MuonCandidates"),
                                    trackCandidates = cms.InputTag("none"),
                                    muonFilteredCandidates = cms.InputTag("hltOniaMuL2Filter"),
                                    trackFilteredCandidates = cms.InputTag("none"),
                                    jpsiCandidates = cms.InputTag("none")
                                    )
process.l3analyzer = cms.EDAnalyzer("OniaHLTAnalyzer",
                                    muonCandidates = cms.InputTag("hltL3MuonCandidates"),
                                    trackCandidates = cms.InputTag("hltPixelTrackCands"),
                                    muonFilteredCandidates = cms.InputTag("hltOniaMuL3Filter"),
                                    trackFilteredCandidates = cms.InputTag("hltOniaPixelTrackCands"),
                                    jpsiCandidates = cms.InputTag("hltOniaPixelMassFilter")
                                    )
process.l35analyzer = cms.EDAnalyzer("OniaHLTAnalyzer",
                                     muonCandidates = cms.InputTag("hltL3MuonCandidates"),
                                     trackCandidates = cms.InputTag("hltOniaCtfTrackCands"),
                                     muonFilteredCandidates = cms.InputTag("hltOniaMuL3Filter"),
                                     trackFilteredCandidates = cms.InputTag("hltOniaCtfTrackCands"),
                                     jpsiCandidates = cms.InputTag("hltOniaCtfMassFilter")
                                     )
process.analyzerPath = cms.Path(
#    process.triggerSummaryAnalyzerRAW+
#    process.triggerSummaryAnalyzerAOD #+
#    process.muonPlusPixelTrack+
#    process.contentAnalyzer+
  process.l2analyzer+
                                process.l3analyzer+process.l35analyzer
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


def customise(process):

  process.oniaTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
    ComponentName = cms.string( "oniaTrajectoryBuilder" ),
    updator = cms.string( "KFUpdator" ),
    propagatorAlong = cms.string( "PropagatorWithMaterial" ),
    propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
    estimator = cms.string( "Chi2" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    MeasurementTrackerName = cms.string( "" ),
    trajectoryFilterName = cms.string( "oniaTrajectoryFilter" ),
    maxCand = cms.int32( 1 ),
    lostHitPenalty = cms.double( 30.0 ),
    intermediateCleaning = cms.bool( True ),
    alwaysUseInvalidHits = cms.bool( False ),
    appendToDataLabel = cms.string( "" )
  )
  process.oniaTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
    ComponentName = cms.string( "oniaTrajectoryFilter" ),
    appendToDataLabel = cms.string( "" ),
    filterPset = cms.PSet(
      chargeSignificance = cms.double( -1.0 ),
      minPt = cms.double( 1.0 ),
      minHitsMinPt = cms.int32( 3 ),
      ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
      maxLostHits = cms.int32( 1 ),
      maxNumberOfHits = cms.int32( 8 ),
      maxConsecLostHits = cms.int32( 1 ),
      nSigmaMinPt = cms.double( 5.0 ),
      minimumNumberOfHits = cms.int32( 5 )
    )
  )
  process.hltL1sOnia = cms.EDFilter( "HLTLevel1GTSeed",
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
#    L1SeedsLogicalExpression = cms.string( "L1_SingleMuOpen OR L1_SingleMu0 OR L1_SingleMu3" ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu3" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" )
  )
  process.hltPreOnia = cms.EDFilter( "HLTPrescaler" )
  process.hltOniaMuL1Filter = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "hltL1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sOnia" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    SelectQualities = cms.vint32(  )
  )
  process.hltOniaMuL2Filter = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "hltOfflineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltOniaMuL1Filter" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 )
  )
  process.hltOniaMuL3Filter = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "hltOfflineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltOniaMuL2Filter" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    SaveTag = cms.untracked.bool( True )
  )
  process.hltPixelTrackCands = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltPixelTracks" ),
    particleType = cms.string( "mu-" )
  )
  process.hltOniaPixelTrackSequence = cms.Sequence( process.HLTDoLocalPixelSequence
                                                    + process.hltPixelTracks + process.hltPixelTrackCands )
  process.hltOniaPixelTrackSelector = cms.EDProducer("QuarkoniaTrackSelector",
    muonCandidates = cms.InputTag("hltL3MuonCandidates"),
    tracks = cms.InputTag("hltPixelTracks"),
    MinMasses = cms.vdouble( 1.6 ),
    MaxMasses = cms.vdouble( 4.6 ),
    checkCharge = cms.bool(False),
    MinTrackPt = cms.double(0.),
    MinTrackP = cms.double(3.),
    MaxTrackEta = cms.double(999.)
  )
  process.hltOniaPixelTrackCands = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltOniaPixelTrackSelector" ),
    particleType = cms.string( "mu-" )
  )
  process.hltOniaPixelMassFilter = cms.EDFilter("HLTMuonTrackMassFilter",
    beamspot = cms.InputTag("hltOfflineBeamSpot"),
    muonCandidates = cms.InputTag("hltL3MuonCandidates"),
    trackCandidates = cms.InputTag("hltOniaPixelTrackCands"),
    previousCandidates = cms.InputTag("hltOniaMuL3Filter"),
    MinMasses = cms.vdouble(1.6),
    MaxMasses = cms.vdouble(4.6),
    checkCharge = cms.bool(False),
    MinTrackPt = cms.double(0.),
    MinTrackP = cms.double(3.),
    MaxTrackEta = cms.double(999.),
    MaxTrackDxy = cms.double(999.),
    MaxTrackDz = cms.double(999.),
    MinTrackHits = cms.int32(3),
    MaxTrackNormChi2 = cms.double(999999999.),
    MaxDzMuonTrack = cms.double(999.),
    saveTag = cms.bool(True)
  )
  process.hltOniaSeeds = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltOniaPixelTrackSelector"),
    TTRHBuilder = cms.string("WithTrackAngle"),
    useProtoTrackKinematics = cms.bool(False)
    )
  process.hltOniaCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltOniaSeeds" ),
    TrajectoryBuilder = cms.string( "oniaTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet(
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False )
  )
  process.hltOniaCtfTracks = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltOniaCtfTracks" ),
    Fitter = cms.string( "FittingSmootherRK" ),
    Propagator = cms.string( "RungeKuttaTrackerPropagator" ),
    src = cms.InputTag( "hltOniaCkfTrackCandidates" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
  )
  process.hltOniaCtfTrackCands = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltOniaCtfTracks" ),
    particleType = cms.string( "mu-" )
  )
  process.hltOniaTrackSequence = cms.Sequence( process.HLTDoLocalStripSequence + process.hltOniaSeeds
                                               + process.hltOniaCkfTrackCandidates
                                               + process.hltOniaCtfTracks + process.hltOniaCtfTrackCands )
  process.hltOniaCtfMassFilter = cms.EDFilter("HLTMuonTrackMassFilter",
    beamspot = cms.InputTag("hltOfflineBeamSpot"),
    muonCandidates = cms.InputTag("hltL3MuonCandidates"),
    trackCandidates = cms.InputTag("hltOniaCtfTrackCands"),
    previousCandidates = cms.InputTag("hltOniaPixelMassFilter"),
    MinMasses = cms.vdouble( 2.8 ),
    MaxMasses = cms.vdouble( 3.4 ),
    checkCharge = cms.bool(True),
    MinTrackPt = cms.double(0.),
    MinTrackP = cms.double(3.),
    MaxTrackEta = cms.double(999.),
    MaxTrackDxy = cms.double(999.),
    MaxTrackDz = cms.double(999.),
    MinTrackHits = cms.int32(5),
    MaxTrackNormChi2 = cms.double(999999999.),
    MaxDzMuonTrack = cms.double(0.5),
    saveTag = cms.bool(True)
  )


  process.HLT_Onia_8E29 = cms.Path( process.HLTBeginSequence + process.hltL1sOnia + process.hltPreOnia
                          + process.hltOniaMuL1Filter + process.HLTL2muonrecoSequence + process.hltOniaMuL2Filter
                          + process.HLTL3muonrecoSequence + process.hltOniaMuL3Filter
                          + process.hltOniaPixelTrackSequence + process.hltPixelTrackCands
                          + process.hltOniaPixelTrackSelector + process.hltOniaPixelTrackCands
                          + process.hltOniaPixelMassFilter
                          + process.hltOniaTrackSequence + process.hltOniaCtfMassFilter
                          + process.HLTEndSequence )

  process.schedule.insert( process.schedule.index(process.HLTriggerFinalPath), process.HLT_Onia_8E29 )

#  process.load("HLTriggerOffline.HeavyFlavor.heavyFlavorValidationSequence_cff")
#  process.load("HLTrigger/HLTanalyzers/hlTrigReport_cfi")
#  process.heavyFlavorValidation_step = cms.Path(process.heavyFlavorValidationSequence+process.hlTrigReport)
#  process.schedule.insert( process.schedule.index(process.endjob_step), process.heavyFlavorValidation_step )
#  process.output.outputCommands = cms.untracked.vstring(
#      'drop *',
#      'keep *_MEtoEDMConverter_*_*'
#  )
#  process.output.fileName = cms.untracked.string('/tmp/heavyFlavorValidation.root')

  return process


process = customise(process)

