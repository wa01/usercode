import FWCore.ParameterSet.Config as cms

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
  process.hltL1sJPsi = cms.EDFilter( "HLTLevel1GTSeed",
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMuOpen OR L1_SingleMu0 OR L1_SingleMu3" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" )
  )
  process.hltPreJPsi = cms.EDFilter( "HLTPrescaler" )
  process.hltJPsiMuL1Filter = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "hltL1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sJPsi" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    SelectQualities = cms.vint32(  )
  )
  process.hltJPsiMuL2Filter = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "hltOfflineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltJPsiMuL1Filter" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 )
  )
  process.hltJPsiMuL3Filter = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "hltOfflineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltJPsiMuL2Filter" ),
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
  process.hltJPsiPixelTrackFilter = cms.EDFilter("HLTOniaTrackFilter",
                                       candTag = cms.InputTag("hltPixelTrackCands"),
                                       MinPt = cms.double(0.),
                                       MinP = cms.double(3.),
                                       MaxEta = cms.double(2.5),
                                       saveTag = cms.bool(True)
                                       )
  process.hltJPsiPixelMassFilter = cms.EDFilter("HLTOniaMuonTrackMassFilter",
                                      muonCandidates = cms.InputTag("hltJPsiMuL3Filter"),
                                      trackCandidates = cms.InputTag("hltJPsiPixelTrackFilter"),
#                                      MinMass = cms.double(2.6),
#                                      MaxMass = cms.double(3.6),
                                      MinMass = cms.double(1.6),
                                      MaxMass = cms.double(4.6),
                                      checkCharge = cms.bool(False),
                                      saveTag = cms.bool(True)
                                      )
  process.hltJPsiPixelTrackSequence = cms.Sequence( process.HLTDoLocalPixelSequence + process.hltPixelTracks + process.hltPixelTrackCands )
  process.hltPixelCandsToTracks = cms.EDProducer("RecoChargedCandidatesToTracks",
                                       candTag = cms.InputTag("hltJPsiPixelMassFilter")
                                       )
  process.hltJPsiSeeds = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
                              InputCollection = cms.InputTag("hltPixelCandsToTracks"),
                              TTRHBuilder = cms.string("WithTrackAngle"),
                              useProtoTrackKinematics = cms.bool(False)
                              )
  process.hltJPsiCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltJPsiSeeds" ),
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
  process.hltJPsiCtfTracks = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltJPsiCtfTracks" ),
    Fitter = cms.string( "FittingSmootherRK" ),
    Propagator = cms.string( "RungeKuttaTrackerPropagator" ),
    src = cms.InputTag( "hltJPsiCkfTrackCandidates" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
  )
  process.hltJPsiCtfTrackCands = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltJPsiCtfTracks" ),
    particleType = cms.string( "mu-" )
  )
  process.hltJPsiTrackSequence = cms.Sequence( process.HLTDoLocalStripSequence + process.hltPixelCandsToTracks + process.hltJPsiSeeds + process.hltJPsiCkfTrackCandidates
                                     + process.hltJPsiCtfTracks + process.hltJPsiCtfTrackCands )

  process.hltJPsiCtfTrackFilter = cms.EDFilter("HLTOniaTrackFilter",
                                     candTag = cms.InputTag("hltJPsiCtfTrackCands"),
                                     MinPt = cms.double(0.),
                                     MinP = cms.double(3.),
                                     MaxEta = cms.double(2.5),
                                     saveTag = cms.bool(True)
                                     )
  process.hltJPsiCtfMassFilter = cms.EDFilter("HLTOniaMuonTrackMassFilter",
                                    muonCandidates = cms.InputTag("hltJPsiMuL3Filter"),
                                    trackCandidates = cms.InputTag("hltJPsiCtfTrackFilter"),
#                                    MinMass = cms.double(2.9),
#                                    MaxMass = cms.double(3.3),
                                    MinMass = cms.double(2.8),
                                    MaxMass = cms.double(3.4),
                                    checkCharge = cms.bool(True),
                                    saveTag = cms.bool(True)
                                    )

  process.HLT_JPsi_8E29 = cms.Path( process.HLTBeginSequence + process.hltL1sJPsi + process.hltPreJPsi
                          + process.hltJPsiMuL1Filter + process.HLTL2muonrecoSequence + process.hltJPsiMuL2Filter
                          + process.HLTL3muonrecoSequence + process.hltJPsiMuL3Filter
                          + process.hltJPsiPixelTrackSequence + process.hltJPsiPixelTrackFilter + process.hltJPsiPixelMassFilter
                          + process.hltJPsiTrackSequence  + process.hltJPsiCtfTrackFilter + process.hltJPsiCtfMassFilter
                          + process.HLTEndSequence )

  process.schedule.insert( process.schedule.index(process.HLTriggerFinalPath), process.HLT_JPsi_8E29 )

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

