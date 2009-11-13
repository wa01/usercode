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
  process.hltL1sOnia = cms.EDFilter( "HLTLevel1GTSeed",
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMuOpen OR L1_SingleMu0 OR L1_SingleMu3" ),
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
  process.HLTOniaPixelTrackSequence = cms.Sequence( process.HLTDoLocalPixelSequence
                                                    + process.hltPixelTracks + process.hltPixelTrackCands )
  process.hltOniaPixelTrackSelector = cms.EDProducer("QuarkoniaTrackSelector",
    muonCandidates = cms.InputTag("hltL3MuonCandidates"),
    tracks = cms.InputTag("hltPixelTracks"),
    MinMasses = cms.vdouble( 2.6 ),
    MaxMasses = cms.vdouble( 3.6 ),
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
    BeamSpotTag = cms.InputTag("hltOfflineBeamSpot"),
    CandTag = cms.InputTag("hltL3MuonCandidates"),
    TrackTag = cms.InputTag("hltOniaPixelTrackCands"),
    PreviousCandTag = cms.InputTag("hltOniaMuL3Filter"),
    MinMasses = cms.vdouble(2.6),
    MaxMasses = cms.vdouble(3.6),
    checkCharge = cms.bool(False),
    MinTrackPt = cms.double(0.),
    MinTrackP = cms.double(3.),
    MaxTrackEta = cms.double(999.),
    MaxTrackDxy = cms.double(999.),
    MaxTrackDz = cms.double(999.),
    MinTrackHits = cms.int32(3),
    MaxTrackNormChi2 = cms.double(999999999.),
    MaxDzMuonTrack = cms.double(999.),
    SaveTag = cms.untracked.bool(True)
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
  process.HLTOniaTrackSequence = cms.Sequence( process.HLTDoLocalStripSequence + process.hltOniaSeeds
                                               + process.hltOniaCkfTrackCandidates
                                               + process.hltOniaCtfTracks + process.hltOniaCtfTrackCands )
  process.hltOniaCtfMassFilter = cms.EDFilter("HLTMuonTrackMassFilter",
    BeamSpotTag = cms.InputTag("hltOfflineBeamSpot"),
    CandTag = cms.InputTag("hltL3MuonCandidates"),
    TrackTag = cms.InputTag("hltOniaCtfTrackCands"),
    PreviousCandTag = cms.InputTag("hltOniaPixelMassFilter"),
    MinMasses = cms.vdouble( 2.9 ),
    MaxMasses = cms.vdouble( 3.3 ),
    checkCharge = cms.bool(True),
    MinTrackPt = cms.double(0.),
    MinTrackP = cms.double(3.),
    MaxTrackEta = cms.double(999.),
    MaxTrackDxy = cms.double(999.),
    MaxTrackDz = cms.double(999.),
    MinTrackHits = cms.int32(5),
    MaxTrackNormChi2 = cms.double(999999999.),
    MaxDzMuonTrack = cms.double(0.5),
    SaveTag = cms.untracked.bool(True)
  )


  process.HLT_Onia_1E31 = cms.Path( process.HLTBeginSequence + process.hltL1sOnia + process.hltPreOnia
                          + process.hltOniaMuL1Filter + process.HLTL2muonrecoSequence + process.hltOniaMuL2Filter
                          + process.HLTL3muonrecoSequence + process.hltOniaMuL3Filter
                          + process.HLTOniaPixelTrackSequence + process.hltPixelTrackCands
                          + process.hltOniaPixelTrackSelector + process.hltOniaPixelTrackCands
                          + process.hltOniaPixelMassFilter
                          + process.HLTOniaTrackSequence + process.hltOniaCtfMassFilter
                          + process.HLTEndSequence )

  process.schedule.insert( process.schedule.index(process.HLTriggerFinalPath), process.HLT_Onia_1E31 )

  return process
