import FWCore.ParameterSet.Config as cms

tracksFromQuadSeeds = cms.EDProducer("TrackProducerFromSeed",
    src = cms.InputTag("photonConvTrajSeedFromQuadruplets","conv2SeedCandidates"),
#    clusterRemovalInfo = cms.InputTag(""),
    beamSpot = cms.InputTag("offlineBeamSpot"),
#    Fitter = cms.string('KFFittingSmootherWithOutliersRejectionAndRK'),
    useHitsSplitting = cms.bool(False),
    alias = cms.untracked.string('mySeedTracks'),
    TrajectoryInEvent = cms.bool(False),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    AlgorithmName = cms.string('undefAlgorithm'),
    Propagator = cms.string('RungeKuttaTrackerPropagator'),

    ### These are paremeters related to the filling of the Secondary hit-patterns
    #set to "", the secondary hit pattern will not be filled (backward compatible with DetLayer=0)
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    MeasurementTracker = cms.string('')
)



