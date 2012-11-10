import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.GlobalTag.globaltag = "START53_V7A::All"
process.GlobalTag.globaltag = "FT_R_44_V11::All"


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/home/adamwo/scratch/C264F70C-1915-E111-874B-00A0D1EEF6B8.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histo.root")
                                   )
# filter on good vertex
process.pvFilter = cms.EDFilter("PrimaryVertexObjectFilter",
                                src = cms.InputTag('offlinePrimaryVertices'),
                                filterParams = cms.PSet(
    minNdof = cms.double(4),
    maxZ = cms.double(24),
    maxRho = cms.double(2)
    )
                                )

process.demo = cms.EDAnalyzer('DiPiAnalyzer',
                              pvSrc = cms.InputTag("pvFilter"),
                              trackSrc = cms.InputTag("generalTracks"),
                              minRawMass = cms.double(1.),
                              maxRawMass = cms.double(999.),
                              maxVtxChi2 = cms.double(10.),
                              minVtxMass = cms.double(2.),
                              maxVtxMass = cms.double(20.),
                              sigTipLimit = cms.double(3.)
                              )


process.p = cms.Path(process.pvFilter*process.demo)

#process.source.fileNames.extend(
#["/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_10_1_sy8.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_16_3_cYl.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_17_3_nMl.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_1_1_mtb.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_21_3_hlL.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_23_3_yZv.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_24_3_FTD.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_25_3_sdJ.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_26_3_SI5.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_28_3_lKH.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_29_3_Nnx.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_2_1_QmM.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_31_3_zcI.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_33_3_ePl.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_34_3_HuR.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_35_3_pyd.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_36_3_Wn3.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_37_3_YTu.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_38_3_NR3.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_39_3_wz4.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_40_3_2qx.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_41_3_Qdf.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_42_3_bq3.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_43_3_8gD.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_44_3_Xpp.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_45_3_vf3.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_46_3_zhn.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_47_3_fRZ.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_4_1_eAO.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_50_3_qdV.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_51_3_utD.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_53_3_hAp.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_55_3_Cl9.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_56_3_tD2.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_57_3_FR8.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_59_3_ryp.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_5_1_eUM.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_61_3_9WE.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_62_3_9RP.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_63_3_BNO.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_64_3_w2D.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_65_3_7yR.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_66_3_QBg.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_67_3_Vu6.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_68_3_XNd.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_69_3_6cC.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_6_1_UUa.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_70_3_qwR.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_71_3_0lF.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_74_3_Eoq.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_75_3_ISK.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_77_3_OL7.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_79_3_1iu.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_7_1_DyU.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_80_3_xuQ.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_8_1_ari.root",
#"/store/user/adamwo/Conversions/MinBias_TuneZ2star_8TeV-pythia6/chib_ConvFilter_SelEvts_9_1_hCo.root"]
#)
