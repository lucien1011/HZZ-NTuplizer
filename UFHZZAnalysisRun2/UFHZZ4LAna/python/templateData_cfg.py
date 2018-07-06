import FWCore.ParameterSet.Config as cms

process = cms.Process("UFHZZ4LAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag='74X_dataRun2_reMiniAOD_v0'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring(DUMMYFILELIST)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )

import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'Run2015_30Oct_JSON.txt').getVLuminosityBlockRange()
process.source.lumisToProcess = LumiList.LumiList(filename = 'Run2015_Nov13_Silver.txt').getVLuminosityBlockRange()

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DUMMYFILENAME.root")
)

# clean muons by segments 
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
				     src = cms.InputTag("slimmedMuons"),
				     preselection = cms.string("track.isNonnull"),
				     passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
				     fractionOfSharedSegments = cms.double(0.499),
				     )

# Electron Calibrations
process.calibratedElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    electrons = cms.InputTag("slimmedElectrons","","PAT"),
    grbForestName = cms.string("gedelectron_p4combination_25ns"),
    isMC = cms.bool(False),
    isSynchronization = cms.bool(False)
)

# Electron MVA ID producer
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
# add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag("calibratedElectrons")

process.mvaSpring15NonTrig25nsV1 = cms.EDProducer("SlimmedElectronMvaIDProducer",
                                     mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                     electronsCollection = cms.InputTag("calibratedElectrons"),
                                     Trig = cms.bool(False),
                                     )
     
# FSR Photons
process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

# Jet Energy Corrections
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.patJetCorrFactorsUpdated.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute',
              'L2L3Residual'
              ],
    payload = 'AK4PFchs' ) 

process.slimmedJetsJEC = process.patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )

# UF Ntuplizer
process.Ana = cms.EDAnalyzer('UFHZZ4LAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("mvaSpring15NonTrig25nsV1","NonTrig"),
                              muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
#                              jetSrc       = cms.untracked.InputTag("slimmedJetsJEC"),
                              metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"), #or selectedVertices 
                              isMC         = cms.untracked.bool(False),
                              isSignal     = cms.untracked.bool(False),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(1.0),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(False),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              rhoSrcSUS    = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              pileupSrc     = cms.untracked.InputTag("slimmedAddPileupInfo"),
                              reweightForPU = cms.untracked.bool(False),
                              triggerSrc = cms.untracked.InputTag("TriggerResults","","HLT"),
                              triggerObjects = cms.InputTag("selectedPatTrigger"),
                              triggerList = cms.untracked.vstring(
                                    'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                    'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
                                    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
                                    'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
                                    'HLT_TripleMu_12_10_5_v',
                                    'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
                                    'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                    'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
                                    'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v',
                                    'HLT_Ele27_WPLoose_Gsf_v', 
                              ),
                              verbose = cms.untracked.bool(False)              
#                              verbose = cms.untracked.bool(True)              
                             )


process.p = cms.Path(process.fsrPhotonSequence*
                     process.boostedMuons*
                     process.calibratedElectrons*
                     process.electronMVAValueMapProducer*
                     process.mvaSpring15NonTrig25nsV1*
#                     process.jetCorrFactors*
#                     process.slimmedJetsJEC*
                     process.Ana
                     )
