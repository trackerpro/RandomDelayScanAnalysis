import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser                                                                                                                                        
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('delayStep',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,
                  'integer numer to identify the delay xml for each partition');

options.register ('ouputFileName',"trackerDPG.root",VarParsing.multiplicity.singleton,VarParsing.varType.string,
                  'name of the outtput root file');

options.register ('jsonFile',"",VarParsing.multiplicity.singleton,VarParsing.varType.string,
                  'json file to apply in case one wants to ....');

options.register ('isRawDAQFile',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                  'when running as input a raw DAQ file instead of a standard FEVT/AOD');

options.register ('isRawEDMFile',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
                  'when running as input a raw EDM file instead of a standard FEVT/AOD')

options.register ('globalTag',"run3_data_express",VarParsing.multiplicity.singleton,VarParsing.varType.string,
                  'Name of the global tag used for recostruct data');

options.register ('nThreads',1,VarParsing.multiplicity.singleton,VarParsing.varType.int,
                  'Number of threads in the processing')

options.register ('inputDirectory',"",VarParsing.multiplicity.singleton,VarParsing.varType.string,
                  'directory where the xml files with the delays are stored')

options.register ('triggerList','',VarParsing.multiplicity.list,
                  VarParsing.varType.string,"List of HLT trigger paths to require")

options.register ('eventRange',"",VarParsing.multiplicity.singleton,VarParsing.varType.string,
                  '');

options.parseArguments()

### start cmssw job
import FWCore.PythonUtilities.LumiList as LumiList
from Configuration.StandardSequences.Eras import eras

process = cms.Process("clusterAnalysis",eras.Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

if not options.isRawDAQFile:
    process.source = cms.Source ("PoolSource", 
                                 fileNames = cms.untracked.vstring(options.inputFiles))
else:
    process.source = cms.Source("NewEventStreamFileReader",
                                fileNames = cms.untracked.vstring(options.inputFiles))


if options.eventRange:
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventRange);

# Conditions (Global Tag is used here):
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:'+options.globalTag,'');

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(options.nThreads),
    numberOfStreams = cms.untracked.uint32(options.nThreads))


#Geometry and magnetic field to be loaded
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Geometry.CommonTopologies.globalTrackingGeometryDB_cfi")
process.load("TrackingTools.RecoGeometry.RecoGeometries_cff")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")
process.load('RecoTracker.DeDx.dedxEstimators_cff')

### to select a specific set of triggers --> if one runs on the streamer files at T0 is no needed
process.hltfiter = cms.EDFilter("HLTHighLevel",
  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
  HLTPaths = cms.vstring(options.triggerList),
  throw = cms.bool(True),
  andOr = cms.bool(True) ## logical OR between trigger bits                               
)


### output file definition
process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.ouputFileName)
)


### when one runs on a AOD or FEVT file and wants to re-fit tracks from already existing ones
if not options.isRawDAQFile and not options.isRawEDMFile: 

    process.load('RecoTracker.DeDx.dedxEstimators_cff')
    process.load("RecoTracker.TrackProducer.TrackRefitter_cfi")
    process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

    process.TrackRefitter.NavigationSchool = ''
    process.TrackRefitter.Fitter = 'FlexibleKFFittingSmoother'
    process.generalTracksFromRefit = process.TrackRefitter.clone(
        src = cms.InputTag("generalTracks")
    )
    
    process.doAlldEdXEstimators    += process.dedxMedian
    process.dedxTruncated40.tracks = "generalTracksFromRefit"
    process.dedxHarmonic2.tracks   = "generalTracksFromRefit"
    process.dedxMedian.tracks      = "generalTracksFromRefit"
    process.dedxHitInfo.tracks     = "generalTracksFromRefit"
    process.dedxPixelHarmonic2.tracks  = "generalTracksFromRefit"
    process.dedxPixelAndStripHarmonic2T085.tracks  = "generalTracksFromRefit"

### run the tracking step on top of raw DAQ files
elif options.isRawDAQFile or options.isRawEDMFile:

    process.load('Configuration.StandardSequences.RawToDigi_cff')    
    process.load('Configuration.StandardSequences.L1Reco_cff')
    process.load('Configuration.StandardSequences.Reconstruction_cff')

    ## save trajectories everywhere since they are used by the TrackerDpgAnalyzer
    process.generalTracks.copyTrajectories = cms.untracked.bool(True);
    process.mergedDuplicateTracks.TrajectoryInEvent = cms.bool(True);
    process.preDuplicateMergingGeneralTracks.copyTrajectories = cms.untracked.bool(True);
    process.earlyGeneralTracks.copyTrajectories = cms.untracked.bool(True);
    process.initialStepTracks.TrajectoryInEvent = cms.bool(True);
    process.initialStepTracksPreSplitting.TrajectoryInEvent = cms.bool(True);
    process.jetCoreRegionalStepTracks.TrajectoryInEvent = cms.bool(True);
    process.lowPtTripletStepTracks.TrajectoryInEvent = cms.bool(True);
    process.highPtTripletStepTracks.TrajectoryInEvent = cms.bool(True);
    process.lowPtQuadStepTracks.TrajectoryInEvent = cms.bool(True);
    process.pixelPairStepTracks.TrajectoryInEvent = cms.bool(True);
    process.detachedTripletStepTracks.TrajectoryInEvent = cms.bool(True);
    process.mixedTripletStepTracks.TrajectoryInEvent = cms.bool(True);
    process.pixelLessStepTracks.TrajectoryInEvent = cms.bool(True);
    process.tobTecStepTracks.TrajectoryInEvent = cms.bool(True);
    process.muonSeededTracksInOut.TrajectoryInEvent = cms.bool(True);
    process.muonSeededTracksOutIn.TrajectoryInEvent = cms.bool(True);
    process.detachedQuadStepTracks.TrajectoryInEvent = cms.bool(True);

### Main Analyzer
process.analysis = cms.EDAnalyzer('TrackerDpgAnalysis',
                                  ClustersLabel = cms.InputTag("siStripClusters"),
                                  TracksLabel   = cms.VInputTag( cms.InputTag("generalTracksFromRefit")), 
                                  vertexLabel   = cms.InputTag('offlinePrimaryVertices'),
                                  beamSpotLabel = cms.InputTag('offlineBeamSpot'),
                                  DeDx1Label    = cms.InputTag('dedxHarmonic2'),
                                  DeDx2Label    = cms.InputTag('dedxTruncated40'),
                                  DeDx3Label    = cms.InputTag('dedxMedian'),
                                  L1Label       = cms.InputTag('gtDigis'),
                                  HLTLabel      = cms.InputTag("TriggerResults::HLT"),
                                  InitalCounter = cms.uint32(1),
                                  keepOntrackClusters  = cms.untracked.bool(True),
                                  keepOfftrackClusters = cms.untracked.bool(False),
                                  keepTracks           = cms.untracked.bool(True),
                                  keepVertices         = cms.untracked.bool(True),
                                  keepEvents           = cms.untracked.bool(True),
                                  DelayFileNames = cms.untracked.vstring(
                                      "TI_27-JAN-2010_2_delayStep_"+str(options.delayStep)+".xml",
                                      "TO_30-JUN-2009_1_delayStep_"+str(options.delayStep)+".xml",
                                      "TP_09-JUN-2009_1_delayStep_"+str(options.delayStep)+".xml",
                                      "TM_09-JUN-2009_1_delayStep_"+str(options.delayStep)+".xml",
                                  ),
)



if  options.isRawDAQFile or options.isRawEDMFile: 
    process.analysis.TracksLabel = cms.VInputTag(cms.InputTag("generalTracks"))

if options.jsonFile: ### to be checked / created by hand when the runs are take
    process.source.lumisToProcess = LumiList.LumiList(filename = options.jsonFile).getVLuminosityBlockRange()

if options.inputDirectory != "":
    for string in range(len(process.analysis.DelayFileNames)) :
        process.analysis.DelayFileNames[string] = options.inputDirectory+"/"+process.analysis.DelayFileNames[string];
    if options.jsonFile: ### to be checked / created by hand when the runs are take
        process.source.lumisToProcess = LumiList.LumiList(filename = options.inputDirectory+"/"+options.jsonFile).getVLuminosityBlockRange()


process.edTask = cms.Task()

### Define the Path to be executed
if options.isRawDAQFile or options.isRawEDMFile:

    process.doAlldEdXEstimators += process.dedxMedian

    process.edTask.add(process.globalreco_trackingTask);
    process.edTask.add(process.localrecoTask);
    process.edTask.add(process.RawToDigiTask);
    process.edTask.add(process.dedxMedian);

    process.p = cms.Path(
        process.hltfiter*  ## HLT skim                                                                                                                                            
        process.analysis,process.edTask)

else: ## in this case one just need to re-fit the tracks

    process.edTask.add(process.generalTracksFromRefit);
    process.edTask.add(process.TrackRefitter);
    process.edTask.add(process.dedxMedian);
    process.edTask.add(process.doAlldEdXEstimatorsTask);

    process.p = cms.Path(
        process.hltfiter*  ## HLT skim
        process.analysis,process.edTask)

process.schedule = cms.Schedule(process.p);


from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from Configuration.DataProcessing.Utils import addMonitoring
process = addMonitoring(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
