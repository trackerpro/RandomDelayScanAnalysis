import glob
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')
options.register(
'inputFileDir',"",VarParsing.multiplicity.singleton,VarParsing.varType.string,'input directory with files'
)
options.parseArguments()

from Configuration.StandardSequences.Eras import eras
process = cms.Process("fedheaderanalysis",eras.Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

if options.inputFileDir:
    for files in glob.glob(options.inputFileDir+"/*root"):
        options.inputFiles.append("file:"+files);

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True)
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,"auto:run3_data_express");

process.load("EventFilter.SiStripRawToDigi.FedChannelDigis_cfi")
process.FedChannelDigis.UnpackBadChannels = cms.bool(True)
process.FedChannelDigis.DoAPVEmulatorCheck = cms.bool(True)
process.FedChannelDigis.ProductLabel = cms.InputTag('rawDataCollector')
process.FedChannelDigis.LegacyUnpacker = cms.bool(False)

process.stripSummaryAnalyzer = cms.EDAnalyzer("CheckStripEventSummary",
        stripEventSummary = cms.InputTag("FedChannelDigis")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("outputStripSummaryAnalyzer.root")
)

process.path = cms.Path(process.FedChannelDigis+
                        process.stripSummaryAnalyzer)

process.schedule = cms.Schedule(process.path);                        


