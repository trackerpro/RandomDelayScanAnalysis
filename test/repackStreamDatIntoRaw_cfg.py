import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('nThreads',1,VarParsing.multiplicity.singleton,VarParsing.varType.int,
                  'Number of threads in the processing')

options.parseArguments()

process = cms.Process("Convert")
process.load("DQM.SiStripCommon.MessageLogger_cfi")
 
process.source = cms.Source("NewEventStreamFileReader",
                            fileNames = cms.untracked.vstring(options.inputFiles),
                            inputFileTransitionsEachEvent = cms.untracked.bool(True)
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(options.nThreads),
    numberOfStreams = cms.untracked.uint32(options.nThreads))

maxEvents = cms.PSet(
    input = cms.untracked.int32(-1)
 )

process.anal = cms.EDAnalyzer("EventContentAnalyzer")
 
process.out = cms.OutputModule("PoolOutputModule",
 fileName = cms.untracked.string("/tmp/output.root"),
 outputCommands = cms.untracked.vstring("drop *", "keep FEDRawDataCollection_*_*_*")
)

process.p = cms.Path(process.anal)
process.e = cms.EndPath(process.out)
