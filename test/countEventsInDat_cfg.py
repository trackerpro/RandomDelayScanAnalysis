import FWCore.ParameterSet.Config as cms
### CMSSW command line parameter parser                                                                                                                                        
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register ('nThreads',1,VarParsing.multiplicity.singleton,VarParsing.varType.int,
                  'Number of threads in the processing')
options.parseArguments()

process = cms.Process("Counter")

process.source = cms.Source("NewEventStreamFileReader",
                            fileNames = cms.untracked.vstring(options.inputFiles),
                            inputFileTransitionsEachEvent = cms.untracked.bool(True)
                            )

process.count = cms.EDFilter("countEvents");
process.path = cms.Path(process.count);

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(options.nThreads),
    numberOfStreams = cms.untracked.uint32(options.nThreads))
