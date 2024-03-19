from CRABClient.UserUtilities import config

config = config()

#### Basic crab config structure
config.General.requestName  = 'delayscan_Express_Run355206'
config.General.transferLogs = False
config.General.transferOutputs = True


config.JobType.pluginName       = 'Analysis'
config.JobType.psetName         = '../../test/trackerdpganalysis_cfg.py'
config.JobType.pyCfgParams      = ['nThreads=4','isRawEDMFile=False','globalTag=run3_data_express','inputDelayFile=TrackerDealyMap_Run355206_pll.csv','triggerList=HLT_ExpressMuons_v1,HLT_IsoMu24_v13,HLT_Physics_v7,HLT_ZeroBias_IsolatedBunches_v5,HLT_ZeroBias_FirstCollisionAfterAbortGap_v5,HLT_ZeroBias_v6']
config.JobType.outputFiles      = ['trackerDPG.root']
config.JobType.numCores         = 4
config.JobType.inputFiles       = ['TrackerDealyMap_Run355206_pll.csv']

config.Data.allowNonValidInputDataset = True
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 4500
config.JobType.numCores    = 4
config.Data.inputDataset  = '/ExpressPhysics/Run2022B-Express-v1/FEVT' 
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'LumiBased'
config.Data.unitsPerJob   = 3
config.Data.outLFNDirBase = '/store/group/dpg_tracker_strip/tracker/Online/RandomDelayScan/Run355206/StreamExpressPhysics'
config.Data.publication   = False
config.Data.lumiMask      = './json_355206.json'

config.Site.storageSite = 'T2_CH_CERN'
