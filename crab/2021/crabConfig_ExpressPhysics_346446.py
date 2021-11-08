from CRABClient.UserUtilities import config

config = config()

#### Basic crab config structure
config.General.requestName  = 'delayscan_Express_Run346446'
config.General.transferLogs = False
config.General.transferOutputs = True


config.JobType.pluginName       = 'Analysis'
config.JobType.psetName         = '../../test/trackerdpganalysis_cfg.py'
config.JobType.pyCfgParams      = ['nThreads=4','delayStep=0','globalTag=run3_data_express','inputDirectory=./','triggerList=HLT_HcalNZS*,HLT_L1ETT_ZeroBias*,HLT_PixelClusters_WP1_ZeroBias*']
config.JobType.outputFiles      = ['trackerDPG.root']
config.JobType.numCores         = 4
config.JobType.inputFiles       = ['TI_27-JAN-2010_2_delayStep_0.xml','TM_09-JUN-2009_1_delayStep_0.xml','TO_30-JUN-2009_1_delayStep_0.xml','TP_09-JUN-2009_1_delayStep_0.xml']

config.Data.allowNonValidInputDataset = True
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 4500
config.JobType.numCores    = 4
config.Data.inputDataset  = '/ExpressPhysics/Commissioning2021-Express-v1/FEVT' 
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'LumiBased'
config.Data.unitsPerJob   = 15
config.Data.outLFNDirBase = '/store/group/dpg_tracker_strip/tracker/Online/RandomDelayScan/Run346446/StreamExpressPhysics'
config.Data.publication   = False
config.Data.lumiMask      = '/afs/cern.ch/user/r/rgerosa/work/TrackerDAQ/DELAY_SCAN/CMSSW_12_0_3_patch1/src/TrackerDAQAnalysis/RandomDelayScanAnalysis/crab/2021/json_346446.json'

config.Site.storageSite = 'T2_CH_CERN'
