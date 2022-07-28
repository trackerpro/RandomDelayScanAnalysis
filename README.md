# Random Delay Scan Analysis 

## CMSSW Setup:

```sh
cmsrel CMSSW_12_4_2 ;
cd CMSSW_12_4_2/src ;
cmsenv;		      
git-cms-init; 
git remote add trackerpro-cmssw git@github.com:trackerpro/cmssw.git;
git fetch trackerpro-cmssw;
git-cms-addpkg DataFormats/SiStripCommon;
git-cms-addpkg EventFilter/SiStripRawToDigi;
git clone git@github.com:trackerpro/RandomDelayScanAnalysis.git TrackerDAQAnalysis/RandomDelayScanAnalysis
scramv1 b -j 4;					 
```

## Run on RAW DAQ files

* Example on how to run on a couple of input files. Location of streamer file is ```/store/t0streamer/Data/``` for the various streams. The job has to run the re-packing, un-packing, local and global tracking steps.
    
```sh
cmsRun trackerdpganalysis_cfg.py isRawDAQFile=True globalTag=run3_data_express nThreads=4 inputDelayFile=TrackerDealyMap_Run346446_pll.csv inputFiles=/store/t0streamer/Data/Express/000/346/446/run346446_ls0001_streamExpress_StorageManager.dat delayStep=0 triggerList="HLT_HcalNZS*","HLT_L1ETT_ZeroBias*","HLT_PixelClusters_WP1_ZeroBias*" maxEvents=100
```

Trigger that can be requried depending on the output stream can be accessed via OMS for each run number. Above examples are taken for run 346446 taken in fall 2021 at Run3 commissioning startup.

* Submit jobs via lxbatch (condorHT) on these streamer files are not published in crab:

```sh
python3 submitBatchJobs.py --inputDIR /eos/cms/store/t0streamer/Data/Express/000/346/446/ --outputDIR /eos/cms/store/group/dpg_tracker_strip/tracker/Online/RandomDelayScan/Run346446/MinimumBias0/ --jsonFile json_346446.json --eventsPerJob 1500 --delayStep 0 --nThreads 2 --triggerList "HLT_HcalNZS*,HLT_L1ETT_ZeroBias*,HLT_PixelClusters_WP1_ZeroBias*,HLT_PixelClusters_WP2_ZeroBias*" --globalTag run3_data_express --delayFileDirectory ../crab/2021/ --jobDIR jobs_minimumbias0 --isRawDAQFile --submit
```	
## Run on RAW EDM files
    
* Example on how to run on some input files. For Commissioning2021, the following queries to DBS can be performed:
    
```sh
dasgoclient --query "dataset=/*MinimumBias*/*2021*/*RAW*"
asgoclient --query "file dataset=/MinimumBias0/Commissioning2021-v1/RAW run=346446"
asgoclient --query "site dataset=/MinimumBias0/Commissioning2021-v1/RAW"
```
	
However typically RAW files are just kept on TAPE or at T0_CERN where jobs cannot be run. However, one can also access directly to the T0 space like in ```/eos/cms/tier0/store/data/Commissioning2021/```

* if files are on-disk, you can use the following commands that performs un-packing, local and global tracking steps:

```sh
cmsRun trackerdpganalysis_cfg.py isRawEDMFile=True globalTag=run3_data nThreads=2 inputDelayFile=TrackerDealyMap_Run346446_pll.csv inputFiles=<file location> triggerList="HLT_PixelClusters_WP2_ZeroBias*","HLT_L1ETT_ZeroBias*" maxEvents=100
```

* Submit jobs via lxbatch (condorHT) on these streamer files are not published in crab:

```sh
python3 submitBatchJobs.py --inputDIR /eos/cms//tier0/store/data/Commissioning2021/MinimumBias0/RAW/v1/000/346/446/ --outputDIR /eos/cms/store/group/dpg_tracker_strip/tracker/Online/RandomDelayScan/Run346446/MinimumBias0/ --jsonFile json_346446.json --eventsPerJob 1500 --delayStep 0 --triggerList "HLT_HcalNZS*,HLT_L1ETT_ZeroBias*,HLT_PixelClusters_WP1_ZeroBias*,HLT_PixelClusters_WP2_ZeroBias*" --globalTag run3_data --delayFileDirectory ../crab/2021/ --jobDIR jobs_minimumbias0 --isRawEDMFile --submit
```

## Run on FEVT files from ExpressStream
    
* Example on how to run on some input files. For Commissioning2021, the following queries to DBS can be performed:

```sh
dasgoclient --query "file dataset=/ExpressPhysics/Commissioning2021-Express-v1/FEVT run=346446"
dasgoclient --query "site dataset=/ExpressPhysics/Commissioning2021-Express-v1/FEVT"
```

* If files are on-disk, you can use the following commands:

```sh
cmsRun trackerdpganalysis_cfg.py globalTag=run3_data_express nThreads=2 inputDelayFile=TrackerDealyMap_Run346446_pll.csv inputFiles=/store/express/Commissioning2021/ExpressPhysics/FEVT/Express-v1/000/346/446/00000/0455c1c4-357f-434c-951f-ab6f2ba98683.root delayStep=0 triggerList="HLT_HcalNZS*","HLT_L1ETT_ZeroBias*","HLT_PixelClusters_WP1_ZeroBias*" maxEvents=100
```
    
* To run crab jobs: 
    
TrackerDAQAnalysis/RandomDelayScanAnalysis/crab/json*json = example of json file to select a particular run. it can be obtained via: ```sh dasgoclient --query "run,lumi dataset=/ExpressPhysics/Commissioning2021-Express-v1/FEVT run=346446" ```
TrackerDAQAnalysis/RandomDelayScanAnalysis/crab/crab_*py = example of json crab config to analyze the delay run
    
## Event Skim:

* Run Locally:

```sh
cd macros/makeSkimTrees/;
root -l;
.L skimTrees.C;
skimTrees(<input file>, <outputfile>, <isBon = rule the selection string written inside the code>);
```

* After running jobs, the output files belonging to a single run must be skimmed applying the desired analysis selection and dropping tracks/vertxes and event information. Once these selections are applied, we just need the delayMap and the clusters tree for the offiline analysis.

* CERN CondorHT submission submission (assumes files are stored on cern EOS):
    
```sh
python scripts/submitTreeSkim.py  --inputDIR <directory with all the files for a given run, produced by crab is ok> --outputDIR <output location on Cern EOS> --outputBaseName <base name for the output root file> --isBOn (in case you want to apply bOn selections) --jobDIR <JOBDIR> --queque <QUEQUE> --submit
```

## Merge trees belonging to a given run (optional):

* Script to automatically merge files into a single ROOT file. The reference file should be one of the files that is going to be merged, from which the readout map is taken.

```sh
cd macros/makeMergeTrees;
root -l;
.L TreeMerge.C;
TreeMerge(<destination file>, < one reference to take the readoutMap>, <directory where all the single files are located>, <if you want to cancel single inputs after merging)
```

## Build readoutmap with PSU information

```sh
cd macros/makePSUMap;
root -l;
.L makePSUDetIdMap.C;
makePSUDetIdMap("<input merged file from which take the readoutmap>", "<PSU-DCU association file>", "<output readout map>");
```

## Fit Cluster Charge or S/N for each detId:

* The code produces a root file with some fits, just to see if they are making sense. Then, a text file is produced to be displayed on the tracker map: <detId,peak of landau conv gaussian shape>. To display it: http://test-stripdbmonitor.web.cern.ch/test-stripdbmonitor/PrintTrackerMap/print_TrackerMap.php

```sh
cd macros/makeFitChargeDistribution/
root -l;
.L fitChargeDistribution.C;
fitChargeDistribution(<merged file>,<output dir>,<observable name (branch name)>, delayMin, delayMax, <apply or not path lenght correction>, <store outputs>)
```
 
## Delay analysis per layer/ring/partition:

* The code produces a set of root files with canvases for the different plots: profile fits and best delya vs partition, rings or layers (TEC divide by thin and thick sensors).

```sh
cd macros/makeDelayAnalysis
root -l;
.L delayValidation.C;
delayValidation(<merged file>,<no correction file stored in ../data/nocorrection.root>,<observable name (branch name)>, <plotParitions: to analyze the delay per partions>, <plotLayers: to find the delay per layer>, <plotSlices: to find the best delay per slices>, <outputDIR: name and path of the output directory>)
```    			    	      

## Delay analysis per module

* Run the delay analysis over a set of different runs with different random delay configuration (fine time calibration per module). Nota bene: the output tree is build from the readoutMap which basic entities are APV chips, hence the number of entries in the final tree will be equal to the number of APVs and not the number of modules. APV within the same module, detid or dcuId, will have the same delay correction to be applied on the PLL.

```sh
cd macros/makeDelayAnalysis
root -l;
.L delayValidationPerModule.C;
delayValidationPerModule(<input directory where all the merged files for different runs are located>,<no correction file stored in ../data/nocorrection.root>,<postfix: substring to be find to be sure to run on the merged files>, <observable name (branch name)>, <outputDIR: name and path of the output directory>,<saveMeanCanvas: store some gaussian fits of mean charge vs delay>, <saveMPVCanvas: save MPV fit canvases vs delay >, <saveCorrectionTree: to save the delay per channel in a TTree format. Can be analyzed then through the tkCommissioner>
```

### Produce the final correction file:

* The code produces a lighter output file with a tree called again delayCorrection. It contains three branches: Detid, fedCh and the correction in units of 25/24ns, to be propagated to the DB. There are two macros, one to produce layer based corrections to be run on the output produced by delayValidation.C, the other produces module based corrections to be run on the output of delayValidationPerModule.C. Example provided only for the second one. The code also eliminates degenerancy in the detid.

```sh
cd macros/makeCorrectionFile
root -l;
.L delayCorrectionPerModule.C;
delayCorrectionPerModule(<file name containing the delay tree tree>, <output directory>, <outout file name containing corrections>, <saveFits: save canvas for % reductionFactor modules showing charge vs delay>, <delayCutForPlotting: max delay in order to display large correction modules on the trackermap>, <observable>);
```

### Upload to P5 database:

* Log into a p5 srv machine as trackerpro, execute ``getConfdb CONFDB_ONLINE``, copy there the text file containing the delay-corrections produced by the previous step, called ``delayCorrection.txt``, and run the following code for the upload as follows:
```sh
cd /opt/cmssw/shifter/erik/AccessDb/Raffaele;
unset LD_LIBRARY_PATH
cd /opt/cmssw/Stable/current/
eval `scramv1 runtime -sh` ;
cd -;
sh compile.sh
python3 createSingleRandomStep.py
``` 
* Before running, edit the ``createSingleRandomStep.py`` in order to give to him, as input, the right text file containg the delays to be applied along with the name of the output summary xml files produced in order to show the values of coarse and fine delays to be applied on the PLL and FED channel registers. 