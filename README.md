# Random Delay Scan Analysis 

## CMSSW Setup:

```sh
## load the CMSSW singularity for el8
cmssw-el8
## checkout latest cmssw release used in p5
cmsrel CMSSW_14_0_1
cd CMSSW_14_0_1/src
cmsenv
git-cms-init
git-cms-addpkg DataFormats/SiStripCommon
git-cms-addpkg EventFilter/SiStripRawToDigi
## add tracker pro cmssw forking as we need some specific fixes implemented in a couple of data formats
git remote add trackerpro-cmssw git@github.com:trackerpro/cmssw.git
git fetch trackerpro-cmssw
git merge trackerpro-cmssw/random_delay_scan_update
## clone the analysis package
git clone git@github.com:trackerpro/RandomDelayScanAnalysis.git TrackerDAQAnalysis/RandomDelayScanAnalysis -b automatic_delay_scan
scram b -j 4				 
```

## Run on RAW EDM files:

The analysis can be run on the RAW EDM files produced, after repacking, at T0. These are usually made available on DBS after few days the run has been taken. If you run the analysis starting from RAW, please be aware that local-reco + track-reco needs to be re-done, so jobs are slower than those running on FEVT events (see next section). However, since you process `MinimumBias` or `ZeroBias` events, you will have a large track statistics.

Example on how to run on some input files. 

For Commissioning2021, the following queries to DBS can be performed to access some RAW files for a given Run:  
```sh
voms-proxy-init -voms cms
dasgoclient --query "dataset=/*MinimumBias*/*2021*/*RAW*"
asgoclient --query "file dataset=/MinimumBias0/Commissioning2021-v1/RAW run=346446"
asgoclient --query "site dataset=/MinimumBias0/Commissioning2021-v1/RAW"
```

The ntuplizer can be run on RAW filess as follows:
```sh
cd test
cmsRun trackerdpganalysis_cfg.py globalTag=run3_data nThreads=2 inputDelayFile=TrackerDealyMap_Run346446_pll.csv inputFiles=<file location from dbs> triggerList="HLT_PixelClusters_WP2_ZeroBias*","HLT_L1ETT_ZeroBias*" maxEvents=100 isRawEDMFile=True
```
* You need to copy from P5 the PLL csv file containing the intiial random delay settings
* The trigger list is used to select events that are passing a specific trigger. It is up to you to decide what you like most (if any)
* The cmssw config takes in input also a json file that can be use to filter runs and lumi-sections. Also an event list or range can be provided. Please look inside it to the specific options.

Now you have two alternatives:

* if **RAW** files are **on DSIK**, you can just process them via **crab** using a config file similar to that reported here [crab config](./crab/2022/crabConfig_ExpressPhysics_355206.py).
* Else if **RAW** files are **on TAPE** you could either wait for crab to perform the **TAPE RECALL** for you or you can try to find the file (if present) on EOS CERN and process them via **lxbatch (CondorHT)**. Please check if files for the run/runs of interest are in `/eos/cms/tier0/store/data/` and, if so, please run:
```sh
cd python
python3 submitBatchJobs.py --inputDIR <folder with input files> --outputDIR <output dir on EOS CERN> --jsonFile json_346446.json --eventsPerJob <n events> --triggerList <list of triggers> --globalTag run3_data --delayFileDirectory ../crab/2022/ --jobDIR <name> --submit
```

## Run on FEVT files from ExpressStream:

The analysis can be run on the FEVT files produced by the stream express workflow. Process these files is typically faster than running on RAW, as you have already tracks available in the inputs. However, since a big prescale is applied in impact to the physics stream express, the statistics is more limited so you need to check if you have enough tracks/hits to perform the timing calibration.

For Commissioning2021, the following queries to DBS can be performed:
```sh
dasgoclient --query "file dataset=/ExpressPhysics/Commissioning2021-Express-v1/FEVT run=346446"
dasgoclient --query "site dataset=/ExpressPhysics/Commissioning2021-Express-v1/FEVT"
```

The ntuplizer can be run on FEVT filess as follows:
```sh
cd test
cmsRun trackerdpganalysis_cfg.py globalTag=run3_data_express nThreads=2 inputDelayFile=TrackerDealyMap_Run346446_pll.csv inputFiles=<file location from dbs> triggerList="HLT_HcalNZS*","HLT_L1ETT_ZeroBias*","HLT_PixelClusters_WP1_ZeroBias*" maxEvents=100 isRawEDMFile=False
```
* You need to copy from P5 the PLL csv file containing the intiial random delay settings
* The trigger list is used to select events that are passing a specific trigger. It is up to you to decide what you like most (if any)
* The cmssw config takes in input also a json file that can be use to filter runs and lumi-sections. Also an event list or range can be provided. Please look inside it to the specific options.

Now you have two alternatives:

* if **RAW** files are **on DSIK**, you can just process them via **crab** using a config file similar to that reported here [crab config](./crab/2022/crabConfig_ExpressPhysics_355206.py).
* Else if **RAW** files are **on TAPE** you could either wait for crab to perform the **TAPE RECALL** for you or you can try to find the file (if present) on EOS CERN and process them via **lxbatch (CondorHT)**. Please check if files for the run/runs of interest are in `/eos/cms/tier0/store/express/` and, if so, please run:
```sh
cd python
python3 submitBatchJobs.py --inputDIR <folder with input files> --outputDIR <output dir on EOS CERN> --jsonFile json_346446.json --eventsPerJob <n events> --triggerList <list of triggers> --globalTag run3_data --delayFileDirectory ../crab/2022/ --jobDIR <name> --submit
```
    
## Event Skim:

* Run Locally:
```sh
cd macros/makeSkimTrees/
root -l -b -q skimTrees.C\(input file>, <outputfile>, <isBon = rule the selection string written inside the code>`)
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
root -l -b -q TreeMerge\(<destination file>, < one reference to take the readoutMap>, <directory where all the single files are located>, <if you want to cancel single inputs after merging\)
```

## Build readoutmap with PSU information:

```sh
cd macros/makePSUMap;
root -l -b -q makePSUDetIdMap.C\("<input merged file from which take the readoutmap>", "<PSU-DCU association file>", "<output readout map>"\)
```
 
## Delay analysis per layer/ring/partition:

* The code produces a set of root files with canvases for the different plots: profile fits and best delay vs partition, rings or layers (TEC divide by thin and thick sensors). Please look inside the code for the default values.

```sh
cd macros/makeDelayAnalysis
root -l -b -q delayValidation.C\(<input DIR with all files>, <no correction file stored in ../data/nocorrection.root>, <name to grep for files in inputDIR>, <observable name (branch)>, <plotParitions: to analyze the delay per partions>, <plotLayers: to find the delay per layer>, <plotSlices: to find the best delay per slices>, <outputDIR: name and path of the output directory>)
```    			    	      

## Delay analysis per module:

* Run the delay analysis over a set of different runs with different random delay configuration (fine time calibration per module). Nota bene: the output tree is build from the readoutMap which basic entities are APV chips, hence the number of entries in the final tree will be equal to the number of APVs and not the number of modules. APV within the same module, detid or dcuId, will have the same delay correction to be applied on the PLL.

```sh
cd macros/makeDelayAnalysis
root -l -b -q delayValidationPerModule.C\(<input directory where all files for different runs are located>,<no correction file stored in ../data/nocorrection.root>,<postfix: substring to be find to be sure to run on the input files>, <observable name (branch name)>, <outputDIR: name and path of the output directory>,<saveMeanCanvas: store some gaussian fits of mean charge vs delay>, <saveMPVCanvas: save MPV fit canvases vs delay >, <saveCorrectionTree: to save the delay per channel in a TTree format. Can be analyzed then through the tkCommissioner>
```

### Produce the final correction file

* The code produces a lighter output file with a tree called again delayCorrection. It contains three branches: Detid, fedCh and the correction in units of 25/24ns, to be propagated to the DB. There are two macros, one to produce layer based corrections to be run on the output produced by delayValidation.C, the other produces module based corrections to be run on the output of delayValidationPerModule.C. Example provided only for the second one. The code also eliminates degenerancy in the detid.

```sh
cd macros/makeCorrectionFile
root -l -b -q delayCorrectionPerModule.C(<file name containing the delay tree tree>, <output directory>, <outout file name containing corrections>, <saveFits: save canvas for % reductionFactor modules showing charge vs delay>, <delayCutForPlotting: max delay in order to display large correction modules on the trackermap>, <observable>);
```

### Other macros and codes:

There are other subfolders with some additional macros in order to:
* `macros/makeChannelChecks/`: make some checks on channels for specific FEDs
* `macros/makePlotsForApproval`: make some pretty plot for approval reasons
* `macros/makePSUMap/`: make the PSU / DCU map if there are updates at P5
* `macros/makeFitChargeDistribution/`: fit the charge distribution per module or per layer to define the best parametrization used then in the delay analysis

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
