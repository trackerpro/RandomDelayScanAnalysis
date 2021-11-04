##################################
### Random Delay Scan Analysis ###
##################################

### CMSSW Setup:

    ```sh
    cmsrel CMSSW_12_0_3_patch1 ;
    cd CMSSW_12_0_3_patch1/src ;
    cmsenv;		      
    git-cms-init; 
    git clone git@github.com:trackerpro/RandomDelayScanAnalysis.git TrackerDAQAnalysis/RandomDelayScanAnalysis
    scramv1 b -j 4;					 
    ```
### To run the TrackerDPGAnalyzer locally:

    ```sh 
    cd  TrackerDAQAnalysis/RandomDelayScan/test/
    cmsRun trackerdpganalysis_cfg.py delayStep=<int corresponding to the delay xml chosen>
    ``` 

### Run on t0streamer files

    Example on how to run on a couple of input files
    
    ```sh
    cmsRun trackerdpganalysis_cfg.py isRawDAQFile=True globalTag=run3_data_express nThreads=2 inputDirectory=../crab/2021/ inputFiles=/store/t0streamer/Data/Express/000/346/446/run346446_ls0001_streamExpress_StorageManager.dat delayStep=0 triggerList="HLT_HcalNZS*","HLT_L1ETT_ZeroBias*","HLT_PixelClusters_WP1_ZeroBias*" maxEvents=100
    cmsRun trackerdpganalysis_cfg.py isRawDAQFile=True globalTag=run3_data nThreads=2 inputDirectory=../crab/2021/ inputFiles=/store/t0streamer/Data/PhysicsMinimumBias0/000/346/446/run346446_ls0001_streamPhysicsMinimumBias0_StorageManager.dat triggerList="HLT_PixelClusters_WP2_ZeroBias*","HLT_L1ETT_ZeroBias*" maxEvents=100
    ```

    Trigger that can be requried depending on the output stream can be accessed via OMS for each run number. Above examples are taken for run 346446 taken in fall 2021 at Run3 commissioning startup.

    Submit jobs via lxbatch (condorHT) on these streamer files are not published in crab:

    ```sh
    ```	

### Run on RAW files
    
    Example on how to run on some input files. For Commissioning2021, the following queries to DBS can be performed:
    
    ```sh
     dasgoclient --query "dataset=/*MinimumBias*/*2021*/*RAW*"
     dasgoclient --query "file dataset=/MinimumBias0/Commissioning2021-v1/RAW run=346446"
     dasgoclient --query "site dataset=/MinimumBias0/Commissioning2021-v1/RAW"
    ```

    If files are on-disk, you can use the following commands:

    ```sh
    cmsRun trackerdpganalysis_cfg.py isRawEDMFile=True globalTag=run3_data nThreads=2 inputDirectory=../crab/2021/ inputFiles=<file location> triggerList="HLT_PixelClusters_WP2_ZeroBias*","HLT_L1ETT_ZeroBias*" maxEvents=100
    ```

    Submit jobs via lxbatch (condorHT) on these streamer files are not published in crab:


### Run on FEVT files
    
    Example on how to run on some input files. For Commissioning2021, the following queries to DBS can be performed:

    ```sh
     dasgoclient --query "file dataset=/ExpressPhysics/Commissioning2021-Express-v1/FEVT run=346446"
     dasgoclient --query "site dataset=/ExpressPhysics/Commissioning2021-Express-v1/FEVT"
    ```

    If files are on-disk, you can use the following commands:

    ```sh
    cmsRun trackerdpganalysis_cfg.py globalTag=run3_data_express nThreads=2 inputDirectory=../crab/2021/ inputFiles=/store/express/Commissioning2021/ExpressPhysics/FEVT/Express-v1/000/346/446/00000/0455c1c4-357f-434c-951f-ab6f2ba98683.root delayStep=0 triggerList="HLT_HcalNZS*","HLT_L1ETT_ZeroBias*","HLT_PixelClusters_WP1_ZeroBias*" maxEvents=100
    ```
    
    To run crab jobs: 
    
    TrackerDAQAnalysis/RandomDelayScanAnalysis/crab/json*json = example of json file to select a particular run. it can be obtained via: ```sh dasgoclient --query "run,lumi dataset=/ExpressPhysics/Commissioning2021-Express-v1/FEVT run=346446" ```
    TrackerDAQAnalysis/RandomDelayScanAnalysis/crab/crab_*py = example of json crab config to analyze the delay run


    
### Event Skim:

Run Locally:

    ```sh
    cd TrackerDAQAnalysis/RandomDelayScan/macros/;
    root -l;
    .L skimTrees.C;
    skimTrees(<input file>, <outputfile>, <isBon = rule the selection string written inside the code>);
    ```

    After running jobs, the output files belonging to a single run must be skimmed applying the desired analysis selection and dropping tracks/vertxes and event information. Once these selections are applied, we just need the delayMap and the clusters tree for the offiline analysis.

    CERN CondorHT submission submission (assumes files are stored on cern EOS):
    
    ```sh
    cd TrackerDAQAnalysis/RandomDelayScan ;
    python scripts/submitTreeSkim.py  --inputDIR <directory with all the files for a given run, produced by crab is ok> --outputDIR <output location on Cern EOS> --outputBaseName <base name for the output root file> --isBOn (in case you want to apply bOn selections) --jobDIR <JOBDIR> --queque <QUEQUE> --submit
    ```

### Copy from EOS to a local machine:

    cd TrackerDAQAnalysis/RandomDelayScan ;
    python scripts/copyFilesEOS.py --inputDIR <directory where skimmed trees are located> --outputDIR <local directory to be copied>

    Once skimmed trees are ready, to copy them to a local machine you could use:
    

### Merge trees belonging to a given run:

    ```sh
    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L TreeMerge.C;
    TreeMerge(<destination file>, < one reference to take the readoutMap>, <directory where all the single files are located>, <if you want to cancel single inputs after merging)
    ```

### Fit Cluster Charge or S/N for each detId:

    ```sh
    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L fitChargeDistribution.C;
    fitChargeDistribution(<merged file>,<output dir>,<observable name (branch name)>, delayMin, delayMax, <apply or not path lenght correction>, <store outputs>)
    ```

    The code produces a root file with some fits, just to see if they are making sense. Then, a text file is produced to be displayed on the tracker map: <detId,peak of landau conv gaussian shape>. To display it: http://test-stripdbmonitor.web.cern.ch/test-stripdbmonitor/PrintTrackerMap/print_TrackerMap.php

 
### Delay analysis per layer/ring/partition:

    ```sh
    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L delayValidation.C;
    delayValidation(<merged file>,<no correction file stored in ../data/nocorrection.root>,<observable name (branch name)>, <plotParitions: to analyze the delay per partions>, <plotLayers: to find the delay per layer>, <plotSlices: to find the best delay per slices>, <outputDIR: name and path of the output directory>)
   ```    			    	      
   The code produces a set of root files with canvases for the different plots: profile fits and best delya vs partition, rings or layers (TEC divide by thin and thick sensors).

### Delay analysis per module

    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L delayValidationPerModule.C;
    delayValidationPerModule(<input directory where all the merged files for different runs are located>,<no correction file stored in ../data/nocorrection.root>,<postfix: substring to be find to be sure to run on the merged files>, <observable name (branch name)>, <outputDIR: name and path of the output directory>,<saveMeanCanvas: store some gaussian fits of mean charge vs delay>, <saveMPVCanvas: save MPV fit canvases vs delay >, <saveCorrectionTree: to save the delay per channel in a TTree format. Can be analyzed then through the tkCommissioner>

Run the delay analysis over a set of different runs with different random delay configuration (fine time calibration per module).

### Delay correction per module

    ```sh
    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L delayCorrectionPerModule.C;
    delayCorrectionPerModule(<input file from delayValidationPerModule.C>,<output directory>,<outout root file name>)
    ```

### Compare different runs:

    ```sh
    cd TrackerDAQAnalysis/RandomDelayScan/macros;
    root -l;
    .L compareClustersVsRun.C;
    compareClustersVsRun(<inputDIR : directory with all the dpg trees for all the runs> <runList: list of run numbers to be considered> <outputDIR: where plots and root files are created> <postfix: string to grep to access to input root files>);
    ```

    The code produces a lighter output file with a tree called again delayCorrection. It contains three branches: Detid, fedCh and te correction in units of 25/24ns, to be propagated to the DB.

