import os
import glob
import math
from array import array
import sys
import time
import subprocess
import ROOT

from optparse import OptionParser
from subprocess import Popen

############################################                                                                                                                                 
parser = OptionParser()
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",      default="", help="to be used to correctly set the working area")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="", help="output dir where workspaces are copied")
parser.add_option('--delayFile', action="store", type="string", dest="delayFile", default = "", help="Location of the delay file")
parser.add_option('--jsonFile',     action="store", type="string", dest="jsonFile",      default="", help="json file to be applied")
parser.add_option('--nThreads',     action="store", type=int,      dest="nThreads",      default=2,  help="number of threads")
parser.add_option('--filesPerJob',  action="store", type=int,      dest="filesPerJob",   default=1,  help="number of files to process in each --> to be used with small files")
parser.add_option('--eventsPerJob', action="store", type=int,      dest="eventsPerJob",  default=1000, help="number of events to process in each job --> to be used with big files")
parser.add_option('--isRawEDMFile',    action="store_true", dest="isRawEDMFile", help="isRawEDMFile")
parser.add_option('--triggerList',  action="store", type="string", dest="triggerList",   default="HLT_ZeroBias*", help="list of triggers to be requested")
parser.add_option('--globalTag',    action="store", type="string", dest="globalTag",     default="run3_data", help="global tag to be used in the job")
parser.add_option('--submit',       action="store_true", dest="submit",    help="submit")
parser.add_option('--jobDIR',       action="store",      type="string", dest="jobDIR", default="",  help="directory for job")
parser.add_option('--queque',       action="store",      type="string", dest="queque", default="longlunch",  help="queque for HTCondor")

(options, args) = parser.parse_args()
############################################                                                                                                                                  

if __name__ == '__main__':

  currentDIR = os.getcwd();
  
  listOfFiles = [];
  for item in glob.glob(options.inputDIR+"/*root", recursive=True):
    if item == "" or item == "\n": continue;
    if not ".root" in item: continue;
    listOfFiles.append(item.replace("\n",""));
    
  print ("number of files ffound : ",len(listOfFiles));

  out = subprocess.run(['rm file_delay.temp','mkdir -p '+options.jobDIR,'mkdir -p '+options.outputDIR], shell=True);

  ## open each file and split according to eventsPerJob
  njobs  = 0;
  file_list    = [];
  event_list   = [];
  event_counter = 0;

  for ifile,filename in enumerate(listOfFiles):      

    if len(file_list) <= njobs:
      file_list.append([]);
      file_list[njobs].append("file:"+filename);
    else:        
      file_list[njobs].append("file:"+filename);


    tfile = ROOT.TFile.Open(filename);
    tree  = tfile.Get("Events");
    nevents = tree.GetEntries();
    ## loop on events and save the starting one
    if nevents == 0: continue;
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus("EventAuxiliary",1);
    for event in range(0,nevents):        
      tree.GetEntry(event);
      if event_counter % options.eventsPerJob == 0:
        event_list.append({"run":  tree.EventAuxiliary.run(),  
                           "lumi":  tree.EventAuxiliary.luminosityBlock(),  
                           "event":  tree.EventAuxiliary.event()});
        njobs = njobs + 1;
        if len(file_list) <= njobs:
          file_list.append([]);
          file_list[njobs].append("file:"+filename);
        else:        
          file_list[njobs].append("file:"+filename);
        event_counter = event_counter + 1;

    if ifile == len(listOfFiles)-1 and event_counter % options.eventsPerJob != 0:
      tree.GetEntry(nevents);
      event_list.append({"run":  tree.EventAuxiliary.run(),
                         "lumi":  tree.EventAuxiliary.luminosityBlock(),
                         "event":  tree.EventAuxiliary.event()});        
      njobs = njobs + 1;


  print("Creating njbos ",njobs);
  print(event_list);
  print(file_list);

  jobscript = open("%s/condor_job.sh"%(options.jobDIR),"w");
  jobscript.write("#!/bin/bash\n");
  jobscript.write('cd '+currentDIR+'\n');
  jobscript.write('eval `scramv1 runtime -sh;`\n')
  jobscript.write("cd -;\n");  
  jobscript.write('scp '+currentDIR+'/../test/trackerdpganalysis_cfg.py ./ \n');
  jobscript.write('scp '+currentDIR+"/"+options.delayFile+'./ \n');
  jobscript.write('scp '+currentDIR+"/"+options.jsonFile+'./ \n');

  for ijob in range(0,njobs-1):
    jobscript.write("if [ $1 -eq "+str(ijob)+" ]; then\n")
    
    inputFiles = "";
    for ientry,entry in enumerate(file_list[ijob]):
      if ientry < len(file_list[ijob])-1:
        inputFiles += str(entry)+",";
      else:
        inputFiles += str(entry)+"";

    if options.isRawEDMFile:
      jobscript.write('cmsRun trackerdpganalysis_cfg.py inputFiles='+inputFiles+' ouputFileName=trackerDPG_'+str(ijob)+".root jsonFile="+options.jsonFile+" globalTag="+options.globalTag+" nThreads="+str(options.nThreads)+" inputDelayFile="+options.delayFile+" triggerList="+options.triggerList+" eventRange="+str(event_list[ijob]["run"])+":"+str(event_list[ijob]["lumi"])+":"+str(event_list[ijob]["event"])+"-"+str(event_list[ijob+1]["run"])+":"+str(event_list[ijob+1]["lumi"])+":"+str(event_list[ijob+1]["event"])+" isRawEDMFile=True\n");
    else:
      jobscript.write('cmsRun trackerdpganalysis_cfg.py inputFiles='+inputFiles+' ouputFileName=trackerDPG_'+str(ijob)+".root jsonFile="+options.jsonFile+" globalTag="+options.globalTag+" nThreads="+str(options.nThreads)+" inputDelayFile="+options.delayFile+" triggerList="+options.triggerList+" eventRange="+str(event_list[ijob]["run"])+":"+str(event_list[ijob]["lumi"])+":"+str(event_list[ijob]["event"])+"-"+str(event_list[ijob+1]["run"])+":"+str(event_list[ijob+1]["lumi"])+":"+str(event_list[ijob+1]["event"])+" isRawEDMFile=False\n");
    jobscript.write("/usr/bin/eos mkdir -p "+options.outputDIR+"\n");
    jobscript.write("xrdcp -f trackerDPG_"+str(ijob)+".root "+options.outputDIR+"/\n");       
    jobscript.write("fi\n")

  jobscript.write("\n");
  jobscript.close();
  out = subprocess.run('chmod a+x %s/condor_job.sh'%(options.jobDIR),shell=True);
  
  condor_job = open("%s/condor_job.sub"%(options.jobDIR),"w");
  condor_job.write("executable = %s/%s/condor_job.sh\n"%(currentDIR,options.jobDIR));
  condor_job.write("arguments = $(ProcId)\n");
  condor_job.write("requestcpus = "+str(options.nThreads)+"\n");
  condor_job.write("output = %s/%s/condor_job_$(ProcId).out\n"%(currentDIR,options.jobDIR));
  condor_job.write("error  = %s/%s/condor_job_$(ProcId).err\n"%(currentDIR,options.jobDIR));
  condor_job.write("log    = %s/%s/condor_job_$(ProcId).log\n"%(currentDIR,options.jobDIR));
  condor_job.write("transfer_output_files=\"\"\n");
  condor_job.write("when_to_transfer_output = ON_EXIT\n");
  condor_job.write("should_transfer_files   = YES\n");
  condor_job.write("universe = vanilla\n");
  condor_job.write("+JobFlavour = \""+options.queque+"\"\n");
  condor_job.write("queue "+str(njobs-1)+"\n");
  condor_job.close();
      
  if options.submit:
    out = subprocess.run("condor_submit %s/condor_job.sub"%(options.jobDIR),shell=True);
