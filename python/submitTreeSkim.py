### A priori assumption: run on cern batch system with files stored on EOS
### To run: python scripts/submitTreeSkim.py --inputDIR /store/user/rgerosa/TrackerDAQ/DELAYSCAN/ExpressPhysics/crab_delaySCAN_Express_Run271332/ --outputDIR /store/user/rgerosa/TrackerDAQ/DELAYSCAN/ExpressPhysics/SkimmedTrees/Run271332 --outputBaseName tree --batchMode --jobDIR JOB_Run_271332 --queque 1nh --submit
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
#            Job steering                  #                                                                                                                                   
############################################                                                                                                                                    
parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                 
parser.add_option('--inputDIR',       action="store", type="string", dest="inputDIR",     default="",   help="input directory where files are contained")
parser.add_option('--outputDIR',      action="store", type="string", dest="outputDIR",    default="",   help="output DIR")
parser.add_option('--outputBaseName', action="store", type="string", dest="outputBaseName", default="", help="output base name")
parser.add_option('--eventSelection', action="store", type="string", dest="eventSelection", default="", help="event based selection")
parser.add_option('--vertexSelection',action="store", type="string", dest="vertexSelection", default="", help="vertex based selection")
parser.add_option('--trackSelection', action="store", type="string", dest="trackSelection", default="", help="track based selection")
parser.add_option('--clusterSelection', action="store", type="string", dest="clusterSelection", default="", help="cluster based selection")
parser.add_option('--nthreads',       action="store", type=int,      dest="nthreads",    default=4,   help="number of threads for each job")
parser.add_option('--filesPerJob',    action="store", type=int,      dest="filesPerJob", default=5,   help="number of files for each job")
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",  default="",        help="directory for job")
parser.add_option('--queque',       action="store", type="string", dest="queque",  default="longlunch",        help="queque for LSF")
parser.add_option('--submit',       action="store_true",           dest="submit",                     help="submit")

(options, args) = parser.parse_args()

if __name__ == '__main__':


   print ("################################");
   print ("##### Start job submission #####");
   print ("################################");
   
   currentDIR = os.getcwd();
   ## generate binary file                                                                                                                                                      
   ROOT.gROOT.ProcessLine(".L ../macros/makeSkimTrees/skimTrees.C");

   out = subprocess.run(["rm -r "+options.jobDIR,"mkdir -p "+options.jobDIR],shell=True);
   
   ## make the file list ... typically all the files of a given run
   fileList = [];
   for line in glob.globa(options.inputDIR+"/*root",recursive=True)
      if line == "": continue;
      if ".root" in line:
         line = line.replace('\n','');
         fileList.append(line);
   out = subprocess.run(["rm file_temp.txt",shell=True);
   print ("----> Found",len(fileList),"inside the input directory");

   #### given n-files per job calculate and create jobs
   if float(len(fileList)/options.filesPerJob).is_integer():
      njobs = int(len(fileList)/options.filesPerJob);
   else:
      njobs = int(len(fileList)/options.filesPerJob)+1;   
   print ("----> Will create and submit",njobs,"jobs");
   
   jobscript = open("%s/condor_job.sh"%(options.jobDIR),"w");
   jobscript.write("#!/bin/bash\n");
   jobscript.write('cd '+currentDIR+'\n');
   jobscript.write('eval `scramv1 runtime -sh;`\n')
   jobscript.write("cd -;\n");
   jobscript.write("\n");

   #### create one macro per job
   for ijob in range(njobs):

      jobfilelist = open('%s/inputfile_job_%d.list'%(options.jobDIR,ijob),'w')

      for ifile in range(ijob*options.filesPerJob,min(len(fileList),ijob*options.filesPerJob+options.filesPerJob)):         
         nameList = fileList[ifile].split("/");
         while True:
            if nameList[len(nameList)-1] != '':
               break
            else:
               nameList.pop();
         fileName = nameList[len(nameList)-1];
         jobfilelist.write("%s \n"%(fileName));
      
      ## write job sh file                                                                                                                                             
      jobmacro = open('%s/job_%d.C'%(options.jobDIR,ijob),'w')
      jobmacro.write("{\n");
      jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/../macros/makeSkimTrees/skimTrees.C\");\n");      
      if options.eventSelection and options.vertexSelection and options.trackSelection and options.options.clusterSelection:
         jobmacro.write("gROOT->ProcessLine(\""+"skimTrees(\\\"inputfile_job_%d.list\\\",\\\"%s\\\",%i,\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\")\");\n"%(ijob,options.outputBaseName+"_"+str(ijob)+".root",options.nthreads,options.eventSelection,options.vertexSelection,options.trackSelection,options.clusterSelection));
      elif options.eventSelection and options.vertexSelection and options.trackSelection and not options.clusterSelection:
         jobmacro.write("gROOT->ProcessLine(\""+"skimTrees(\\\"inputfile_job_%d.list\\\",\\\"%s\\\",%i,\\\"%s\\\",\\\"%s\\\",\\\"%s\\\")\");\n"%(ijob,options.outputBaseName+"_"+str(ijob)+".root",options.nthreads,options.eventSelection,options.vertexSelection,options.trackSelection));
      elif options.eventSelection and options.vertexSelection and not options.trackSelection and not options.clusterSelection:
         jobmacro.write("gROOT->ProcessLine(\""+"skimTrees(\\\"inputfile_job_%d.list\\\",\\\"%s\\\",%i,\\\"%s\\\",\\\"%s\\\")\");\n"%(ijob,options.outputBaseName+"_"+str(ijob)+".root",options.nthreads,options.eventSelection,options.vertexSelection));
      elif options.eventSelection and options.vertexSelection:
         jobmacro.write("gROOT->ProcessLine(\""+"skimTrees(\\\"inputfile_job_%d.list\\\",\\\"%s\\\",%i,\\\"%s\\\",\\\"%s\\\")\");\n"%(ijob,options.outputBaseName+"_"+str(ijob)+".root",options.nthreads,options.eventSelection,options.vertexSelection));
      elif options.eventSelection:         
         jobmacro.write("gROOT->ProcessLine(\""+"skimTrees(\\\"inputfile_job_%d.list\\\",\\\"%s\\\",%i,\\\"%s\\\")\");\n"%(ijob,options.outputBaseName+"_"+str(ijob)+".root",options.nthreads,options.eventSelection));
      else:
         jobmacro.write("gROOT->ProcessLine(\""+"skimTrees(\\\"inputfile_job_%d.list\\\",\\\"%s\\\",%i)\");\n"%(ijob,options.outputBaseName+"_"+str(ijob)+".root",options.nthreads));
      jobmacro.write("}\n");
      jobmacro.close();

      jobscript.write("if [ $1 -eq "+str(ijob)+" ]; then\n")
      for ifile in range(ijob*options.filesPerJob,min(len(fileList),ijob*options.filesPerJob+options.filesPerJob)):         
         jobscript.write("xrdcp -f "+fileList[ifile]+" ./\n")
      jobscript.write('scp '+currentDIR+'/%s/job_%d.C ./ \n'%(options.jobDIR,ijob))
      jobscript.write('scp '+currentDIR+'/%s/inputfile_job_%d.list ./ \n'%(options.jobDIR,ijob))
      jobscript.write('root -l -b -q job_%d.C\n'%(ijob));
      jobscript.write("mkdir -p "+options.outputDIR+"\n");
      jobscript.write("xrdcp -f "+options.outputBaseName+"_"+str(ijob)+".root "+options.outputDIR+"/\n");       
      jobscript.write("fi\n");

   jobscript.write("\n");
   jobscript.close();
   out = subprocess.run('chmod a+x %s/condor_job.sh'%(options.jobDIR),shell=True);
   
   condor_job = open("%s/condor_job.sub"%(options.jobDIR),"w");
   condor_job.write("executable = %s/%s/condor_job.sh\n"%(currentDIR,options.jobDIR));
   condor_job.write("arguments = $(ProcId)\n");
   condor_job.write("output = %s/%s/condor_job_$(ProcId).out\n"%(currentDIR,options.jobDIR));
   condor_job.write("error  = %s/%s/condor_job_$(ProcId).err\n"%(currentDIR,options.jobDIR));
   condor_job.write("log    = %s/%s/condor_job_$(ProcId).log\n"%(currentDIR,options.jobDIR));
   condor_job.write("transfer_output_files=\"\"\n");
   condor_job.write("when_to_transfer_output = ON_EXIT\n");
   condor_job.write("should_transfer_files   = YES\n");
   condor_job.write("universe = vanilla\n");
   condor_job.write("+JobFlavour = \""+options.queque+"\"\n");
   condor_job.write("queue "+str(njobs)+"\n");
   condor_job.write("request_cpus   = "+str(options.nthreads)+"\n");
   condor_job.write("request_memory = "+str(2000*options.nthreads)+"\n");
   condor_job.close();
      
   if options.submit:
      out = subprocess.run("condor_submit %s/condor_job.sub"%(options.jobDIR),shell=True);
