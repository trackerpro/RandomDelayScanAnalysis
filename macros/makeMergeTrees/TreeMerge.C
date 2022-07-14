// Merge all the trees beloging to a delay setting: destination = merged file, rsources are files to be merged
void TreeMerge(const string & destinationFile, 
	       const string & inputDirectoryPath, 
	       const string & grepNameDir = "", 
	       const bool & cancelInputFiles = false, 
	       const int  & numberOfThreads  = 4) {

  system(("rm -rf "+destinationFile).c_str());

  // make list of all inputFIles before creting the merge oen
  if(not grepNameDir.empty())
    system(("find "+inputDirectoryPath+" | grep root | grep "+grepNameDir+" > file.temp").c_str());
  else
    system(("find "+inputDirectoryPath+" | grep root > file.temp").c_str());

  ifstream infile;
  string line;
  vector<string> fileList;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line != "" and TString(line).Contains(".root"))
	fileList.push_back(line);
    }
  }
  system("rm file.temp");

  // make hadd from a command line
  string inputFiles;
  for(auto file : fileList)
    inputFiles += file+" ";
  cout<<inputFiles<<endl;
  gSystem->Exec(Form("hadd -j %d -f -d %s -O %s %s",numberOfThreads,inputDirectoryPath.c_str(),destinationFile.c_str(),inputFiles.c_str()));

  std::cout << "Reindexing the element of the merged tree..."  << std::endl;
  // just take the first tree of those merged as reference
  string referenceFile = fileList.front();
  TFile* f = TFile::Open(destinationFile.c_str(),"update");
  TFile* f2 = TFile::Open(referenceFile.c_str());
  TTree* tree = (TTree*) f2->FindObjectAny("psumap");

  TTree *newtree = NULL;
  if(tree == 0 or tree == NULL){
    cout<<"[TreeMerge] no psumap file --> please check"<<endl;
  }
  else{
    f->cd();
    newtree = tree->CloneTree();
    newtree->BuildIndex("dcuId"); // deti-id index
    newtree->Write("psumap",TObject::kOverwrite);
  }

  // Look at the redoutMap
  tree = (TTree*)f2->FindObjectAny("readoutMap");
  if(tree == 0 or tree == NULL){
    cout<<"[TreeMerge] no readoutMap file --> please check"<<endl;
  }
  else{  
    f->cd();
    newtree = tree->CloneTree();
    newtree->BuildIndex("detid"); // det-id index
    newtree->Write("readoutMap",TObject::kOverwrite);
    f2->Close();
  }

  //close
  std::cout << "Saving merged file." << std::endl;
  f->Write();
  f->Close();

  //cancel all the single files, keeping only the unmerged one
  if(cancelInputFiles){
    for(auto fileName : fileList){
      cout<<"rm "+inputDirectoryPath+"/"+fileName<<endl;
      system(("rm "+inputDirectoryPath+"/"+fileName).c_str());
    }
  }
}
