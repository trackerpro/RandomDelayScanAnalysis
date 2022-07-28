#include "AccessDb.h"


static void getDbConfiguration(string &name,string &passwd,string &path){
  char* cpath=getenv("CONFDB");
  if (cpath == NULL) return;
  string confdb(cpath);
  int ipass = confdb.find("/");
  int ipath = confdb.find("@");
  name.clear();
  name = confdb.substr(0,ipass);
  passwd.clear();
  passwd = confdb.substr(ipass+1,ipath-ipass-1);
  path.clear();
  path = confdb.substr(ipath+1,confdb.size()-ipath);
}

AccessDb::AccessDb(std::string partition, int reference) : 
  Partition_(partition), 
  ReferenceRun_(reference),
  FecVersionMajorId_(0), 
  FecVersionMinorId_(0), 
  FedVersionMajorId_(-1),
  FedVersionMinorId_(-1),
  ConnectionVersionMajorId_(0), 
  ConnectionVersionMinorId_(0),
  DcuInfoVersionMajorId_(0),
  DcuInfoVersionMinorId_(0),
  DcuPsuMapVersionMajorId_(0),
  DcuPsuMapVersionMinorId_(0),
  AnalysisVersionMapPointerId_(0),
  MaskVersionMajorId_(0),
  MaskVersionMinorId_(0),
  DownloadFec_(true),
  DownloadFed_(true),
  DownloadConnection_(true),
  DownloadPsu_(true),
  DownloadDetId_(true),
  CheckBadPll_(true),
  CheckMissingDevices_(true),
  CheckMissingFedChannels_(true),
  CheckFedConnection_(true),
  StripListName_("NONE"),
  DisableStrips_(false),
  UploadFeds_(false),
  UploadFecs_(false),
  UploadConnections_(false)
{

  string login="nil";
  string passwd ="nil";
  string path ="nil";

  getDbConfiguration(login,passwd,path);

  deviceFactory_ = new DeviceFactory(login,passwd,path);

  std::cout<<"Device Factory created" <<std::endl;

  if (ReferenceRun_ !=-1){
    std::cout<<"Reference run "<<ReferenceRun_<<std::endl;
    tkRunVector vRun_ = deviceFactory_->getAllRuns();    
    for (unsigned int j=0;j<vRun_.size();j++){
      TkRun* run = vRun_[j];
      if (run->getRunNumber()== (unsigned int) ReferenceRun_){
	std::cout<<"Downloading reference run informations "<<std::endl;
	run->display();
	FecVersionMajorId_ = run->getFecVersionMajorId ( );  
	FecVersionMinorId_ = run->getFecVersionMinorId ( ) ; 
	FedVersionMajorId_ = run->getFedVersionMajorId ( ) ; 
	FedVersionMinorId_ = run->getFedVersionMinorId ( ) ; 
	ConnectionVersionMajorId_ = run->getConnectionVersionMajorId ( ) ; 
	ConnectionVersionMinorId_ = run->getConnectionVersionMinorId ( ) ; 
	DcuInfoVersionMajorId_  = run->getDcuInfoVersionMajorId ( ) ; 
	DcuInfoVersionMinorId_  = run->getDcuInfoVersionMinorId ( ) ; 
	DcuPsuMapVersionMajorId_ = run->getDcuPsuMapVersionMajorId ( ) ; 
	DcuPsuMapVersionMinorId_ = run->getDcuPsuMapVersionMinorId ( ) ; 
	AnalysisVersionMapPointerId_ =run->getAnalysisVersionMapPointerId() ; 
	MaskVersionMajorId_  =run->getMaskVersionMajorId ( ) ; 
	MaskVersionMinorId_  =run->getMaskVersionMinorId ( ) ; 
      }
    }
  }
}

void AccessDb::Download(){

  if (DownloadPsu_){
    deviceFactory_->getDcuPsuMapPartition(Partition_);
    psus_= deviceFactory_->getPowerGroupDcuPsuMaps();
  }

  if (DownloadDetId_){
    deviceFactory_->addDetIdPartition(Partition_,DcuInfoVersionMajorId_,DcuInfoVersionMinorId_);
    std::cout<<"find Det ids : Partition " <<Partition_<<" DCUInfo: major "<<DcuInfoVersionMajorId_<<" minor "<<DcuInfoVersionMinorId_<<std::endl;
    mapInfo_=deviceFactory_->getInfos();
  }

  // Get FEC devices
  if (DownloadFec_){
    deviceFactory_->getFecDeviceDescriptions(Partition_,vDevice_,FecVersionMajorId_,FecVersionMinorId_,MaskVersionMajorId_,MaskVersionMinorId_,false,false) ;
    std::cout<<"find Device descriptions" <<FecVersionMajorId_<<"."<<FecVersionMinorId_<<" Mask "<<MaskVersionMajorId_<<"."<<MaskVersionMinorId_<<std::endl;
  }

  if (DownloadConnection_){
    deviceFactory_->getConnectionDescriptions(Partition_,vCon_,ConnectionVersionMajorId_,ConnectionVersionMinorId_,MaskVersionMajorId_,MaskVersionMinorId_) ;
    std::cout<<"find Connection descriptions" <<ConnectionVersionMajorId_<<"."<<ConnectionVersionMinorId_<<" Mask "<<MaskVersionMajorId_<<"."<<MaskVersionMinorId_<<std::endl;
  }

  if (DownloadFed_){
    fedVector_ = deviceFactory_->getFed9UDescriptions(Partition_,FedVersionMajorId_,FedVersionMinorId_,MaskVersionMajorId_,MaskVersionMinorId_) ;
    std::cout<<"find FED descriptions" <<FedVersionMajorId_<<"."<<FedVersionMinorId_<<" Mask "<<MaskVersionMajorId_<<"."<<MaskVersionMinorId_<<std::endl;    
  }
}


void AccessDb::DelayDetIds(std::string filename,std::string xmlname){

  std::ofstream out(xmlname.c_str());
  std::ifstream ti; 
  ti.open(filename.c_str());

  int detid=0, dcuid=0;
  float delay;
  vselected_.clear();

  while ( ti>>detid>>delay ){
    if (delay == 0){
      std::cout<<"Delay is 0 --> No need to delay the module with detid: "<<detid<<std::endl;
      continue;
    }
    dcuid = 0;
    // find the dcuid
    for(Sgi::hash_map<unsigned long, TkDcuInfo *>::iterator it=mapInfo_.begin();it!=mapInfo_.end();it++){
      if (detid==it->second->getDetId()){
	dcuid=it->second->getDcuHardId();
	break;
      }
    }
    if (dcuid==0){
      std::cout<<"Cannot find the dcuid for detid "<<detid<<" for the tracker partition "<<Partition_<<std::endl;
      continue;
    }
    DelayModule(dcuid,delay,out);
  }
  if (ti.eof()) {
    ti.close();
  } else {
    ti.close();
    std::string error = "Unable to read file " ;
    throw std::runtime_error(error);
  }
  return;
}


void AccessDb::DelayModule(unsigned int dcuid, float delay,std::ofstream & outfile){

  // Find the PLL
  bool goodPLL = false;

  for (unsigned int i=0;i<vCon_.size();i++){

    if (vCon_[i]->getDcuHardId() != dcuid ) continue;

    // Now Loop on Device and find the PLL
    for (unsigned int j=0; j < vDevice_.size();j++){
      
      deviceDescription *d = vDevice_[j] ;
      
      if (d->getDeviceType()!=PLL) continue;	
      if (d->getFecSlot()!=vCon_[i]->getFecSlot()) continue;
      if (d->getRingSlot()!=vCon_[i]->getRingSlot()) continue;
      if (d->getCcuAddress()!=vCon_[i]->getCcuAddress()) continue;
      if (d->getChannel()!=vCon_[i]->getI2cChannel()) continue;
      
      pllDescription* p = (pllDescription*) d;
      float plldelay = p->getDelayCoarse()*25.+p->getDelayFine()*25./24;
      float newdelay = plldelay+delay; // 
      int coarse = int(newdelay/25.);
      int fine = rint((newdelay-coarse*25.)*24./25.);
      if (fine == 24) { fine=0;coarse+=1;}
      //std::cout<<"### New Module PLL "<<std::hex<<dcuid<<std::dec<<p->getChannel()<<" Old Settings "<<plldelay<<" "<<(int) p->getDelayCoarse()<<"/"<<(int) p->getDelayFine();	      
      //std::cout<<" ### New Settings "<<newdelay<<" "<<coarse*25.+fine*25./24.<<" "<<coarse<<"/"<<fine<<std::endl;
      float newplldelay = coarse*25+fine*25./24.;
      outfile<<"<PLL dcuid=\""<<std::hex<<dcuid<<"\" ocoarse=\""<<std::dec<<(int) p->getDelayCoarse()<<"\" coarse=\""<<coarse<<"\" ofine=\""<<(int) p->getDelayFine()<<"\" fine=\""<<fine<<"\" odelay=\""<<plldelay<<"\" delay=\""<<newplldelay<<"\" difpll=\""<<std::fixed << std::setprecision(2)<<newplldelay-plldelay<<"\" required=\""<<std::fixed << std::setprecision(2)<< delay<<"\" >"<<std::endl;
      if (coarse<0 || coarse>15){
	std::cout<<" Invalid coarse delay "<<std::endl;
	exit(0);
      }
      if (fine<0 || fine>23){
	std::cout<<" Invalid fine delay "<<std::endl;
	exit(0);
      }      
      p->setDelayCoarse(coarse);
      p->setDelayFine(fine);
      goodPLL = true;
      break;
    }
    break;
  }
  
  // Repeat the loop for fed delay
  if(goodPLL){
    for (unsigned int i=0;i<vCon_.size();i++){
      if (vCon_[i]->getDcuHardId()!=dcuid ) continue;
      for (unsigned int k=0;k<fedVector_->size();k++){
	if ((*fedVector_)[k]->getFedId()!=vCon_[i]->getFedId()) continue;
	Fed9U::Fed9UAddress addr;
	addr.setFedChannel(static_cast<Fed9U::u8>(vCon_[i]->getFedChannel()));
	int cur_finedelay=(*fedVector_)[k]->getFineDelay(addr);
	int cur_coarsedelay=(*fedVector_)[k]->getCoarseDelay(addr);
	float feddelay=cur_coarsedelay*25.-cur_finedelay;
	float newdelay=feddelay-delay; // Subtract the delay on FED channel 
	//std::cout<<"\t### FED "<<vCon_[i]->getFedId()<<":"<<vCon_[i]->getFedChannel()<<" "<<cur_finedelay<<" "<<cur_coarsedelay<<" "<<feddelay<<" "<<newdelay<<std::endl;
       	int coarse = int(newdelay/25.)+1;
	int fine   = rint(coarse*25 -newdelay);
	if (fine ==25){ fine = 0; coarse-=1; }
	//std::cout<<" \t### "<<newdelay<<" "<<coarse*25.-fine<<" "<<coarse<<"/"<<fine<<std::endl;
	outfile<<"\t <FedChannel id=\""<<vCon_[i]->getFedId()<<"\" channel=\""<<vCon_[i]->getFedChannel()<<"\" ofine=\""<<cur_finedelay<<"\" fine=\""<<fine<<"\" ocoarse=\""<<cur_coarsedelay<<"\" coarse=\""<<coarse<<"\" odelay=\""<<feddelay<<"\" delay=\""<<coarse*25.-fine<<"\" dif=\""<<std::fixed << std::setprecision(2) <<coarse*25.-fine - feddelay<<"\" />"<< std::endl;
	
	(*fedVector_)[k]->setFineDelay(addr,fine);
	(*fedVector_)[k]->setCoarseDelay(addr,coarse);
	
      }     
      outfile<<"</PLL>"<<std::endl;
    }
  }
}


void  AccessDb::Check(){
  unsigned int nDcu =0;
  unsigned int nApv =0;
  unsigned int badpll=0;
  unsigned int apvError=0;
  for (unsigned int j=0;j<vDevice_.size();j++){
    deviceDescription      *d = vDevice_[j] ;
    if (d->getDeviceType()==DCU) {
      dcuDescription* dcu = (dcuDescription*) d;
      if (dcu->getDcuType() != DCUCCU) nDcu++;      
    }
    if (d->getDeviceType()==APV25) {
      nApv++;
      apvDescription* a = (apvDescription*) d;
      if (a->getApvError()>1) apvError++;
    }
    if (CheckBadPll_){
      if (d->getDeviceType()==PLL){
	pllDescription* p = (pllDescription*) d;
	if (p->getDelayCoarse()>=15)
	  {
	    badpll++;
	    std::cout<<" Bad PLL " <<dec<<
	      d->getFecSlot()<<" "<<
	      d->getRingSlot()<<" "<<
	      d->getCcuAddress()<<" "<<
	      d->getChannel()<<" ->"<<(int) p->getDelayCoarse()<<std::endl;
	  }
      }
    }
  }

  
  std::cout<<"Number of DETIDs "<<dec<<mapInfo_.size() <<std::endl;
  std::cout<<dec<<nDcu<<" Dcu Devices found "<<std::endl;
  std::cout<<dec<<nApv<<" Apv Devices found "<<std::endl;
  std::cout<<dec<<badpll<<" bad PLL devices "<<std::endl;
  std::cout<<dec<<apvError<<" Apv in Error "<<std::endl;
  std::cout<<dec<<vDevice_.size()<<" Devices found  with"<<vCon_.size()<<"connections and " <<fedVector_->size()<<" Feds" <<std::endl; 

  int missingDevices=0;
  int missingFedChannels=0;
  int nondisabledFedChannels=0;


  if (CheckMissingDevices_){
    for (unsigned int i=0;i<vCon_.size();i++){
      bool found=false;
      for (unsigned int j=0;j<vDevice_.size();j++){
	  deviceDescription      *d = vDevice_[j] ;
	  if (d->getDeviceType()!=APV25) continue;
	  if (d->getFecSlot()!=vCon_[i]->getFecSlot()) continue;
	  if (d->getRingSlot()!=vCon_[i]->getRingSlot()) continue;
	  if (d->getCcuAddress()!=vCon_[i]->getCcuAddress()) continue;
	  if (d->getChannel()!=vCon_[i]->getI2cChannel()) continue;
	  if (d->getAddress()!=vCon_[i]->getApvAddress()) continue;
	  found=true;
	  break;
      }
      if (!found){
	missingDevices++;
	std::cout<<dec<<vCon_[i]->getFecSlot()<<":"
		 <<vCon_[i]->getRingSlot()<<":"
		 <<vCon_[i]->getCcuAddress()<<":"
		 <<vCon_[i]->getI2cChannel()<<":"
		 <<vCon_[i]->getApvAddress()<<" is Missing "<< std::endl;
	for (unsigned int k=0;k<fedVector_->size();k++){
	  if ((*fedVector_)[k]->getFedId()!=vCon_[i]->getFedId()) continue;
	  Fed9U::Fed9UAddress addr;
	  addr.setFedChannel(static_cast<Fed9U::u8>(vCon_[i]->getFedChannel()));
	  addr.setChannelApv(0);         
	  std::cout <<dec<<"Status of FED  "<<vCon_[i]->getFedId()<<" Channel "<< vCon_[i]->getFedChannel()<<" is ";
	  std::cout<<(*fedVector_)[k]->getApvDisable(addr)<<std::endl;
	  printf("exec PkgMaskModules.MaskFedModules('%s',%d,%d,0);\n",Partition_.c_str(),vCon_[i]->getFedId(),vCon_[i]->getFedChannel());
	}
      } 
      
      // Now check if the channel is enable
      for (unsigned int k=0;k<fedVector_->size();k++){
	if ((*fedVector_)[k]->getFedId()!=vCon_[i]->getFedId()) continue;
	Fed9U::Fed9UAddress addr;
	addr.setFedChannel(static_cast<Fed9U::u8>(vCon_[i]->getFedChannel()));
	addr.setChannelApv(0);
	if ((*fedVector_)[k]->getApvDisable(addr)){
	  missingFedChannels++;
	  std::cout <<dec<<"Status of FED  "<<vCon_[i]->getFedId()<<" Channel "<< vCon_[i]->getFedChannel()<<" is ";
	  std::cout<<(*fedVector_)[k]->getApvDisable(addr)<<std::endl;
	}          
      }                
    }
  }

  // NOw loop on Fed channels and verify the connection exist
  bool disable=true;
  for (unsigned int k=0;k<fedVector_->size();k++){
    for (unsigned int l=0;l<96;l++){
      Fed9U::Fed9UAddress addr;
      addr.setFedChannel(static_cast<Fed9U::u8>(l));
      addr.setChannelApv(0);
      if (!(*fedVector_)[k]->getApvDisable(addr)){
	bool found=false;
	for (unsigned int i=0;i<vCon_.size();i++){
	  if ((*fedVector_)[k]->getFedId()!=vCon_[i]->getFedId()) continue;
	  if (vCon_[i]->getFedChannel()!=l) continue;
	  found=true;
	  break;
	}
	if (!found){
	  std::cout<<"Fed "<<(*fedVector_)[k]->getFedId()<<" Channel "<<l<<" is enabled but not connected"<<std::endl;	  
	  if (disable){ 
	    (*fedVector_)[k]->setApvDisable(addr,true);
	    Fed9U::Fed9UAddress addr1;
	    addr1.setFedChannel(static_cast<Fed9U::u8>(l));
	    addr1.setChannelApv(1);
	    (*fedVector_)[k]->setApvDisable(addr1,true);	    
	  }
	  nondisabledFedChannels++;
	}
      }
    }
  }

  int nstillenable=0;
  for (unsigned int k=0;k<fedVector_->size();k++){
    for (unsigned int l=0;l<96;l++){
      Fed9U::Fed9UAddress addr;
      addr.setFedChannel(static_cast<Fed9U::u8>(l));
      addr.setChannelApv(0);
      if (!(*fedVector_)[k]->getApvDisable(addr)){
	bool found=false;
	for (unsigned int i=0;i<vCon_.size();i++){
	  if ((*fedVector_)[k]->getFedId()!=vCon_[i]->getFedId()) continue;
	  if (vCon_[i]->getFedChannel()!=l) continue;
	  found=true;
	  break;
	}
	if (!found){
	  std::cout<<"Fed "<<(*fedVector_)[k]->getFedId()<<" Channel "<<l<<" is enabled but not connected"<<std::endl;	  
	  nstillenable++;
	}
      }
    }
  }
  
  std::cout<<"Connected Missing Devices"<< missingDevices<<std::endl;
  std::cout<<"Connected missing FED channels"<< missingFedChannels<<std::endl;
  std::cout<<"Non disable unconnected FED channels"<< nondisabledFedChannels<<std::endl;
  std::cout<<"Non disable unconnected FED channels after disbale"<< nstillenable<<std::endl;
  
}

void AccessDb::Upload()
{
  if (getUploadFeds()){
    unsigned short versionMajor,versionMinor;
    deviceFactory_->setFed9UDescriptions((*fedVector_),Partition_,&versionMajor,&versionMinor,1);
    std::cout << "fed9U uploaded to new vesion " << " "<< versionMajor <<"."<<versionMinor << std::endl;
  }
  if (getUploadFecs()){
    unsigned int major,minor;
    bool ismajor=false;
    deviceFactory_->setFecDeviceDescriptions (vDevice_, Partition_, &major, &minor, ismajor, false) ;      
    std::cout << "fec devices uploaded to new vesion " << " "<< major <<"."<<minor << std::endl;
  } // Upload the devices 
}

// Setters
void AccessDb::setReferenceRun(  int t){ ReferenceRun_=t;}
void AccessDb::setFecVersionMajorId(  unsigned int t){ FecVersionMajorId_=t;}
void AccessDb::setFecVersionMinorId(  unsigned int t){ FecVersionMinorId_=t;}
void AccessDb::setFedVersionMajorId(  int t){ FedVersionMajorId_=t;}
void AccessDb::setFedVersionMinorId(  int t){ FedVersionMinorId_=t;}
void AccessDb::setConnectionVersionMajorId(  unsigned int t){ ConnectionVersionMajorId_=t;} 
void AccessDb::setConnectionVersionMinorId(  unsigned int t){ ConnectionVersionMinorId_=t;}
void AccessDb::setDcuInfoVersionMajorId(  unsigned int t){ DcuInfoVersionMajorId_=t;}
void AccessDb::setDcuInfoVersionMinorId(  unsigned int t){ DcuInfoVersionMinorId_=t;}
void AccessDb::setDcuPsuMapVersionMajorId(  unsigned int t){ DcuPsuMapVersionMajorId_=t;}
void AccessDb::setDcuPsuMapVersionMinorId(  unsigned int t){ DcuPsuMapVersionMinorId_=t;}
void AccessDb::setAnalysisVersionMapPointerId(  unsigned int t){ AnalysisVersionMapPointerId_=t;}
void AccessDb::setMaskVersionMajorId(  unsigned int t){ MaskVersionMajorId_=t;}
void AccessDb::setMaskVersionMinorId(  unsigned int t){ MaskVersionMinorId_=t;}
void AccessDb::setDownloadFec(  bool t){ DownloadFec_=t;}
void AccessDb::setDownloadFed(  bool t){ DownloadFed_=t;}
void AccessDb::setDownloadConnection(  bool t){ DownloadConnection_=t;}
void AccessDb::setDownloadPsu(  bool t){ DownloadPsu_=t;}
void AccessDb::setDownloadDetId(  bool t ){ DownloadDetId_=t;}
void AccessDb::setCheckBadPll(  bool t){ CheckBadPll_=t;}
void AccessDb::setCheckMissingDevices(  bool t){ CheckMissingDevices_=t;}
void AccessDb::setCheckMissingFedChannels(  bool t){ CheckMissingFedChannels_=t;}
void AccessDb::setCheckFedConnection(  bool t){ CheckFedConnection_=t;}
void AccessDb::setStripListName(   std::string t){ StripListName_=t;}
void AccessDb::setDisableStrips(  bool t){ DisableStrips_=t;}
void AccessDb::setUploadFeds(  bool t){ UploadFeds_=t;}
void AccessDb::setUploadFecs(  bool t){ UploadFecs_=t;}
void AccessDb::setUploadConnections(  bool t){ UploadConnections_=t;}
void AccessDb::setPartition(  std::string t ){ Partition_=t;}


// Getters
int AccessDb::getReferenceRun(){ return ReferenceRun_;}
unsigned int AccessDb::getFecVersionMajorId(){ return FecVersionMajorId_;} 
unsigned int AccessDb::getFecVersionMinorId(){ return FecVersionMinorId_;} 
int AccessDb::getFedVersionMajorId(){ return FedVersionMajorId_;} 
int AccessDb::getFedVersionMinorId(){ return FedVersionMinorId_;} 
unsigned int AccessDb::getConnectionVersionMajorId(){ return ConnectionVersionMajorId_;} 
unsigned int AccessDb::getConnectionVersionMinorId(){ return ConnectionVersionMinorId_;} 
unsigned int AccessDb::getDcuInfoVersionMajorId(){ return DcuInfoVersionMajorId_;} 
unsigned int AccessDb::getDcuInfoVersionMinorId(){ return DcuInfoVersionMinorId_;} 
unsigned int AccessDb::getDcuPsuMapVersionMajorId(){ return DcuPsuMapVersionMajorId_;} 
unsigned int AccessDb::getDcuPsuMapVersionMinorId(){ return DcuPsuMapVersionMinorId_;} 
unsigned int AccessDb::getAnalysisVersionMapPointerId(){ return AnalysisVersionMapPointerId_;} 
unsigned int AccessDb::getMaskVersionMajorId(){ return MaskVersionMajorId_;} 
unsigned int AccessDb::getMaskVersionMinorId(){ return MaskVersionMinorId_;}
bool AccessDb::getDownloadFec(){ return DownloadFec_;}
bool AccessDb::getDownloadFed(){ return DownloadFed_;}
bool AccessDb::getDownloadConnection(){ return DownloadConnection_;}
bool AccessDb::getDownloadPsu(){ return DownloadPsu_;}
bool AccessDb::getDownloadDetId(){ return DownloadDetId_;}
bool AccessDb::getCheckBadPll(){ return CheckBadPll_;}
bool AccessDb::getCheckMissingDevices(){ return CheckMissingDevices_;}
bool AccessDb::getCheckMissingFedChannels(){ return CheckMissingFedChannels_;}
bool AccessDb::getCheckFedConnection(){ return CheckFedConnection_;}
std::string AccessDb::getStripListName(){ return StripListName_;}
bool AccessDb::getDisableStrips(){ return DisableStrips_;}
bool AccessDb::getUploadFeds(){ return UploadFeds_;}
bool AccessDb::getUploadFecs(){ return UploadFecs_;}
bool AccessDb::getUploadConnections(){ return UploadConnections_;}
std::string AccessDb::getPartition(){ return Partition_;}


