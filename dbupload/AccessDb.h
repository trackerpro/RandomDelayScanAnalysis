#ifndef _AccessDb_
#define _AccessDb_
#include "DeviceFactory.h"
#include "ConnectionFactory.h"
#include "TkDcuPsuMapFactory.h"
#include "PedestalsAnalysisDescription.h"
#include "TimingAnalysisDescription.h"
#include <string>
#include <iostream>
#include <fstream>

class AccessDb
{

 public:

  AccessDb(std::string partition,int nrun=0);
  void Download();
  void Check();
  void Upload();

  
  // Setters
  void setReferenceRun(  int t);
  void setFecVersionMajorId( unsigned int t); 
  void setFecVersionMinorId( unsigned int t); 
  void setFedVersionMajorId( int t); 
  void setFedVersionMinorId( int t); 
  void setConnectionVersionMajorId( unsigned int t); 
  void setConnectionVersionMinorId( unsigned int t); 
  void setDcuInfoVersionMajorId( unsigned int t); 
  void setDcuInfoVersionMinorId( unsigned int t); 
  void setDcuPsuMapVersionMajorId( unsigned int t); 
  void setDcuPsuMapVersionMinorId( unsigned int t); 
  void setAnalysisVersionMapPointerId( unsigned int t); 
  void setMaskVersionMajorId( unsigned int t); 
  void setMaskVersionMinorId( unsigned int t);
  void setDownloadFec( bool t);
  void setDownloadFed( bool t);
  void setDownloadConnection( bool t);
  void setDownloadPsu( bool t);
  void setDownloadDetId( bool t);
  void setCheckBadPll( bool t);
  void setCheckMissingDevices( bool t);
  void setCheckMissingFedChannels( bool t);
  void setCheckFedConnection( bool t);
  void setStripListName(  std::string t);
  void setDisableStrips( bool t);
  void setUploadFeds( bool t);
  void setUploadFecs( bool t);
  void setUploadConnections( bool t);
  void setPartition( std::string t);
 // Getters
  int getReferenceRun();
  unsigned int getFecVersionMajorId(); 
  unsigned int getFecVersionMinorId(); 
  int getFedVersionMajorId(); 
  int getFedVersionMinorId(); 
  unsigned int getConnectionVersionMajorId(); 
  unsigned int getConnectionVersionMinorId(); 
  unsigned int getDcuInfoVersionMajorId(); 
  unsigned int getDcuInfoVersionMinorId(); 
  unsigned int getDcuPsuMapVersionMajorId(); 
  unsigned int getDcuPsuMapVersionMinorId(); 
  unsigned int getAnalysisVersionMapPointerId(); 
  unsigned int getMaskVersionMajorId(); 
  unsigned int getMaskVersionMinorId();
  bool getDownloadFec();
  bool getDownloadFed();
  bool getDownloadConnection();
  bool getDownloadPsu();
  bool getDownloadDetId();
  bool getCheckBadPll();
  bool getCheckMissingDevices();
  bool getCheckMissingFedChannels();
  bool getCheckFedConnection();
  std::string getStripListName();
  bool getDisableStrips();
  bool getUploadFeds();
  bool getUploadFecs();
  bool getUploadConnections();
  std::string getPartition();
  
  tkRunVector& getRunVector(){ return vRun_;}
  tkDcuPsuMapVector getPsuVector() { return psus_;}
  Sgi::hash_map<unsigned long, TkDcuInfo *>& getMapInfo() {return mapInfo_;}
  deviceVector& getDeviceVector(){return  vDevice_;}
  ConnectionVector& getConnectionVector(){ return vCon_;}
  std::vector<Fed9U::Fed9UDescription*>* getFedVector() {return fedVector_ ;}


  void DelayModule(unsigned int dcuid, float delay,std::ofstream &o);
  void DelayDetIds(std::string filename,std::string xmlname="delay_Modules.xml");
 
 private:
  void PrintConnection(unsigned int i);
  DeviceFactory* deviceFactory_;
  std::string Partition_;
  int ReferenceRun_;
  unsigned int FecVersionMajorId_; 
  unsigned int FecVersionMinorId_; 
  int FedVersionMajorId_; 
  int FedVersionMinorId_; 
  unsigned int ConnectionVersionMajorId_; 
  unsigned int ConnectionVersionMinorId_; 
  unsigned int DcuInfoVersionMajorId_; 
  unsigned int DcuInfoVersionMinorId_; 
  unsigned int DcuPsuMapVersionMajorId_; 
  unsigned int DcuPsuMapVersionMinorId_; 
  unsigned int AnalysisVersionMapPointerId_; 
  unsigned int MaskVersionMajorId_; 
  unsigned int MaskVersionMinorId_;
  bool DownloadFec_;
  bool DownloadFed_;
  bool DownloadConnection_;
  bool DownloadPsu_;
  bool DownloadDetId_;
  bool CheckBadPll_;
  bool CheckMissingDevices_;
  bool CheckMissingFedChannels_;
  bool CheckFedConnection_;
  std::string StripListName_;
  bool DisableStrips_;
  bool UploadFeds_;
  bool UploadFecs_;
  bool UploadConnections_;
  bool DumpPsu_;
  // Vectors
  tkRunVector vRun_;
  tkDcuPsuMapVector psus_;
  Sgi::hash_map<unsigned long, TkDcuInfo *> mapInfo_;
  deviceVector vDevice_;
  ConnectionVector vCon_;
  std::vector<Fed9U::Fed9UDescription*> *fedVector_ ;

  std::vector<unsigned int> vselected_;

};

#endif
