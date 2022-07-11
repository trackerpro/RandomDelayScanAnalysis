#ifndef countEvents_h
#define countEvents_h

#include <iostream>
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class countEvents : public edm::stream::EDFilter<> {
 
 public:
 
  //! ctor
  explicit countEvents (const edm::ParameterSet&);
 
  //! dtor 
  ~countEvents();
 
 private:
 
  //! the actual filter method 
  bool filter(edm::Event&, const edm::EventSetup&);
 
 private:
};

#endif
