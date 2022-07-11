#include "TrackerDAQAnalysis/RandomDelayScanAnalysis/interface/countEvents.h"

//! ctor
countEvents::countEvents(const edm::ParameterSet& iConfig) {
}

// ----------------------------------------------------------------

//! dtor
countEvents::~countEvents()
{}

// ----------------------------------------------------------------


//! loop over the reco particles and count leptons
bool countEvents::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout<<"countEvents::filter,"<<iEvent.eventAuxiliary().run()<<","<<iEvent.eventAuxiliary().luminosityBlock()<<","<<iEvent.eventAuxiliary().event()<<std::endl;
  return true;
}

DEFINE_FWK_MODULE(countEvents);
