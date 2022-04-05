#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "CalibFormats/SiStripObjects/interface/SiStripFecCabling.h"
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"
#include "CondFormats/DataRecord/interface/SiStripFedCablingRcd.h"
#include "DataFormats/SiStripCommon/interface/SiStripEnumsAndStrings.h"
#include "DataFormats/SiStripCommon/interface/SiStripEventSummary.h"
#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include "DataFormats/SiStripCommon/interface/SiStripFedKey.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <ctime>
#include <iomanip>
#include <memory>
#include <sstream>

#include <TTree.h>

using namespace sistrip; 

class CheckStripEventSummary : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {

public:
  explicit CheckStripEventSummary(const edm::ParameterSet&);
  ~CheckStripEventSummary();

  virtual void beginJob() override;
  virtual void beginRun(edm::Run const&, const edm::EventSetup&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const edm::EDGetTokenT<SiStripEventSummary> stripEventSummaryToken_;
  const edm::ESGetToken<SiStripFedCabling, SiStripFedCablingRcd> fedCablingToken_;
  SiStripFedCabling* fedCabling_;
  SiStripFecCabling* fecCabling_;

  TTree* tree;
  unsigned event;
  unsigned runnumber;
  unsigned lumisection;
  unsigned runtype;
  unsigned dcuid;
  unsigned fedreadout;
  unsigned pllcoarse;
  unsigned pllfine;
  unsigned delayrange;
  unsigned delaystepsize;
  unsigned delaystep;
  
};

CheckStripEventSummary::CheckStripEventSummary(const edm::ParameterSet& pset):
  stripEventSummaryToken_(consumes<SiStripEventSummary>(pset.getParameter<edm::InputTag>("stripEventSummary"))),
  fedCablingToken_(esConsumes<edm::Transition::BeginRun>()),
  fecCabling_(nullptr){
  usesResource("TFileService"); 
}


CheckStripEventSummary::~CheckStripEventSummary(){
  if (fedCabling_) {
    fedCabling_ = nullptr;
  }
  if (fecCabling_) {
    delete fecCabling_;
    fecCabling_ = nullptr;
  }
}

void CheckStripEventSummary::beginJob() {
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree","tree");
  tree->Branch("event", &event, "event/i");
  tree->Branch("runnumber", &runnumber, "runnumber/i");
  tree->Branch("lumisection", &lumisection, "lumisection/i");
  tree->Branch("runtype", &runtype, "runtype/i");
  tree->Branch("fedreadout", &fedreadout, "fedreadout/i");
  tree->Branch("dcuid", &dcuid, "dcuid/i");
  tree->Branch("pllcoarse", &pllcoarse, "pllcoarse/i");
  tree->Branch("pllfine", &pllfine, "pllfine/i");
  tree->Branch("delayrange", &delayrange, "delayrange/i");
  tree->Branch("delaystepsize", &delaystepsize, "delaystepsize/i");
}

void CheckStripEventSummary::beginRun(edm::Run const& run, const edm::EventSetup& setup) {
  const auto& fed_cabling = setup.getData(fedCablingToken_);
  fedCabling_ = const_cast<SiStripFedCabling*>(&fed_cabling);
  fecCabling_ = new SiStripFecCabling(fed_cabling);
}


void CheckStripEventSummary::analyze(const edm::Event& ievent, const edm::EventSetup& setup) {

  edm::Handle<SiStripEventSummary> summary;
  ievent.getByToken(stripEventSummaryToken_, summary);

  runnumber = ievent.id().run();
  lumisection = ievent.luminosityBlock();
  event = ievent.id().event();
  runtype = summary->runType();
  fedreadout = summary->fedReadoutMode();
  dcuid = summary->dcuId();
  pllcoarse = summary->pllCoarse();
  pllfine = summary->pllFine();
  delayrange = summary->delayRange();
  delaystepsize = summary->delayStepSize();

  tree->Fill();
}

void CheckStripEventSummary::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void CheckStripEventSummary::endRun(edm::Run const&, edm::EventSetup const&) {
}



DEFINE_FWK_MODULE(CheckStripEventSummary);
