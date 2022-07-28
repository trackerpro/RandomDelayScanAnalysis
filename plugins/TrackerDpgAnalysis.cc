// -*- C++ -*-
//
// Package:    DPGAnalysis
// Class:      TrackerDpgAnalysis
//
/**\class TrackerDpgAnalysis TrackerDpgAnalysis.cc DPGAnalysis/SiStripTools/plugins/TrackerDpgAnalysis.cc
 Description: analysis of the clusters and digis in the tracker

 Implementation:
      analysis of the clusters and digis in the tracker
*/
//
// Original Author:  Christophe DELAERE
//         Created:  Tue Sep 23 02:11:44 CEST 2008
//         Revised:  Thu Nov 26 10:00:00 CEST 2009
// part of the code was inspired by http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/YGao/LhcTrackAnalyzer/
// part of the code was inspired by
// other inputs from Andrea Giammanco, Gaelle Boudoul, Andrea Venturi, Steven Lowette, Gavril Giurgiu
// $Id: TrackerDpgAnalysis.cc,v 1.14 2013/02/27 19:49:47 wmtan Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <algorithm>
#include <functional>
#include <string>
#include <sstream>
#include <fstream>

// root include files
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/transform.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/TrackerMap/interface/TrackerMap.h"
#include "CondFormats/SiStripObjects/interface/FedChannelConnection.h"
#include "CondFormats/DataRecord/interface/SiStripFedCablingRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiStripCommon/interface/SiStripEventSummary.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerMultiRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit2DLocalPos.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit1D.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"
#include "DPGAnalysis/SiStripTools/interface/EventShape.h"

//
// class decleration
//
typedef math::XYZPoint Point;

class TrackerDpgAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns>  {
public:
  explicit TrackerDpgAnalysis(const edm::ParameterSet&);
  ~TrackerDpgAnalysis() override;
  
protected:

  // vertex function
  std::map<size_t,int> inVertex(const reco::VertexCollection&, unsigned int);
  double sumPtSquared(const reco::Vertex&);
  
  // strip cluster ontrack
  std::vector<int> onTrack(edm::Handle<edmNew::DetSetVector<SiStripCluster> >&,const reco::TrackCollection&, unsigned int );
  std::vector<double> onTrackAngles(edm::Handle<edmNew::DetSetVector<SiStripCluster> >&,const std::vector<Trajectory>& );

  void insertMeasurement(std::multimap<const unsigned int,std::pair<LocalPoint,double> >&, const TransientTrackingRecHit*,double);
  void insertMeasurement(std::multimap<const unsigned int,std::pair<int,int> >&, const TrackingRecHit*,int);
  void insertMeasurement(std::multimap<const unsigned int,std::pair<LocalPoint,std::pair<double,double> > >&, const TransientTrackingRecHit*,double,double);
  void insertMeasurement(std::multimap<const unsigned int,std::pair<std::pair<float, float>,int> >&, const TrackingRecHit*,int);

  // for extracting the delay information
  std::string toStringName(unsigned int, const TrackerTopology*);
  std::string toStringId(unsigned int);
  std::map<unsigned int,float> delay(const std::string &);
  
private:
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void beginJob() override ;  
  virtual void endJob() override ;  
  
  SiStripClusterInfo siStripClusterInfo_;
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster> > clusterToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > dedx1Token_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > dedx2Token_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > dedx3Token_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> L1Token_;
  edm::EDGetTokenT<SiStripEventSummary> stripEventSummaryToken_;  
  edm::InputTag HLTTag_;
  edm::EDGetTokenT<edm::TriggerResults> HLTToken_;
  edm::EDGetTokenT<reco::TrackCollection>  trackTokens_;
  edm::EDGetTokenT<std::vector<Trajectory> >  trajectoryTokens_;
  edm::EDGetTokenT<TrajTrackAssociationCollection> trajTrackAssoTokens_;
  const SiStripFedCabling* cabling_;
  const TrackerGeometry* tracker_;
  std::multimap<const unsigned int,const FedChannelConnection*> connections_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoToken_;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  edm::ESGetToken<SiStripFedCabling, SiStripFedCablingRcd> fedCablingToken_;

  bool functionality_offtrackClusters_;
  bool functionality_ontrackClusters_;
  bool functionality_tracks_;
  bool functionality_vertices_;
  bool functionality_events_;

  edm::ParameterSet pset_; 
  std::vector<std::string>  hltNames_;  // name of each HLT algorithm
  std::string delayFileName_; // name of file with PLL delya map
  std::map<std::string,int> triggerPathsMap;

  // Output trees and their branches
  TTree* clusters_;
  TTree* pixclusters_;
  TTree* tracks_;
  TTree* vertices_;
  TTree* event_;
  TTree* readoutmap_;

  bool onTrack_, isValid_pvtx_, isFake_pvtx_;
  std::string moduleName_, moduleId_, PSUname_;
  std::vector<bool> hltTriggerBits_;
  std::vector<std::string> hltTriggerNames_;
  unsigned int runtype_, fedreadout_, delayrange_, delaystepsize_, delaystep_;
  unsigned int fecCrate_, fecSlot_, fecRing_, ccuAdd_, ccuChan_, lldChannel_, fedId_, fedCh_, fiberLength_;
  unsigned int vertexid_, eventid_, runid_, globalvertexid_;
  unsigned int detid_, dcuId_, type_;
  unsigned int nclusters_, nclustersOntrack_, dedxNoM_, quality_, foundhits_, lostHits_, ndof_;
  unsigned int nVertices_, nLayers_,foundhitsStrips_,foundhitsPixels_,losthitsStrips_,losthitsPixels_;
  unsigned int nTracks_pvtx_;
  unsigned int clSize_, clSizeX_, clSizeY_;
  unsigned int orbit_, orbitL1_, bx_, store_, time_;
  unsigned int lumisec_, physicsDeclared_;
  unsigned int globaltrackid_, trackid_, ntracks_, ntrajs_;
  float globalX_, globalY_, globalZ_;
  float measX_, measY_, errorX_, errorY_;
  float angle_, maxCharge_, maxChargeCorrected_;
  float clCorrectedCharge_, clCorrectedSignalOverNoise_;
  float clNormalizedCharge_, clNormalizedNoise_, clSignalOverNoise_;
  float clBareNoise_, clBareCharge_;
  float clWidth_, clPosition_, thickness_, stripLength_, distance_;
  float eta_, phi_, chi2_;
  float dedx1_, dedx2_, dedx3_;
  float fBz_, clPositionX_, clPositionY_, alpha_, beta_, chargeCorr_;
  float recx_pvtx_, recy_pvtx_, recz_pvtx_, recx_err_pvtx_, recy_err_pvtx_, recz_err_pvtx_, sumptsq_pvtx_;
  float pterr_, etaerr_, phierr_;
  float dz_, dzerr_, dzCorr_, dxy_, dxyerr_, dxyCorr_;
  float qoverp_, xPCA_, yPCA_, zPCA_, trkWeightpvtx_;
  float charge_, p_, pt_;
  float bsX0_, bsY0_, bsZ0_, bsSigmaZ_, bsDxdz_, bsDydz_;
  float delayrandom_, delay_;
  float lowPixelProbabilityFraction_;
  std::map<unsigned int, float> delayDetIdMap_;
};

//
// constructors and destructor
//
TrackerDpgAnalysis::TrackerDpgAnalysis(const edm::ParameterSet& iConfig):
  siStripClusterInfo_(consumesCollector(), std::string("")),
  magFieldToken_(esConsumes()),
  tTopoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
  tkGeomToken_(esConsumes<edm::Transition::BeginRun>()),
  fedCablingToken_(esConsumes<edm::Transition::BeginRun>()){

  // members
  usesResource();
  usesResource("TFileService");
  // store PSET
  pset_ = iConfig;

  // enable/disable functionalities
  functionality_offtrackClusters_ = iConfig.getUntrackedParameter<bool>("keepOfftrackClusters",true);
  functionality_ontrackClusters_  = iConfig.getUntrackedParameter<bool>("keepOntrackClusters",true);
  functionality_tracks_           = iConfig.getUntrackedParameter<bool>("keepTracks",true);
  functionality_vertices_         = iConfig.getUntrackedParameter<bool>("keepVertices",true);
  functionality_events_           = iConfig.getUntrackedParameter<bool>("keepEvents",true);
  
  // parameters
  clusterToken_    = consumes<edmNew::DetSetVector<SiStripCluster> >(iConfig.getParameter<edm::InputTag>("ClustersLabel"));
  trackTokens_     = consumes<reco::TrackCollection> (iConfig.getParameter<edm::InputTag>("TracksLabel"));
  trajectoryTokens_  = consumes<std::vector<Trajectory> > (iConfig.getParameter<edm::InputTag>("TracksLabel"));
  trajTrackAssoTokens_  = consumes<TrajTrackAssociationCollection> (iConfig.getParameter<edm::InputTag>("TracksLabel"));
  dedx1Token_      = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("DeDx1Label"));
  dedx2Token_      = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("DeDx2Label"));
  dedx3Token_      = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("DeDx3Label"));
  vertexToken_     = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexLabel"));
  bsToken_         = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotLabel"));
  L1Token_         = consumes<L1GlobalTriggerReadoutRecord>(iConfig.getParameter<edm::InputTag>("L1Label"));
  stripEventSummaryToken_ = consumes<SiStripEventSummary>(iConfig.getParameter<edm::InputTag>("StripEventSummary"));
  HLTTag_          = iConfig.getParameter<edm::InputTag>("HLTLabel");
  HLTToken_        = consumes<edm::TriggerResults>(HLTTag_);
  hltNames_        = iConfig.getParameter<std::vector<std::string> >("HLTNames");
  delayFileName_   = iConfig.getParameter<std::string>("DelayFileName");

}

TrackerDpgAnalysis::~TrackerDpgAnalysis()
{}


// ------------ method called to for each event  ------------
void
TrackerDpgAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace reco;
  using namespace std;
  using reco::TrackCollection;
  
  // Sistrip summary
  edm::Handle<SiStripEventSummary> summary;
  iEvent.getByToken(stripEventSummaryToken_, summary);
  
  // Trigger info
  edm::Handle<L1GlobalTriggerReadoutRecord> gtrr_handle;
  iEvent.getByToken(L1Token_, gtrr_handle);
  const L1GlobalTriggerReadoutRecord* gtrr = gtrr_handle.product();
  L1GtFdlWord fdlWord = gtrr->gtFdlWord();
  
  edm::Handle<edm::TriggerResults> trh;
  iEvent.getByToken(HLTToken_, trh);
  
  // load beamspot
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken_,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;
  
  // load primary vertex
  edm::Handle<reco::VertexCollection> vertexCollectionHandle;
  iEvent.getByToken(vertexToken_,vertexCollectionHandle);
  const reco::VertexCollection vertexColl = *(vertexCollectionHandle.product());
  
  // load dedx info
  Handle<ValueMap<DeDxData> > dEdx1Handle;
  Handle<ValueMap<DeDxData> > dEdx2Handle;
  Handle<ValueMap<DeDxData> > dEdx3Handle;
  try {iEvent.getByToken(dedx1Token_, dEdx1Handle);} catch ( cms::Exception& ) {;}
  try {iEvent.getByToken(dedx2Token_, dEdx2Handle);} catch ( cms::Exception& ) {;}
  try {iEvent.getByToken(dedx3Token_, dEdx3Handle);} catch ( cms::Exception& ) {;}
  const ValueMap<DeDxData> dEdxTrack1 = *dEdx1Handle.product();
  const ValueMap<DeDxData> dEdxTrack2 = *dEdx2Handle.product();
  const ValueMap<DeDxData> dEdxTrack3 = *dEdx3Handle.product();

  // load track collections
  reco::TrackCollection trackCollection;
  edm::Handle<reco::TrackCollection> trackCollectionHandle;
  iEvent.getByToken(trackTokens_,trackCollectionHandle);
  trackCollection = *trackCollectionHandle.product();
  ntracks_ = trackCollection.size();
  
  // load the trajectory collections
  std::vector<Trajectory> trajectoryCollection;
  edm::Handle<std::vector<Trajectory> > trajectoryCollectionHandle;
  iEvent.getByToken(trajectoryTokens_,trajectoryCollectionHandle);
  trajectoryCollection = *trajectoryCollectionHandle.product();
  ntrajs_ = trajectoryCollection.size();
  
  // load the tracks/traj association maps
  TrajTrackAssociationCollection TrajToTrackMap;
  Handle<TrajTrackAssociationCollection> trajTrackAssociationHandle;
  iEvent.getByToken(trajTrackAssoTokens_,trajTrackAssociationHandle);
  TrajToTrackMap = *trajTrackAssociationHandle.product();
  
   // load the clusters
  edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusters;
  iEvent.getByToken(clusterToken_,clusters);
  
  // sanity check
  if(!(!trackCollection.empty() && !trajectoryCollection.empty())){
    edm::LogError("BadConfig") << " the track and trajectory collections are empty --> please check ";
    return;
  }
  
  // SiStrip cluster
  siStripClusterInfo_.initEvent(iSetup);
  
  eventid_     = iEvent.id().event();
  runid_       = iEvent.id().run();
  lumisec_     = iEvent.eventAuxiliary().luminosityBlock();
  
  // determine if each cluster is on a track or not, and record the trackid
  std::vector<int> stripClusterOntrackIndices = onTrack(clusters,trackCollection,globaltrackid_+1);
  nclustersOntrack_ = std::count_if(stripClusterOntrackIndices.begin(),stripClusterOntrackIndices.end(),bind2nd(not_equal_to<int>(), -1));
  
  // load event info
  bx_          = iEvent.eventAuxiliary().bunchCrossing();
  orbit_       = iEvent.eventAuxiliary().orbitNumber();
  store_       = iEvent.eventAuxiliary().storeNumber();
  time_        = iEvent.eventAuxiliary().time().value();
  fBz_         = fabs(iSetup.getData(magFieldToken_).inTesla(GlobalPoint(0, 0, 0)).z());
  
  if(recoBeamSpotHandle.isValid()) {
    bsX0_ = bs.x0();
    bsY0_ = bs.y0();
    bsZ0_ = bs.z0();
    bsSigmaZ_ = bs.sigmaZ();
    bsDxdz_ = bs.dxdz();
    bsDydz_ = bs.dydz();
  } else {
    bsX0_ = 0.;
    bsY0_ = 0.;
    bsZ0_ = 0.;
    bsSigmaZ_ = 0.;
    bsDxdz_ = 0.;
    bsDydz_ = 0.;
  }
  
  runtype_     = static_cast<unsigned int>(summary->runType());
  fedreadout_  = summary->fedReadoutMode();
  delayrange_  = summary->randomDelayRange();
  delaystepsize_ = summary->randomDelayStepSize();
  delaystep_   = summary->randomDelayStep();
  
  orbitL1_         = fdlWord.orbitNr();
  physicsDeclared_ = fdlWord.physicsDeclared();     
  hltTriggerBits_.clear();
  hltTriggerNames_.clear();
  for(auto key : triggerPathsMap){
    hltTriggerNames_.push_back(key.first);
    hltTriggerBits_.push_back(trh->accept(key.second));
  }
   
  nVertices_ = 0;
  for(reco::VertexCollection::const_iterator v=vertexColl.begin(); v!=vertexColl.end(); ++v) {
    if(v->isValid() && !v->isFake()) 
      ++nVertices_;
  }
  
  // iterate over vertices
  std::map<size_t,int> trackVertices = inVertex(vertexColl, globalvertexid_+1);
  for(reco::VertexCollection::const_iterator v=vertexColl.begin(); v!=vertexColl.end(); ++v) {
    nTracks_pvtx_ = v->tracksSize();
    sumptsq_pvtx_ = sumPtSquared(*v);
    isValid_pvtx_ = int(v->isValid());
    isFake_pvtx_ =  int(v->isFake());
    recx_pvtx_ = v->x();
    recy_pvtx_ = v->y();
    recz_pvtx_ = v->z();
    recx_err_pvtx_ = v->xError();
    recy_err_pvtx_ = v->yError();
    recz_err_pvtx_ = v->zError();
    globalvertexid_++;
    if(functionality_vertices_)
      vertices_->Fill();
  }

  // iterate over tracks
  const Point beamSpot = recoBeamSpotHandle.isValid() ? Point(recoBeamSpotHandle->x0(), recoBeamSpotHandle->y0(), recoBeamSpotHandle->z0()) : Point(0, 0, 0);
  
  // Save track info
  lowPixelProbabilityFraction_ = 0;
  unsigned int n_hits_barrel = 0;
  unsigned int n_hits_lowprob = 0;
  for(TrajTrackAssociationCollection::const_iterator it = TrajToTrackMap.begin(); it!=TrajToTrackMap.end(); ++it) {       
    reco::TrackRef itTrack  = it->val;
    edm::Ref<std::vector<Trajectory> > traj  = it->key; // bug to find type of the key
    //associate with PV
    std::map<size_t,int>::const_iterator theV = trackVertices.find(itTrack.key());
    vertexid_ = (theV!=trackVertices.end()) ? theV->second : 0;
    globaltrackid_++;
    try { // not all track collections have the dedx info... indeed at best one.
      dedxNoM_ = dEdxTrack1[itTrack].numberOfMeasurements();
      dedx1_ = dEdxTrack1[itTrack].dEdx();
      dedx2_ = dEdxTrack2[itTrack].dEdx();
      dedx3_ = dEdxTrack3[itTrack].dEdx();
    } catch ( cms::Exception& ) {
      dedxNoM_ = 0;
      dedx1_ = 0.;
      dedx2_ = 0.;
      dedx3_ = 0.;
    }
    p_     = itTrack->p();
    pt_    = itTrack->pt();
    eta_   = itTrack->eta();
    phi_   = itTrack->phi();
    charge_    = itTrack->charge();
    quality_   = itTrack->qualityMask();
    foundhits_ = itTrack->found();
    lostHits_  = itTrack->lost();
    foundhitsStrips_ =  itTrack->hitPattern().numberOfValidStripHits();
    foundhitsPixels_ =  itTrack->hitPattern().numberOfValidPixelHits();
    losthitsStrips_  =  itTrack->hitPattern().numberOfLostStripHits(reco::HitPattern::TRACK_HITS);
    losthitsPixels_  =  itTrack->hitPattern().numberOfLostPixelHits(reco::HitPattern::TRACK_HITS);
    nLayers_         =  itTrack->hitPattern().trackerLayersWithMeasurement();
    chi2_    = itTrack->chi2();
    ndof_    = itTrack->ndof();
    dz_      = itTrack->dz();
    dzerr_   = itTrack->dzError();
    dzCorr_  = itTrack->dz(beamSpot);
    dxy_     = itTrack->dxy();
    dxyerr_  = itTrack->dxyError();
    dxyCorr_ = itTrack->dxy(beamSpot);
    pterr_   = itTrack->ptError();
    etaerr_  = itTrack->etaError();
    phierr_  = itTrack->phiError();
    qoverp_  = itTrack->qoverp();
    xPCA_    = itTrack->vertex().x();
    yPCA_    = itTrack->vertex().y();
    zPCA_    = itTrack->vertex().z();
    
    try { // only one track collection (at best) is connected to the main vertex
      if(!vertexColl.empty() && !vertexColl.begin()->isFake()) {
	trkWeightpvtx_ = vertexColl.begin()->trackWeight(itTrack);
      } else
	trkWeightpvtx_ = 0.;
    } catch ( cms::Exception& ) {
      trkWeightpvtx_ = 0.;
    }
     // compute the fraction of low probability pixels... will be added to the event tree
    for(trackingRecHit_iterator it = itTrack->recHitsBegin(); it!=itTrack->recHitsEnd(); ++it) {
      const TrackingRecHit* hit = &(**it);
      const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(hit);
      if(pixhit) {
	DetId detId = pixhit->geographicalId();
	if(detId.subdetId() == PixelSubdetector::PixelBarrel) {
	  ++n_hits_barrel;
	  double proba = pixhit->clusterProbability(0);
	  if(proba<=0.0) ++n_hits_lowprob;
	}
      }
    }
    if(functionality_tracks_)
      tracks_->Fill();
  }
  lowPixelProbabilityFraction_ = n_hits_barrel>0 ? (float)n_hits_lowprob/n_hits_barrel : -1.;
  
  // iterate over clusters
  unsigned int localCounter = 0;
  nclusters_ = 0;
  std::vector<double> clusterOntrackAngles = onTrackAngles(clusters,trajectoryCollection);
  std::vector<double>::const_iterator angleIt = clusterOntrackAngles.begin();     

  for (edmNew::DetSetVector<SiStripCluster>::const_iterator DSViter=clusters->begin(); DSViter!=clusters->end();DSViter++ ) {
    edmNew::DetSet<SiStripCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiStripCluster>::const_iterator end  =DSViter->end();
    unsigned int detid = DSViter->id();
    nclusters_ += DSViter->size();
     
    for(edmNew::DetSet<SiStripCluster>::const_iterator iter=begin;iter!=end;++iter,++angleIt,++localCounter) {
      siStripClusterInfo_.setCluster(*iter, detid);
      // general quantities
      trackid_ = stripClusterOntrackIndices[localCounter];
      onTrack_    = (trackid_ != (unsigned int)-1);
      clWidth_    = siStripClusterInfo_.width();
      clPosition_ = siStripClusterInfo_.baryStrip();
      angle_      = *angleIt;
      thickness_  = ((((DSViter->id()>>25)&0x7f)==0xd) ||
		     ((((DSViter->id()>>25)&0x7f)==0xe) && (((DSViter->id()>>5)&0x7)>4))) ? 500 : 300;
      stripLength_ = static_cast<const StripGeomDetUnit*>(tracker_->idToDet(detid))->specificTopology().stripLength();
      int nstrips  = static_cast<const StripGeomDetUnit*>(tracker_->idToDet(detid))->specificTopology().nstrips();
      maxCharge_   = siStripClusterInfo_.maxCharge();
      // signal and noise with gain corrections
      clNormalizedCharge_ = siStripClusterInfo_.charge() ;
      clNormalizedNoise_  = siStripClusterInfo_.noiseRescaledByGain() ;
      clSignalOverNoise_  = siStripClusterInfo_.signalOverNoise() ;
      // signal and noise with gain corrections and angle corrections
      clCorrectedCharge_ = clNormalizedCharge_ * fabs(cos(angle_)); // corrected for track angle
      clCorrectedSignalOverNoise_ = clSignalOverNoise_ * fabs(cos(angle_)); // corrected for track angle
      maxChargeCorrected_ = maxCharge_*fabs(cos(angle_));
      // signal and noise without gain corrections
      clBareNoise_  = siStripClusterInfo_.noise();
      clBareCharge_ = clSignalOverNoise_*clBareNoise_;
       // global position
      const StripGeomDetUnit* sgdu = static_cast<const StripGeomDetUnit*>(tracker_->idToDet(detid));
      Surface::GlobalPoint gp = sgdu->surface().toGlobal(sgdu->specificTopology().localPosition(MeasurementPoint(clPosition_,0)));
      globalX_ = gp.x();
      globalY_ = gp.y();
      globalZ_ = gp.z();
       // cabling
      detid_ = detid;
      lldChannel_ = 1+(int(floor(iter->barycenter()))/256);
      if(lldChannel_==2 && nstrips==512) lldChannel_=3;
      // Not yet started with the first random map
      if(delayrange_ == 0)
	delay_ = -1;
      else{
	delay_ = delayDetIdMap_[detid_]+delaystepsize_*delaystep_*(25./24.);
	if(std::round(delay_) > delayrange_) {
	  delay_ -= (2*delayrange_+delaystepsize_)*(25./24.);
	}
      }
      if((functionality_offtrackClusters_&&!onTrack_)||(functionality_ontrackClusters_&&onTrack_)) 
	clusters_->Fill();
    }
  }
  
  // fill event tree
  if(functionality_events_)
    event_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------

void TrackerDpgAnalysis::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){  

  //Retrieve tracker topology from geometry
  const TrackerTopology* tTopo = &iSetup.getData(tTopoToken_);

  //geometry
  tracker_ = &iSetup.getData(tkGeomToken_);
  
  // cabling I (readout)
  cabling_ = &iSetup.getData(fedCablingToken_);
  auto feds = cabling_->fedIds() ;

  // read the delay offsets for each device from input files
  // this is only for the so-called "random delay" run
  std::map<unsigned int,float> delayMap = delay(delayFileName_);
  if(delayMap.empty()){
    edm::LogError("BadConfig") << " The delay file does not exist. The delay map will not be filled properly."                                                                                     
			       << std::endl << " Looking for " << delayFileName_ << "."                                                                                                             
			       << std::endl << " Please specify valid filenames through the DelayFileName FileInPath parameter.";                                                                  
    return;
  }

  edm::ParameterSet pset;
  pset.addUntrackedParameter<bool>("logScale",false);
  TrackerMap themap(pset);

  for(auto fedid = feds.begin();fedid<feds.end();++fedid) { 
    auto connections = cabling_->fedConnections(*fedid);
    for(auto conn=connections.begin();conn<connections.end();++conn) {
      // Fill the "old" map to be used for lookup during analysis
      if(conn->isConnected())
	connections_.insert(std::make_pair(conn->detId(),new FedChannelConnection(*conn)));
      // Fill the standalone tree (once for all)
      if(conn->isConnected()) {
	detid_ = conn->detId();
	moduleName_ = toStringName(detid_,tTopo);
	moduleId_ = toStringId(detid_);
	lldChannel_ = conn->lldChannel();
	dcuId_ = conn->dcuId();
	fecCrate_ = conn->fecCrate();
	fecSlot_ = conn->fecSlot();
	fecRing_ = conn->fecRing();
	ccuAdd_ = conn->ccuAddr();
	ccuChan_ = conn->ccuChan();
	fedId_ = conn->fedId();
	fedCh_ = conn->fedCh();
	fiberLength_ = conn->fiberLength();
	delayrandom_ = delayMap[dcuId_];
	delayDetIdMap_[detid_] = delayrandom_;
	const StripGeomDetUnit* sgdu = static_cast<const StripGeomDetUnit*>(tracker_->idToDet(detid_));
	Surface::GlobalPoint gp = sgdu->surface().toGlobal(LocalPoint(0,0));
	globalX_ = gp.x();
	globalY_ = gp.y();
	globalZ_ = gp.z();
	readoutmap_->Fill();
	themap.fill_current_val(detid_,delayrandom_);
      }
    }
  }

  // Naming convention from DAQ TK is TrackerMapDelay_Run<number>_pll.csv or TrackerMapDelay_Run<number>_fed.csv .. extract run number
  std::stringstream delayName (delayFileName_);
  std::string segment;
  std::vector<std::string> seglist;
  while(std::getline(delayName, segment, '/'))
    seglist.push_back(segment);
  std::stringstream fileName (seglist.back());
  seglist.clear();
  while(std::getline(fileName,segment,'_'))
    seglist.push_back(segment);
  TString mapname (seglist.at(1));
  mapname = Form("DelayMap_"+mapname);
  themap.setTitle("Delays");
  themap.save(false,-10,10,std::string(mapname+".png"),1400,800);
  
  // HLT
  bool flag = false;
  HLTConfigProvider hltConfig;
  hltConfig.init(iRun, iSetup, HLTTag_.process(), flag);
  if(not hltNames_.empty()){
    for(size_t itrig = 0; itrig < hltNames_.size(); itrig++){
      auto it = std::find(hltConfig.triggerNames().begin(),hltConfig.triggerNames().end(),hltNames_.at(itrig));
      if(it != hltConfig.triggerNames().end())
	triggerPathsMap[*it] = (it-hltConfig.triggerNames().begin());
    }
  }
  else{
    for(size_t i = 0; i < hltConfig.triggerNames().size(); i++){
      std::string pathName = hltConfig.triggerNames()[i];
      triggerPathsMap[pathName] = i;      
    }
  } 
}

void TrackerDpgAnalysis::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}

void TrackerDpgAnalysis::beginJob(){

  // create output
  edm::Service<TFileService> fileService;
  TFileDirectory* dir = new TFileDirectory(fileService->mkdir("trackerDPG"));

  globaltrackid_  = pset_.getParameter<unsigned int>("InitalCounter");
  globalvertexid_ = pset_.getParameter<unsigned int>("InitalCounter");
  
  // create a tree for the events
  event_ = dir->make<TTree>("events","event information");
  event_->Branch("eventid",&eventid_,"eventid/i");
  event_->Branch("runid",&runid_,"runid/i");
  event_->Branch("lumi",&lumisec_,"lumi/i");
  event_->Branch("orbit",&orbit_,"orbit/i");
  event_->Branch("orbitL1",&orbitL1_,"orbitL1/i");
  event_->Branch("bx",&bx_,"bx/i");
  event_->Branch("store",&store_,"store/i");
  event_->Branch("time",&time_,"time/i");
  event_->Branch("physicsDeclared",&physicsDeclared_,"physicsDeclared/s");
  event_->Branch("hltTriggerBits","std::vector<bool>",&hltTriggerBits_);
  event_->Branch("hltTriggerNames","std::vector<std::string>",&hltTriggerNames_);
  event_->Branch("nclusters",&nclusters_,"nclusters/i");
  event_->Branch("nclustersOntrack",&nclustersOntrack_,"nclustersOntrack/i");
  event_->Branch("bsX0",&bsX0_,"bsX0/F");
  event_->Branch("bsY0",&bsY0_,"bsY0/F");
  event_->Branch("bsZ0",&bsZ0_,"bsZ0/F");
  event_->Branch("bsSigmaZ",&bsSigmaZ_,"bsSigmaZ/F");
  event_->Branch("bsDxdz",&bsDxdz_,"bsDxdz/F");
  event_->Branch("bsDydz",&bsDydz_,"bsDydz/F");
  event_->Branch("nVertices",&nVertices_,"nVertices/i");
  event_->Branch("MagneticField",&fBz_,"MagneticField/F");
  event_->Branch("ntracks",&ntracks_,"ntracks/i");
  event_->Branch("ntrajs",&ntrajs_,"ntrajs/i");
  event_->Branch("lowPixelProbabilityFraction",&lowPixelProbabilityFraction_,"lowPixelProbabilityFraction/F");
  event_->Branch("runtype",&runtype_,"runtype/i");
  event_->Branch("fedreadout",&fedreadout_,"fedreadout/i");
  event_->Branch("delayrange",&delayrange_,"delayrange/i");
  event_->Branch("delaystepsize",&delaystepsize_,"delaystepsize/i");
  event_->Branch("delaystep",&delaystep_,"delaystep/i");
  
  // create a tree for the vertices
  vertices_ = dir->make<TTree>("vertices","vertex information");
  vertices_->Branch("vertexid",&globalvertexid_,"vertexid/i");
  vertices_->Branch("eventid",&eventid_,"eventid/i");
  vertices_->Branch("runid",&runid_,"runid/i");
  vertices_->Branch("lumi",&lumisec_,"lumi/i");
  vertices_->Branch("nTracks",&nTracks_pvtx_,"nTracks/i");
  vertices_->Branch("sumptsq",&sumptsq_pvtx_,"sumptsq/F");
  vertices_->Branch("isValid",&isValid_pvtx_,"isValid/O");
  vertices_->Branch("isFake",&isFake_pvtx_,"isFake/O");
  vertices_->Branch("recx",&recx_pvtx_,"recx/F");
  vertices_->Branch("recy",&recy_pvtx_,"recy/F");
  vertices_->Branch("recz",&recz_pvtx_,"recz/F");
  vertices_->Branch("recx_err",&recx_err_pvtx_,"recx_err/F");
  vertices_->Branch("recy_err",&recy_err_pvtx_,"recy_err/F");
  vertices_->Branch("recz_err",&recz_err_pvtx_,"recz_err/F");
  
  // create a TTree for clusters
  clusters_ = dir->make<TTree>("clusters","cluster information");
  clusters_->Branch("eventid",&eventid_,"eventid/i");
  clusters_->Branch("runid",&runid_,"runid/i");
  clusters_->Branch("lumi",&lumisec_,"lumi/i");
  clusters_->Branch("trackid",&trackid_,"trackid/i");
  clusters_->Branch("onTrack",&onTrack_,"onTrack/O");
  clusters_->Branch("clWidth",&clWidth_,"clWidth/F");
  clusters_->Branch("clPosition",&clPosition_,"clPosition/F");
  clusters_->Branch("clglobalX",&globalX_,"clglobalX/F");
  clusters_->Branch("clglobalY",&globalY_,"clglobalY/F");
  clusters_->Branch("clglobalZ",&globalZ_,"clglobalZ/F");
  clusters_->Branch("angle",&angle_,"angle/F");
  clusters_->Branch("thickness",&thickness_,"thickness/F");
  clusters_->Branch("maxCharge",&maxCharge_,"maxCharge/F");
  clusters_->Branch("maxChargeCorrected",&maxChargeCorrected_,"maxChargeCorrected/F");
  clusters_->Branch("clNormalizedCharge",&clNormalizedCharge_,"clNormalizedCharge/F");
  clusters_->Branch("clNormalizedNoise",&clNormalizedNoise_,"clNormalizedNoise/F");
  clusters_->Branch("clSignalOverNoise",&clSignalOverNoise_,"clSignalOverNoise/F");
  clusters_->Branch("clCorrectedCharge",&clCorrectedCharge_,"clCorrectedCharge/F");
  clusters_->Branch("clCorrectedSignalOverNoise",&clCorrectedSignalOverNoise_,"clCorrectedSignalOverNoise/F");
  clusters_->Branch("clBareCharge",&clBareCharge_,"clBareCharge/F");
  clusters_->Branch("clBareNoise",&clBareNoise_,"clBareNoise/F");
  clusters_->Branch("stripLength",&stripLength_,"stripLength/F");
  clusters_->Branch("detid",&detid_,"detid/i");
  clusters_->Branch("lldChannel",&lldChannel_,"lldChannel/i");
  clusters_->Branch("delay",&delay_,"delay/F");
  
  // create a tree for tracks   
  tracks_ = dir->make<TTree>("tracks","tracks information");
  tracks_->Branch("eventid",&eventid_,"eventid/i");
  tracks_->Branch("runid",&runid_,"runid/i");
  tracks_->Branch("lumi",&lumisec_,"lumi/i");
  tracks_->Branch("vertexid",&vertexid_,"vertexid/i");
  tracks_->Branch("trackid",&globaltrackid_,"trackid/i");
  tracks_->Branch("chi2",&chi2_,"chi2/F");
  tracks_->Branch("eta",&eta_,"eta/F");
  tracks_->Branch("etaerr",&etaerr_,"etaerr/F");
  tracks_->Branch("phi",&phi_,"phi/F");
  tracks_->Branch("phierr",&phierr_,"phierr/F");
  tracks_->Branch("dedx1",&dedx1_,"dedx1/F");
  tracks_->Branch("dedx2",&dedx2_,"dedx2/F");
  tracks_->Branch("dedx3",&dedx3_,"dedx3/F");
  tracks_->Branch("dedxNoM",&dedxNoM_,"dedxNoM/i");
  tracks_->Branch("charge",&charge_,"charge/F");
  tracks_->Branch("quality",&quality_,"quality/i");
  tracks_->Branch("foundhits",&foundhits_,"foundhits/i");
  tracks_->Branch("lostHits",&lostHits_,"lostHits/i");
  tracks_->Branch("foundhitsStrips",&foundhitsStrips_,"foundhitsStrips/i");
  tracks_->Branch("foundhitsPixels",&foundhitsPixels_,"foundhitsPixels/i");
  tracks_->Branch("losthitsStrips",&losthitsStrips_,"losthitsStrips/i");
  tracks_->Branch("losthitsPixels",&losthitsPixels_,"losthitsPixels/i");
  tracks_->Branch("p",&p_,"p/F");
  tracks_->Branch("pt",&pt_,"pt/F");
  tracks_->Branch("pterr",&pterr_,"pterr/F");
  tracks_->Branch("ndof",&ndof_,"ndof/i");
  tracks_->Branch("dz",&dz_,"dz/F");
  tracks_->Branch("dzerr",&dzerr_,"dzerr/F");
  tracks_->Branch("dzCorr",&dzCorr_,"dzCorr/F");
  tracks_->Branch("dxy",&dxy_,"dxy/F");
  tracks_->Branch("dxyerr",&dxyerr_,"dxyerr/F");
  tracks_->Branch("dxyCorr",&dxyCorr_,"dxyCorr/F");
  tracks_->Branch("qoverp",&qoverp_,"qoverp/F");
  tracks_->Branch("xPCA",&xPCA_,"xPCA/F");
  tracks_->Branch("yPCA",&yPCA_,"yPCA/F");
  tracks_->Branch("zPCA",&zPCA_,"zPCA/F");
  tracks_->Branch("nLayers",&nLayers_,"nLayers/i");
  tracks_->Branch("trkWeightpvtx",&trkWeightpvtx_,"trkWeightpvtx/F");
  
  // cabling
  readoutmap_ = dir->make<TTree>("readoutMap","cabling map");
  readoutmap_->Branch("detid",&detid_,"detid/i");
  readoutmap_->Branch("dcuId",&dcuId_,"dcuId/i");
  readoutmap_->Branch("fecCrate",&fecCrate_,"fecCrate/i");
  readoutmap_->Branch("fecSlot",&fecSlot_,"fecSlot/i");
  readoutmap_->Branch("fecRing",&fecRing_,"fecRing/i");
  readoutmap_->Branch("ccuAdd",&ccuAdd_,"ccuAdd/i");
  readoutmap_->Branch("ccuChan",&ccuChan_,"ccuChan/i");
  readoutmap_->Branch("lldChannel",&lldChannel_,"lldChannel/i");
  readoutmap_->Branch("fedId",&fedId_,"fedId/i");
  readoutmap_->Branch("fedCh",&fedCh_,"fedCh/i");
  readoutmap_->Branch("fiberLength",&fiberLength_,"fiberLength/i");
  readoutmap_->Branch("moduleName",&moduleName_);
  readoutmap_->Branch("moduleId",&moduleId_);
  readoutmap_->Branch("delaymap",&delayrandom_,"delaymap/F");
  readoutmap_->Branch("globalX",&globalX_,"globalX/F");
  readoutmap_->Branch("globalY",&globalY_,"globalY/F");
  readoutmap_->Branch("globalZ",&globalZ_,"globalZ/F");
}


// ------------ method called once each job just after ending the event loop  ------------
void TrackerDpgAnalysis::endJob() {}  

  // build reverse map track -> vertex
std::map<size_t,int> TrackerDpgAnalysis::inVertex(const reco::VertexCollection& vertices, unsigned int firstVertex){
  std::map<size_t,int> output;
  unsigned int vertexid = firstVertex;
  for(reco::VertexCollection::const_iterator v = vertices.begin(); v!=vertices.end(); ++v,++vertexid) {
    reco::Vertex::trackRef_iterator it = v->tracks_begin();
    reco::Vertex::trackRef_iterator lastTrack = v->tracks_end();
    for(;it!=lastTrack;++it) {
      output[it->key()] = vertexid;
    }
  }
  return output;
}

void TrackerDpgAnalysis::insertMeasurement(std::multimap<const unsigned int,std::pair<LocalPoint,double> >& collection,const TransientTrackingRecHit* hit , double tla){
  if(!hit) return;
  const SiTrackerMultiRecHit* multihit=dynamic_cast<const SiTrackerMultiRecHit*>(hit);
  const SiStripRecHit2D* singlehit=dynamic_cast<const SiStripRecHit2D*>(hit);
  const SiStripRecHit1D* hit1d=dynamic_cast<const SiStripRecHit1D*>(hit);
  if(hit1d) { //...->33X
    collection.insert(std::make_pair(hit1d->geographicalId().rawId(),std::make_pair(hit1d->localPosition(),tla)));
  } else if(singlehit) { // 41X->...
    collection.insert(std::make_pair(singlehit->geographicalId().rawId(),std::make_pair(singlehit->localPosition(),tla)));
  }
  else if(multihit){
    std::vector< const TrackingRecHit * > childs = multihit->recHits();
    for(std::vector<const TrackingRecHit*>::const_iterator it=childs.begin();it!=childs.end();++it) {
      insertMeasurement(collection,dynamic_cast<const TrackingRecHit*>(*it),tla);
    }
  }
}


void TrackerDpgAnalysis::insertMeasurement(std::multimap<const unsigned int,std::pair<int, int> >& collection,const TrackingRecHit* hit , int trackid){
  if(!hit) return;
  const SiTrackerMultiRecHit* multihit=dynamic_cast<const SiTrackerMultiRecHit*>(hit);
  const SiStripRecHit2D* singlehit=dynamic_cast<const SiStripRecHit2D*>(hit);
  const SiStripRecHit1D* hit1d=dynamic_cast<const SiStripRecHit1D*>(hit);
  if(hit1d) { // 41X->...
    collection.insert(std::make_pair(hit1d->geographicalId().rawId(),std::make_pair(int(hit1d->cluster()->barycenter()),trackid)));
  } else if(singlehit) { //...->33X
    collection.insert(std::make_pair(singlehit->geographicalId().rawId(),std::make_pair(int(singlehit->cluster()->barycenter()),trackid)));
  }
  else if(multihit){
    std::vector< const TrackingRecHit * > childs = multihit->recHits();
    for(std::vector<const TrackingRecHit*>::const_iterator it=childs.begin();it!=childs.end();++it) {
      insertMeasurement(collection,*it,trackid);
    }
  }
}

std::vector<int> TrackerDpgAnalysis::onTrack(edm::Handle<edmNew::DetSetVector<SiStripCluster> >& clusters,
					     const reco::TrackCollection& trackVec, unsigned int firstTrack) {
  std::vector<int> result;
  // first, build a list of positions and trackid on tracks
  std::multimap<const unsigned int,std::pair<int,int> > onTrackPositions;
  unsigned int trackid = firstTrack;
  for(reco::TrackCollection::const_iterator itTrack = trackVec.begin(); itTrack!=trackVec.end();++itTrack,++trackid) {
    for(trackingRecHit_iterator it = itTrack->recHitsBegin(); it!=itTrack->recHitsEnd(); ++it) {
      const TrackingRecHit* hit = &(**it);
      insertMeasurement(onTrackPositions,hit,trackid);
    }
  }
  // then loop over the clusters to check
  int thetrackid = -1;
  for (edmNew::DetSetVector<SiStripCluster>::const_iterator DSViter=clusters->begin(); DSViter!=clusters->end();DSViter++ ) {
    edmNew::DetSet<SiStripCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiStripCluster>::const_iterator end  =DSViter->end();
    std::pair< std::multimap<unsigned int,std::pair<int,int> >::const_iterator,
               std::multimap<unsigned int,std::pair<int,int> >::const_iterator> range =
      onTrackPositions.equal_range(DSViter->id());
    for(edmNew::DetSet<SiStripCluster>::const_iterator iter=begin;iter!=end;++iter) {
      thetrackid = -1;
      for(std::multimap<unsigned int,std::pair<int,int> >::const_iterator cl = range.first; cl!= range.second; ++cl) {
        if(fabs(cl->second.first-iter->barycenter())<2) {
	  thetrackid = cl->second.second;
	}
      }
      result.push_back(thetrackid);
    }
  }
  return result;
}

std::vector<double> TrackerDpgAnalysis::onTrackAngles(edm::Handle<edmNew::DetSetVector<SiStripCluster> >& clusters,
						      const std::vector<Trajectory>& trajVec ){
  std::vector<double> result;
  // first, build a list of positions and angles on trajectories
  std::multimap<const unsigned int,std::pair<LocalPoint,double> > onTrackPositions;
  for(std::vector<Trajectory>::const_iterator traj = trajVec.begin(); traj< trajVec.end(); ++traj) {
    Trajectory::DataContainer measurements = traj->measurements();
    for(Trajectory::DataContainer::iterator meas = measurements.begin(); meas!= measurements.end(); ++meas) {
      double tla = meas->updatedState().localDirection().theta();
      insertMeasurement(onTrackPositions,&(*(meas->recHit())),tla);
    }
  }
  // then loop over the clusters to check
  double angle = 0.;
  for (edmNew::DetSetVector<SiStripCluster>::const_iterator DSViter=clusters->begin(); DSViter!=clusters->end();DSViter++ ) {
    edmNew::DetSet<SiStripCluster>::const_iterator begin=DSViter->begin();
    edmNew::DetSet<SiStripCluster>::const_iterator end  =DSViter->end();
    std::pair< std::multimap<unsigned int,std::pair<LocalPoint,double> >::const_iterator,
               std::multimap<unsigned int,std::pair<LocalPoint,double> >::const_iterator> range =
      onTrackPositions.equal_range(DSViter->id());
    const GeomDetUnit* gdu = static_cast<const GeomDetUnit*>(tracker_->idToDet(DSViter->id()));
    for(edmNew::DetSet<SiStripCluster>::const_iterator iter=begin;iter!=end;++iter) {
      angle = 0.;
      for(std::multimap<unsigned int,std::pair<LocalPoint,double> >::const_iterator cl = range.first; cl!= range.second; ++cl) {
        if(fabs(gdu->topology().measurementPosition(cl->second.first).x()-iter->barycenter())<2) {
	  angle = cl->second.second;
	}
      }
      result.push_back(angle);
    }
  }
  return result;
}


std::string TrackerDpgAnalysis::toStringName(unsigned int rawid, const TrackerTopology* tTopo) {
  SiStripDetId detid(rawid);
  std::string out;
  std::stringstream output;
  switch(detid.subDetector()) {
  case 3:
    {
      output << "TIB";
      output << (tTopo->tibIsZPlusSide(rawid) ? "+" : "-");
      output << " layer ";
      output << tTopo->tibLayer(rawid);
      output << ", string ";
      output << tTopo->tibString(rawid);
      output << (tTopo->tibIsExternalString(rawid) ? " external" : " internal");
      output << ", module ";
      output << tTopo->tibModule(rawid);
      if(tTopo->tibIsDoubleSide(rawid)) {
	output << " (double)";
      } else {
	output << (tTopo->tibIsRPhi(rawid) ? " (rphi)" : " (stereo)");
      }
      break;
    }
  case 4:
    {
      output << "TID";       
      output << (tTopo->tidIsZPlusSide(rawid) ? "+" : "-");
      output << " disk ";
      output << tTopo->tidWheel(rawid);
      output << ", ring ";
      output << tTopo->tidRing(rawid);
      output << (tTopo->tidIsFrontRing(rawid) ? " front" : " back");
      output << ", module ";
      output << tTopo->tidModule(rawid);
      if(tTopo->tidIsDoubleSide(rawid)) {
	output << " (double)";
      } else {
	output << (tTopo->tidIsRPhi(rawid) ? " (rphi)" : " (stereo)");
      }
      break;
    }
  case 5:
    {
      output << "TOB";
      output << (tTopo->tobIsZPlusSide(rawid) ? "+" : "-");
      output << " layer ";
      output << tTopo->tobLayer(rawid);
      output << ", rod ";
      output << tTopo->tobRod(rawid);
      output << ", module ";
      output << tTopo->tobModule(rawid);
      if(tTopo->tobIsDoubleSide(rawid)) {
	output << " (double)";
      } else {
	output << (tTopo->tobIsRPhi(rawid) ? " (rphi)" : " (stereo)");
      }
      break;
    }
  case 6:
    {
      output << "TEC";
      output << (tTopo->tecIsZPlusSide(rawid) ? "+" : "-");
      output << " disk ";
      output << tTopo->tecWheel(rawid);
      output << " sector ";
      output << tTopo->tecPetalNumber(rawid);
      output << (tTopo->tecIsFrontPetal(rawid) ? " Front Petal" : " Back Petal");
      output << ", module ";
      output << tTopo->tecRing(rawid);
      output << tTopo->tecModule(rawid);
      if(tTopo->tecIsDoubleSide(rawid)) {
	output << " (double)";
      } else {
	output << (tTopo->tecIsRPhi(rawid) ? " (rphi)" : " (stereo)");
      }
      break;
    }
  default:
    {
      output << "UNKNOWN";
    }
  }
  out = output.str();
  return out;
}

std::string TrackerDpgAnalysis::toStringId(unsigned int rawid) {
  std::string out;
  std::stringstream output;
  output << rawid << " (0x" << std::hex << rawid << std::dec << ")";
  out = output.str();
  return out;
}

double TrackerDpgAnalysis::sumPtSquared(const reco::Vertex & v)  {
  double sum = 0.;
  double pT;
  for (reco::Vertex::trackRef_iterator it = v.tracks_begin(); it != v.tracks_end(); it++) {
    pT = (**it).pt();
    sum += pT*pT;
  }
  return sum;
}

// parser of delay file written for the manual delay scan when DB state is uploaded
std::map<unsigned int,float> TrackerDpgAnalysis::delay(const std::string & file) {

  // prepare output
  unsigned int dcuid;
  float delay;
  std::map<unsigned int,float> delayMap;
  
  std::ifstream delayInput(file);
  if(delayInput.is_open()) {
    std::string line;
    while(std::getline(delayInput, line)){
      std::istringstream s (line);
      std::string field;
      std::vector<std::string> entries;
      while(std::getline(s,field,','))
	entries.push_back(field);
      std::istringstream dcuidstr (entries.front());
      std::istringstream delaystr (entries.back());
      dcuidstr >> dcuid;
      delaystr >> delay;
      delayMap[dcuid] = delay;
    }
  }
  else
    edm::LogError("BadConfig") << " Delay file cannot be parsed properly .. failing opening therefore please cross-check the setup";
  return delayMap;
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrackerDpgAnalysis);
