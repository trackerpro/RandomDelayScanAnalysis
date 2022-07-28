#include <vector>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEventList.h>
#include <TCut.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLeaf.h>

#include "../delayUtils.h"

#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"

using namespace RooFit;
 
// parametrize profile with a gaussian shape
static bool isGaussian = true;
// reduce the number of events by
static int reductionFactor = 1;
// min number of filled bins
static int minFilledBin    = 4;
// max allowed delay
static float peakBoundary  = 9.36;
// min amplitude
static float amplitudeMin  = 40.; 
// min signficance delay/sigma_delay
static float significance  = 2;
// min delay to apply significance cut
static float minDelayForSignificance = 2;
// max sigma allowed
static float maxSigma = 20;
// decide whether use histograms values or fitted ones
static bool doFitOfClusterShape = true;
// decide whether use histograms values or fitted ones
static bool storeFitOfClusterShape = true;
// fit range in quantiles
static std::pair<double,double> fit_quantiles = std::make_pair<double,double>(0.001,0.995);
 
/// function that runs on the evnet and produce profiles for layers
void ChannnelPlots(const std::vector<std::string> & fileNames,
		   const std::vector<long int> & nEvents,
		   TTree* corrections,		   
		   std::map<unsigned int,std::map<int,TH1F* > > channelMap,
		   std::map<unsigned int,TH1F* > & channelMapMean,
		   std::map<unsigned int,TH1F* > & channelMapMPV,
		   std::map<unsigned int,int > & fitStatusMean,
		   std::map<unsigned int,int > & fitStatusMPV,
		   const string & observable,
		   const string & outputDIR) {

  std::cout<<"###############################################################"<<std::endl;
  std::cout << "Preparing Layer plot for the for all the different channels "<< std::endl; 
  std::cout<<"##############################################################"<<std::endl;

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;

  if(observable == "clCorrectedSignalOverNoise" or observable == "clSignalOverNoise")
    amplitudeMin = 3.;

  // create vectors for the different Profiles
  float yMin   = 0, yMax = 0;
  int   nBinsY = 0;
  vector<double> delayBins;
  setLimitsAndBinning("delay",delayBins);
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);

  int correction;
  corrections->SetBranchStatus("*",kFALSE);
  corrections->SetBranchStatus("detid",kTRUE);
  corrections->SetBranchStatus("correction",kTRUE);
  corrections->SetBranchAddress("correction",&correction);

  for(size_t ifile = 0; ifile < fileNames.size(); ifile++){
    std::cout<<"Looping on clusters of file "<<ifile<<" out of "<<fileNames.size()<<std::endl;
    auto file = TFile::Open(fileNames.at(ifile).c_str(),"READ");
    auto tree = (TTree*) file->FindObjectAny("clusters");
    // set branches for the cluster, readoutmap and no corrections trees
    unsigned int detid;
    float    clCorrectedSignalOverNoise, clSignalOverNoise, obs, delay;
    tree->SetBranchStatus("*",kFALSE);
    tree->SetBranchStatus("detid",kTRUE);
    tree->SetBranchStatus("clCorrectedSignalOverNoise",kTRUE);
    tree->SetBranchStatus("clSignalOverNoise",kTRUE);
    tree->SetBranchStatus("delay",kTRUE);
    tree->SetBranchStatus(observable.c_str(),kTRUE);
    tree->SetBranchAddress("detid",&detid);
    tree->SetBranchAddress("clCorrectedSignalOverNoise",&clCorrectedSignalOverNoise);
    tree->SetBranchAddress("clSignalOverNoise",&clSignalOverNoise);
    tree->SetBranchAddress("delay",&delay);
    tree->SetBranchAddress(observable.c_str(),&obs);

    for(long int iEvent = 0; iEvent < nEvents.at(ifile)/reductionFactor; iEvent++){

      cout.flush();
      if(iEvent % 100000 == 0) 
	cout<<"\r"<<"iEvent "<<100*double(iEvent)/(nEvents.at(ifile)/reductionFactor)<<" % ";
      // take the event                                                                                                                                                  
      tree->GetEntry(iEvent);
      // take the correction from the detid
      corrections->GetEntryWithIndex(detid);
      int delaybin = std::round(delay) - -std::round(correction);
      if(channelMap[detid][delaybin] == 0 or channelMap[detid][std::round(delay)-std::round(correction)] == NULL){
	channelMap[detid][delaybin] =  new TH1F(Form("detid_%d_delay_%d",detid,delaybin),"",nBinsY,yMin,yMax);
	channelMap[detid][delaybin]->SetDirectory(0);
	channelMap[detid][delaybin]->Sumw2();
      }
      float value = 0;
      if(observable == "maxCharge")
	value = obs*(clCorrectedSignalOverNoise)/(clSignalOverNoise);
      else 	
	value = obs;
      channelMap[detid][delaybin]->Fill(value);
    }
    std::cout<<std::endl;
    file->Close(); // to prevent high RAM occupancy    
  }
  std::cout<<"Loop on events terminated"<<std::endl;  

  std::cout<<"Analyze each single channel --> MPV and Mean "<<std::endl;
  long int iChannel = 0;
  TFile* outfile = NULL;
  if(storeFitOfClusterShape){
    outfile = new TFile((outputDIR+"/distributions_"+observable+".root").c_str(),"RECREATE");
    outfile->cd();
  }

  TCanvas* canvas = new TCanvas("canvas","",600,625);  
  // make the convuluated fit as suggested in https://root.cern.ch/root/html/tutorials/fit/langaus.C.html
  long int iBadChannelFit = 0;
  for(auto iMap : channelMap){
    iChannel++;
    cout.flush();
    if(iChannel % 100 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(channelMap.size())<<" % ";
    channelMapMPV[iMap.first]  = new TH1F(Form("detid_%d_mpv",iMap.first),"",delayBins.size()-1,&delayBins[0]);    
    channelMapMean[iMap.first] = new TH1F(Form("detid_%d_mean",iMap.first),"",delayBins.size()-1,&delayBins[0]);    
    channelMapMPV[iMap.first]->SetDirectory(0);
    channelMapMean[iMap.first]->SetDirectory(0);

    for(std::map<int,TH1F*>::const_iterator itHisto = iMap.second.begin(); itHisto != iMap.second.end(); itHisto++){
      // to fill the Mean --> use always the histogram properties given its binning
      channelMapMean[iMap.first]->SetBinContent(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMean());
      // as error use the error on the mean as a proxy
      channelMapMean[iMap.first]->SetBinError(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMeanError());

      if(not doFitOfClusterShape){ // take the MPV as the maximum of the histograms and store histogram in case
	channelMapMPV[iMap.first]->SetBinContent(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin()));
	// as error use the error on the mean as a proxy
	channelMapMPV[iMap.first]->SetBinError(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMeanError());	
	// to fill the MPV: delay found by means of FindBin, take the bin center of the maximum bin --> is ok as soon as the histogram has enough bins
	if(storeFitOfClusterShape and iChannel % 20 == 0){
	  if(not outfile->GetDirectory(Form("Delay_%d",itHisto->first)))
	    outfile->mkdir(Form("Delay_%d",itHisto->first));
	  outfile->cd(Form("Delay_%d",itHisto->first));
	  itHisto->second->Write();
	  outfile->cd();	            	  
	}
      }  
      else{
	
	float xMin   = 0, xMax = 0;
	int   nBinsX = 0;
	setLimitsAndBinning(observable,xMin,xMax,nBinsX);	

	// make the observable                                                                                                                                                                 
	RooRealVar charge ("charge","",xMin,xMax);
	RooArgList vars (charge);
	double x_qmin = 0, x_qmax = 0;
	itHisto->second->Scale(1./itHisto->second->Integral());
	itHisto->second->GetQuantiles(1,&x_qmin,&fit_quantiles.first);
	itHisto->second->GetQuantiles(1,&x_qmax,&fit_quantiles.second);
	charge.setRange("fitRange",x_qmin,x_qmax);
	
	// covert histogram in a RooDataHist                                                                                                                                                        
	RooDataHist chargeDistribution (itHisto->second->GetName(),"",vars,itHisto->second);

	// Build a Landau Conv Novo PDF                                                                                                                                                         
	RooRealVar ml ("ml","",itHisto->second->GetMean(),1.,itHisto->second->GetMean()*10);
	RooRealVar sl ("sl","",itHisto->second->GetRMS(),1.,itHisto->second->GetRMS()*10);
	RooLandau  landau ("landau","",charge,ml,sl);
	RooRealVar mg ("mg","",0.,-50,50);
	RooRealVar sg ("sg","",10,0.01,50);
	RooRealVar ag ("ag","",0.1,-10,10);
	RooNovosibirsk novo ("novo","",charge,mg,sg,ag);
	RooAbsPdf* total_pdf = new RooFFTConvPdf("landauNovosibirsk","",charge,landau,novo);

	// Perform the fit
	RooFitResult* fitResult = total_pdf->fitTo(chargeDistribution,RooFit::Save(kTRUE),RooFit::SumW2Error(kTRUE),RooFit::Offset(kTRUE));
	fitResult = total_pdf->fitTo(chargeDistribution,RooFit::Range("fitRange"),RooFit::Save(kTRUE),RooFit::SumW2Error(kTRUE),RooFit::Offset(kTRUE));

	if(fitResult->status() != 0 or fitResult->covQual() <= 0){ 
	  iBadChannelFit++;	    
	  // as error use the error on the mean as a proxy
	  channelMapMPV[iMap.first]->SetBinContent(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetBinCenter(itHisto->second->GetMaximumBin()));
	  channelMapMPV[iMap.first]->SetBinError(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMeanError());	  
	}
	else{

	  RooPlot* rooframe = charge.frame();
	  chargeDistribution.plotOn(rooframe,RooFit::Name("data"));
	  total_pdf->plotOn(rooframe,RooFit::Name("total_pdf"),RooFit::Range("fitRange",kFALSE));
	  RooCurve* pdf_graph = rooframe->getCurve("total_pdf");
	  TF1* pdf_func = new TF1("pdf_func",[&](double*x, double *p){ return pdf_graph->Eval(x[0]);},x_qmin,x_qmax,0);
	  	  
 	  channelMapMPV[iMap.first]->SetBinContent(channelMapMPV[iMap.first]->FindBin(itHisto->first),pdf_func->Mean(x_qmin,x_qmax));
	  channelMapMPV[iMap.first]->SetBinError(channelMapMPV[iMap.first]->FindBin(itHisto->first),itHisto->second->GetMeanError());	  
	  
	  if(storeFitOfClusterShape and iChannel %20 == 0){
	    if(not outfile->GetDirectory(Form("Delay_%d",itHisto->first)))
	      outfile->mkdir(Form("Delay_%d",itHisto->first));
	    outfile->cd(Form("Delay_%d",itHisto->first));
	    
	    // plot to store in a root file                                                                                                                                 
	    canvas->cd();
	    
	    TH1* frame = (TH1*) itHisto->second->Clone(Form("frame_%s",itHisto->second->GetName()));
	    frame->Reset();
	    frame->GetYaxis()->SetTitle("a.u.");
	    if(TString(observable).Contains("maxCharge"))
	      frame->GetXaxis()->SetTitle("Leading strip charge (ADC)");
	    else if(observable == "clCorrectedSignalOverNoise")
	      frame->GetXaxis()->SetTitle("Cluster S/N");
	    else if(observable == "clCorrectedCharge")
	      frame->GetXaxis()->SetTitle("Cluster charge (ADC)");
	    
	    frame->GetXaxis()->SetTitleOffset(1.1);
	    frame->GetYaxis()->SetTitleOffset(1.1);
	    frame->GetYaxis()->SetTitleSize(0.040);
	    frame->GetXaxis()->SetTitleSize(0.040);
	    frame->GetYaxis()->SetLabelSize(0.035);
	    frame->GetXaxis()->SetLabelSize(0.035);
	    frame->Draw();

	    itHisto->second->SetMarkerColor(kBlack);
	    itHisto->second->SetMarkerSize(1);
	    itHisto->second->SetMarkerStyle(20);
	    itHisto->second->SetLineColor(kBlack);

	    pdf_func->SetLineColor(kRed);
	    pdf_func->SetLineWidth(2);
	    pdf_func->Draw("Csame");
	    itHisto->second->Draw("EPsame");
	    frame->GetYaxis()->SetRangeUser(0,itHisto->second->GetMaximum()*1.5);
	    CMS_lumi(canvas,"",false);
	    canvas->Modified();
	    canvas->Write(itHisto->second->GetName());
	    if(frame) delete frame;
	    if(rooframe) delete rooframe;
	    if(pdf_func) delete pdf_func;
	  }
	}
	if(fitResult) delete fitResult;
	if(total_pdf) delete total_pdf;
      }      
    }    
  }
  std::cout<<std::endl;
  if(doFitOfClusterShape)
    std::cout<<"Bad channel fit from Landau+Gaus fit "<<iBadChannelFit<<std::endl;  
  if(storeFitOfClusterShape)
    outfile->Close();

  // correct profiles and fit them with a gaussian
  std::cout<<"Analyze each single channel --> Profile "<<std::endl;
  long int iFit = 0;
  long int badFits = 0;
  long int noFitResult = 0;
  for(auto iprof : channelMapMean){
    cout.flush();
    if(iFit % 100 == 0) cout<<"\r"<<"iFit Profile "<<100*double(iFit)/(channelMapMean.size())<<" % ";
    if(iprof.second->Integral() == 0) {
      noFitResult++;
      continue;
    }
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q",false);
    if(result.Get()){
      fitStatusMean[iprof.first] = result->Status()+result->CovMatrixStatus() ;
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
    	badFits++;
      iFit++;
    }
  }
 
  std::cout<<std::endl;
  std::cout<<"NFits performed on Mean "<<iFit<<" bad ones : "<<100*double(badFits)/iFit<<" % "<<" noFit result "<<100*double(noFitResult)/iFit<<std::endl;

  std::cout<<"Start making MPV fits "<<std::endl;
  iFit = 0;
  badFits = 0;
  noFitResult = 0;
  for(auto iprof : channelMapMPV){
    cout.flush();
    if(iFit % 100 == 0) cout<<"\r"<<"iFit Profile "<<100*double(iFit)/(channelMapMPV.size())<<" % ";
    if(iprof.second->Integral() == 0) {
      noFitResult++;
      continue;
    }
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      fitStatusMPV[iprof.first] = result->Status()+result->CovMatrixStatus() ;
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
	badFits++;
      iFit++;
    }
  }
  std::cout<<std::endl;
  std::cout<<"NFits performed on MPV "<<iFit<<" bad ones : "<<100*double(badFits)/iFit<<" % "<<" noFit result "<<100*double(noFitResult)/iFit<<std::endl;
  
  std::cout<<"###################################"<<std::endl;
  std::cout<<"##### End of Channel Analysis #####"<<std::endl;
  std::cout<<"###################################"<<std::endl; 

}


void saveOutputTree(TTree* readoutMap, map<unsigned int,TH1F* > & mapHist, map<unsigned int,int> & fitStatus, const string & outputDIR, const string & postfix, const string & observable){

  TFile*  outputFile;
  TTree*  outputTree;

  std::cout<<"### Dumpt output tree for tkCommissioner "<<std::endl;
  outputFile = new TFile((outputDIR+"/tree_delayCorrection_"+postfix+".root").c_str(),"RECREATE");
  outputFile->cd();
  // branches definition for the output tree
  readoutMap->SetBranchStatus("*",kTRUE);
  readoutMap->SetBranchStatus("moduleName",kFALSE);
  readoutMap->SetBranchStatus("moduleId",kFALSE);
  readoutMap->SetBranchStatus("delaymap",kFALSE);  
  outputTree = (TTree*) readoutMap->CopyTree("");
  outputTree->SetName("delayCorrection");
  
  // set branch
  float    measuredDelay, measuredDelayUnc;
  float    measuredMeanAmplitude, measuredMeanAmplitudeUnc;
  float    measuredSigma, measuredSigmaUnc;
  float    delayCorr;
  int      fitRejected;
  int      nFilledBinRejected;
  int      sigmaRejected;
  int      amplitudeRejected;
  int      significanceRejected;
  int      notFound;
  int      peakOutBoundaryRejected;
  vector<float> amplitude;        
  vector<float> amplitudeUnc;        
  unsigned int detid;
  
  outputFile->cd();
  outputTree->SetBranchStatus("*",kTRUE);
  outputTree->GetBranch("detid")->SetName("Detid");
  outputTree->GetLeaf("detid")->SetName("Detid");
  
  TBranch* bmeasuredDelay    = outputTree->Branch("measuredDelay",&measuredDelay,"measuredDelay/F");
  TBranch* bmeasuredDelayUnc = outputTree->Branch("measuredDelayUnc",&measuredDelayUnc,"measuredDelayUnc/F");
  TBranch* bmeasuredMeanAmplitude    = outputTree->Branch("measuredMeanAmplitude",&measuredMeanAmplitude,"measuredMeanAmplitude/F");
  TBranch* bmeasuredMeanAmplitudeUnc = outputTree->Branch("measuredMeanAmplitudeUnc",&measuredMeanAmplitudeUnc,"measuredMeanAmplitudeUnc/F");
  TBranch* bmeasuredSigma    = outputTree->Branch("measuredSigma",&measuredSigma,"measuredSigma/F");
  TBranch* bmeasuredSigmaUnc = outputTree->Branch("measuredSigmaUnc",&measuredSigmaUnc,"measuredSigmaUnc/F");
  TBranch* bamplitude        = outputTree->Branch("amplitude","std::vector<float>",&amplitude);
  TBranch* bamplitudeUnc     = outputTree->Branch("amplitudeUnc","std::vector<float>",&amplitudeUnc);
  TBranch* bnotFound             = outputTree->Branch("notFound",&notFound,"notFound/I");
  TBranch* bsignificanceRejected = outputTree->Branch("significanceRejected",&significanceRejected,"significanceRejected/I");
  TBranch* bamplitudeRejected    = outputTree->Branch("amplitudeRejected",&amplitudeRejected,"amplitudeRejected/I");
  TBranch* bsigmaRejected        = outputTree->Branch("sigmaRejected",&sigmaRejected,"sigmaRejected/I");
  TBranch* bnFilledBinRejected   = outputTree->Branch("nFilledBinRejected",&nFilledBinRejected,"nFilledBinRejected/I");
  TBranch* bpeakOutBoundaryRejected   = outputTree->Branch("peakOutBoundaryRejected",&peakOutBoundaryRejected,"peakOutBoundaryRejected/I");
  TBranch* bfitRejected          = outputTree->Branch("fitRejected",&fitRejected,"fitRejected/I");
  TBranch* bdelayCorr        = outputTree->Branch("delayCorr",&delayCorr,"delayCorr/F");
  
  readoutMap->SetBranchAddress("detid",&detid);
  outputTree->SetBranchAddress("Detid",&detid);
  
  long int notFoundChannels = 0;
  
  cout<<"### Start the loop on the channel map "<<endl;      
  for(long int iChannel = 0; iChannel < readoutMap->GetEntries(); iChannel++){
    cout.flush();
    if(iChannel % 100 == 0) cout<<"\r"<<"iChannel "<<100*double(iChannel)/(readoutMap->GetEntries())<<" % ";
    readoutMap->GetEntry(iChannel);
        
    measuredDelay    = 0.;
    measuredDelayUnc = 0.;
    measuredMeanAmplitude    = 0.;
    measuredMeanAmplitudeUnc = 0.;
    measuredSigma    = 0.;
    measuredSigmaUnc = 0.;
    delayCorr        = 0.;
    notFound             = 0;
    significanceRejected = 0;
    amplitudeRejected    = 0;
    sigmaRejected        = 0;
    nFilledBinRejected   = 0;
    fitRejected          = 0;
    peakOutBoundaryRejected = 0;
    amplitude.clear();
    amplitudeUnc.clear();

    if(mapHist[detid] == 0 or mapHist[detid] == NULL or mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName())) == 0 or mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName())) == NULL){
      notFound = 1;       
      notFoundChannels++;	
    }    
    else{
      // fill the amplitude branch
      for(int iBin = 0; iBin < mapHist[detid]->GetNbinsX(); iBin++){
	amplitude.push_back(mapHist[detid]->GetBinContent(iBin));
	amplitudeUnc.push_back(mapHist[detid]->GetBinError(iBin));
      }
      
      measuredDelay    = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1);
      measuredDelayUnc = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParError(1);
      measuredMeanAmplitude    = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(0);
      measuredMeanAmplitudeUnc = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParError(0);
      measuredSigma    = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(2);
      measuredSigmaUnc = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParError(2);	
      
      /// apply selections for good fit --> check number of filled bins
      int nFiledBins = getFilledBins(mapHist[detid]);
      if(nFiledBins < minFilledBin){	
	nFilledBinRejected = 1;
	//store it in the output canvas file
	if(not outputFile->GetDirectory("binEntries"))
	  outputFile->mkdir("binEntries");
	outputFile->cd("binEntries");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  TCanvas* c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	  delete c1b;
	}
	outputFile->cd();	
      }
      
      // fit status: 0 for the fit, 3 for the covariance matrix
      if(fitStatus[detid] != 3){	  
	fitRejected = 1;
	//store it in the output canvas file
	if(not outputFile->GetDirectory("fitStatus"))
	  outputFile->mkdir("fitStatus");
	outputFile->cd("fitStatus");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  TCanvas* c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	  delete c1b;
	}
	outputFile->cd();	
      }
      
      // mean value outside boundary
      if(fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1)) > peakBoundary){	    
	peakOutBoundaryRejected = 1;	    
	//store it in the output canvas file
	if(not outputFile->GetDirectory("peakOutRange"))
	  outputFile->mkdir("peakOutRange");
	outputFile->cd("peakOutRange");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  TCanvas* c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	  delete c1b;
	}
	outputFile->cd();	
      } 
      
      // amplitude cut
      if(fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(0)) < amplitudeMin){	    
	amplitudeRejected = 1;	    
	//store it in the output canvas file
	if(not outputFile->GetDirectory("smallAmplitude"))
	  outputFile->mkdir("smallAmplitude");
	outputFile->cd("smallAmplitude");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  TCanvas* c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	  delete c1b;
	}
	outputFile->cd();	
      } 
      
      // amplitude cut
      if(fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1))/fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParError(1)) < significance and fabs(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1)) > minDelayForSignificance){	    
	significanceRejected = 1;	    
	//store it in the output canvas file
	if(not outputFile->GetDirectory("nonSignificantLargeShift"))
	  outputFile->mkdir("nonSignificantLargeShift");
	outputFile->cd("nonSignificantLargeShift");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  TCanvas* c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	  delete c1b;
	}
	outputFile->cd();	
      } 
      
      // sigma cut
      if(mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(2) > maxSigma){	    
	sigmaRejected = 1;
	//store it in the output canvas file                                                                                                                               
	if(not outputFile->GetDirectory("largeGaussSigma"))
	  outputFile->mkdir("largeGaussSigma");
	outputFile->cd("largeGaussSigma");
	if(gDirectory->GetListOfKeys()->FindObject(Form("detid_%d",detid)) == 0){
	  TCanvas* c1b = prepareCanvas(Form("detid_%d",detid),observable);
	  plotAll(c1b,mapHist[detid]);
	  c1b->Write();
	  delete c1b;
	}
	outputFile->cd();
      }
      
      if(not notFound and not significanceRejected and not amplitudeRejected and not sigmaRejected and not nFilledBinRejected and not fitRejected and not peakOutBoundaryRejected)	
	delayCorr = mapHist[detid]->GetFunction(Form("Gaus_%s",mapHist[detid]->GetName()))->GetParameter(1);
    }
    

    bmeasuredDelay->Fill();
    bmeasuredDelayUnc->Fill();
    bmeasuredMeanAmplitude->Fill();
    bmeasuredMeanAmplitudeUnc->Fill();
    bmeasuredSigma->Fill();
    bmeasuredSigmaUnc->Fill();
    bamplitude->Fill();      
    bamplitudeUnc->Fill();      
    bnotFound->Fill();
    bsignificanceRejected->Fill();
    bamplitudeRejected->Fill();
    bsigmaRejected->Fill();
    bnFilledBinRejected->Fill();
    bfitRejected->Fill();      
    bpeakOutBoundaryRejected->Fill();
    bdelayCorr->Fill();    
  }    

  cout<<endl;
  cout<<"### Write output tree for tkCommissioner: not found channels i.e. no clusters "<<100*float(notFoundChannels)/outputTree->GetEntries()<<" % "<<endl;    
  outputFile->cd();
  outputTree->BuildIndex("Detid");
  outputTree->Write(outputTree->GetName(),TObject::kOverwrite);
}


// by dedault some fit results are store in a root file, a text dump: detid fitted delay is produced, as well as deid delay uncertainty
void delayValidationPerModule(string inputDIR, 
			      string file1        = "../data/nocorrection.root", // no correction file
			      string postfix      = "merged", // sub string to be used to find merged files
			      string observable   = "maxCharge",
			      string outputDIR    = "/home/rgerosa/TrackerDAQ/DELAY_SCAN/SkimmedTrees/DelayPerModule_test_maxCharge/",
			      bool   saveMeanCanvas     =  false,
			      bool   saveMPVCanvas      =  false,
			      bool   saveCorrectionTree =  true){

  // prepare style and load macros
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());

  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming
  std::cout<<"######################################"<<std::endl;
  std::cout<<"###### delayValidationPerModule ######"<<std::endl;
  std::cout<<"######################################"<<std::endl;


  std::cout<<"### Make input file list"<<std::endl;
  system(("find "+inputDIR+" -name \"*"+postfix+"*.root\" > file.temp").c_str());
  std::ifstream infile;
  string line;
  vector<string> fileNames;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line != "" and TString(line).Contains(".root"))
        fileNames.push_back(line);
    }
  }  
  system("rm file.temp");

  std::sort(fileNames.begin(),fileNames.end());

  std::cout<<"### Build the input TTree Vector"<<std::endl;
  std::vector<long int> nEvents;
  for(auto ifile : fileNames){
    auto file = TFile::Open(ifile.c_str(),"READ");
    auto tree = (TTree*) file->FindObjectAny("clusters");
    nEvents.push_back(tree->GetEntries());
    cout<<"file "<<ifile<<" nevents "<<nEvents.back()<<endl;
  }

  TFile* _file1 = TFile::Open(file1.c_str());
  TTree* delayCorrections = (TTree*)_file1->FindObjectAny("delayCorrections");
  
  //map with key the detId number, one profile associated to it
  std::map<unsigned int,std::map<int,TH1F* > > channelMap; // for each delay value (i.e. run) gram for each detid storing the observable distribution
  std::map<unsigned int,TH1F* > channelMapMean;
  std::map<unsigned int,TH1F* > channelMapMPV;
  std::map<unsigned int,int > fitStatusMean;
  std::map<unsigned int,int > fitStatusMPV;
  ChannnelPlots(fileNames,nEvents,delayCorrections,channelMap,channelMapMean,channelMapMPV,fitStatusMean,fitStatusMPV,observable,outputDIR);  
  
  if(saveMeanCanvas)
    saveAll(channelMapMean,outputDIR,observable,"Mean"); // save all the plots into a single root files
  if(saveMPVCanvas)
    saveAll(channelMapMPV,outputDIR,observable,"MPV");
  
  //// outptut tree with analysis result
  if(saveCorrectionTree){
    TFile* file0 = TFile::Open(fileNames.front().c_str(),"READ");
    TTree* readoutMap = (TTree*) file0->FindObjectAny("readoutMap");
    saveOutputTree(readoutMap,channelMapMean,fitStatusMean,outputDIR,"Mean",observable);
    saveOutputTree(readoutMap,channelMapMPV,fitStatusMPV,outputDIR,"MPV",observable);
    file0->Close();
  }
  else{
    cout<<"### Dump peak in a text file for Mean"<<endl;
    ofstream dumpMean((outputDIR+"/dumpBestDelay_Mean_"+observable+".txt").c_str());
    int nullPointers = 0;
    for(auto imap : channelMapMean){    
      if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
	dumpMean<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParameter(1)<<"\n";
      else
	nullPointers++;
    }
    dumpMean.close();
    std::cout<<"### Null pointers "<<100*float(nullPointers)/channelMap.size()<<" %s "<<endl;

    cout<<"### Dump peak error in a text file for Mean"<<endl;
    ofstream dumpMeanError((outputDIR+"/dumpBestDelayError_Mean_"+observable+".txt").c_str());
    for(auto imap : channelMapMean){
      if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
	dumpMeanError<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParError(1)<<"\n";
    }  
    dumpMeanError.close();

    cout<<"### Dump peak in a text file for MPV"<<endl;
    ofstream dumpMPV((outputDIR+"/dumpBestDelay_MPV_"+observable+".txt").c_str());
    nullPointers = 0;
    for(auto imap : channelMapMPV){    
      if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
	dumpMPV<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParameter(1)<<"\n";
      else
	nullPointers++;
    }
    dumpMPV.close();
    std::cout<<"### Null pointers "<<100*float(nullPointers)/channelMap.size()<<" %s "<<endl;

    cout<<"### Dump peak error in a text file for MPV"<<endl;
    ofstream dumpMPVError((outputDIR+"/dumpBestDelayError_MPV_"+observable+".txt").c_str());
    for(auto imap : channelMapMPV){
      if(imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName())))
	dumpMPVError<<imap.first<<"   "<<imap.second->GetFunction(Form("Gaus_%s",imap.second->GetName()))->GetParError(1)<<"\n";
    }  
    dumpMPVError.close();
  }
  cout<<"### Clear and Close "<<endl;  
  channelMap.clear();
  channelMapMPV.clear();
  channelMapMean.clear();  
}

