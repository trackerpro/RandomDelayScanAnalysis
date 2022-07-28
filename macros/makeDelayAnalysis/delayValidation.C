#include <vector>
#include <iostream>
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
#include "../delayUtils.h"

using namespace std; 

// parametrize profile with a gaussian shape
static bool isGaussian = true;
// reduce the number of events by
static int  reductionFactor = 1;
// output files and canvases
static TFile*   outputFitFile = NULL;
static TCanvas* outputCanvasFit = NULL;
// to gauge the mean uncertainty
static int  ntoys = 100;

int makeClusterDistributionFit(TH1F* histoToFit, 
			       TH1F* histoToFill, 
			       const string & subdetector, 
			       const float & delay, 
			       const string & observable, 
			       const string & outputDIR, 
			       const string & postfix = "", 		      
			       const double & x_quantile_low  = 0.001, 
			       const double & x_quantile_high = 0.995){
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001); 
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","migradimproved");
 
  float xMin   = 0, xMax = 0;
  int   nBinsX = 0;
  setLimitsAndBinning(observable,xMin,xMax,nBinsX);

  if(outputCanvasFit == NULL or outputCanvasFit == 0)
    outputCanvasFit = new TCanvas("canvas","",600,625);
  if(outputFitFile == NULL or outputFitFile == 0)
    outputFitFile = new TFile((outputDIR+"/outputFitCanvases_"+observable+"_"+postfix+".root").c_str(),"RECREATE");
  outputFitFile->cd();

  // prepare the histogram dividing by the bin-width
  histoToFit->Scale(1./histoToFit->Integral());

  // Define observable and histogram
  RooRealVar charge ("charge","",xMin,xMax);
  RooArgList vars (charge);
  double x_qmin = 0, x_qmax = 0;
  histoToFit->GetQuantiles(1,&x_qmin,&x_quantile_low);
  histoToFit->GetQuantiles(1,&x_qmax,&x_quantile_high);   
  charge.setRange("fitRange",x_qmin,x_qmax);

  RooDataHist chargeDistribution (histoToFit->GetName(),"",vars,histoToFit);  

  // Build a Landau Conv CB PDF                                                                                                                                                                     
  RooRealVar ml ("ml","",histoToFit->GetMean(),1.,histoToFit->GetMean()*10);
  RooRealVar sl ("sl","",histoToFit->GetRMS(),1.,histoToFit->GetRMS()*10);
  RooLandau  landau ("landau","",charge,ml,sl);  
  RooRealVar mg ("mg","",0.,-50,50);
  RooRealVar sg ("sg","",10,0.01,50);
  RooRealVar ag ("ag","",0.1,-10,10);
  RooNovosibirsk novo ("novo","",charge,mg,sg,ag);
  RooAbsPdf* total_pdf = new RooFFTConvPdf("landauNovosibirsk","",charge,landau,novo);
  
  // Perform the FIT                                                                                                                                                                                 
  RooFitResult* fitResult = total_pdf->fitTo(chargeDistribution,RooFit::Save(kTRUE),RooFit::SumW2Error(kTRUE),RooFit::Offset(kTRUE));
  fitResult = total_pdf->fitTo(chargeDistribution,RooFit::Range("fitRange"),RooFit::Save(kTRUE),RooFit::SumW2Error(kTRUE),RooFit::Offset(kTRUE));

  // Increase number of bins and get postfit shape  
  int status = fitResult->status();  
  RooArgList pars = fitResult->floatParsFinal();

  vector<float> maximum_values;
  if(fitResult->status() != 0 or fitResult->covQual() <= 0){
    histoToFill->SetBinContent(histoToFill->FindBin(delay),histoToFit->GetBinCenter(histoToFit->GetMaximumBin()));
    histoToFill->SetBinError(histoToFill->FindBin(delay),histoToFit->GetMeanError());
    status = -1;
  }
  else{    

    RooPlot* rooframe = charge.frame();
    chargeDistribution.plotOn(rooframe,RooFit::Name("data"));
    total_pdf->plotOn(rooframe,RooFit::Name("total_pdf"),RooFit::Range("fitRange",kFALSE));
    RooCurve* pdf_graph = rooframe->getCurve("total_pdf");
    TF1* pdf_func = new TF1("pdf_func",[&](double*x, double *p){ return pdf_graph->Eval(x[0]);},x_qmin,x_qmax,0);

    histoToFill->SetBinContent(histoToFill->FindBin(delay),pdf_func->Mean(x_qmin,x_qmax));    
    histoToFill->SetBinError(histoToFill->FindBin(delay),histoToFit->GetMeanError());

    // plot to store in a root file                          
    outputCanvasFit->cd();
    TH1* frame = (TH1*) histoToFit->Clone(Form("frame_%s",histoToFit->GetName()));
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
    
    histoToFit->SetMarkerColor(kBlack);
    histoToFit->SetMarkerSize(1);
    histoToFit->SetMarkerStyle(20);
    histoToFit->SetLineColor(kBlack);
    
    pdf_func->SetLineColor(kRed);
    pdf_func->SetLineWidth(2);
    pdf_func->Draw("Csame");
    histoToFit->Draw("EPsame");
    frame->GetYaxis()->SetRangeUser(0,histoToFit->GetMaximum()*1.5);
    CMS_lumi(outputCanvasFit,"",false);
    
    outputCanvasFit->Modified();
    outputCanvasFit->Write(histoToFit->GetName());

    if(frame) delete frame;
    if(rooframe) delete rooframe;
    if(pdf_func) delete pdf_func;
  }
  
  if(fitResult) delete fitResult;  
  if(total_pdf) delete total_pdf;
  return status;
}

static std::map<unsigned int,std::map<int,TH1F* > > TIBlayers; // map layer:delay:distribution
static std::map<unsigned int,std::map<int,TH1F* > > TOBlayers;
static std::map<unsigned int,std::map<int,TH1F* > > TIDlayers;
static std::map<unsigned int,std::map<int,TH1F* > > TECPTlayers;
static std::map<unsigned int,std::map<int,TH1F* > > TECPtlayers;
static std::map<unsigned int,std::map<int,TH1F* > > TECMTlayers;
static std::map<unsigned int,std::map<int,TH1F* > > TECMtlayers;

static std::map<unsigned int,TH1F* > TIBlayersMean; // map delay:distribution
static std::map<unsigned int,TH1F* > TOBlayersMean;
static std::map<unsigned int,TH1F* > TIDlayersMean;
static std::map<unsigned int,TH1F* > TECPTlayersMean;
static std::map<unsigned int,TH1F* > TECPtlayersMean;
static std::map<unsigned int,TH1F* > TECMTlayersMean;
static std::map<unsigned int,TH1F* > TECMtlayersMean;

static std::map<unsigned int,TH1F* > TIBlayersMPV; // map delay:distribution
static std::map<unsigned int,TH1F* > TOBlayersMPV;
static std::map<unsigned int,TH1F* > TIDlayersMPV;
static std::map<unsigned int,TH1F* > TECPTlayersMPV;
static std::map<unsigned int,TH1F* > TECPtlayersMPV;
static std::map<unsigned int,TH1F* > TECMTlayersMPV;
static std::map<unsigned int,TH1F* > TECMtlayersMPV;



/// function that runs on the evnet and produce profiles for layers
void LayerPlots(const std::vector<std::string> & fileNames,
		const std::vector<unsigned int> & nEvents,
		TTree* corrections,
		const string & observable,
		const string & outputDIR) {

  std::cout<<"######################################################"<<std::endl;
  std::cout<< "Preparing Layer plot for the different subdetectors "<< std::endl; 
  std::cout<<"######################################################"<<std::endl;

  int   correction;
  corrections->SetBranchStatus("*",kFALSE);
  corrections->SetBranchStatus("detid",kTRUE);
  corrections->SetBranchStatus("correction",kTRUE);
  corrections->SetBranchAddress("correction",&correction);

  // create vectors for the different Profiles
  float yMin   = 0, yMax = 0;
  int   nBinsY = 0;
  vector<double> delayBins;
  setLimitsAndBinning("delay",delayBins);
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);

  for(size_t ifile = 0; ifile < fileNames.size(); ifile++){
    std::cout<<"Looping on clusters of file "<<ifile<<" out of "<<fileNames.size()<<std::endl;
    TFile* file = TFile::Open(fileNames.at(ifile).c_str(),"READ");
    TTree* tree = (TTree*) file->FindObjectAny("clusters");

    unsigned int detid;
    float  clglobalX, clglobalY, clglobalZ, thickness, obs, delay;
    tree->SetBranchStatus("*",kFALSE);
    tree->SetBranchStatus("detid",kTRUE);
    tree->SetBranchStatus("clSignalOverNoise",kTRUE);
    tree->SetBranchStatus("clglobalX",kTRUE);
    tree->SetBranchStatus("clglobalY",kTRUE);
    tree->SetBranchStatus("clglobalZ",kTRUE);
    tree->SetBranchStatus("thickness",kTRUE);
    tree->SetBranchStatus("delay",kTRUE);
    tree->SetBranchStatus(observable.c_str(),kTRUE);
    tree->SetBranchAddress("detid",&detid);
    tree->SetBranchAddress("clglobalX",&clglobalX);
    tree->SetBranchAddress("clglobalY",&clglobalY);
    tree->SetBranchAddress("clglobalZ",&clglobalZ);
    tree->SetBranchAddress("thickness",&thickness);
    tree->SetBranchAddress("delay",&delay);
    tree->SetBranchAddress(observable.c_str(),&obs);

    for(int iEvent = 0; iEvent < nEvents.at(ifile)/reductionFactor; iEvent++){    
      cout.flush();
      if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(nEvents.at(ifile)/reductionFactor)<<" % ";
      // take the event
      tree->GetEntry(iEvent);
      // take the correction from the detid
      corrections->GetEntryWithIndex(detid);

      unsigned int subdetid    = int((detid-0x10000000)/0x2000000);
      unsigned int barrellayer = int((detid%33554432)/0x4000);
      unsigned int TIDlayer    = int((detid%33554432)/0x800)%4;
      unsigned int TECPlayer   = int((detid%33554432)/0x4000)-32;
      unsigned int TECMlayer   = int((detid%33554432)/0x4000)-16;
      float R = sqrt(clglobalX*clglobalX+clglobalY*clglobalY+clglobalZ*clglobalZ);      
      float value = obs;
      
      int delaybin = std::round(delay)-std::round(correction);

      // fill the maps
      if(subdetid == 3){
	if(TIBlayers[barrellayer][delaybin] == 0 or TIBlayers[barrellayer][delaybin] == NULL){
	  TIBlayers[barrellayer][delaybin] = new TH1F(Form("TIB_layer_%d_delay_%d",barrellayer,delaybin),"",nBinsY,yMin,yMax);
	  TIBlayers[barrellayer][delaybin]->SetDirectory(0);
	  TIBlayers[barrellayer][delaybin]->Sumw2();
	}
	TIBlayers[barrellayer][delaybin]->Fill(value);
      }
      else if(subdetid == 5){
	if(TOBlayers[barrellayer][delaybin] == 0 or TOBlayers[barrellayer][delaybin] == NULL) {
	  TOBlayers[barrellayer][delaybin] = new TH1F(Form("TOB_layer_%d_delay_%d",barrellayer,delaybin),"",nBinsY,yMin,yMax);
	  TOBlayers[barrellayer][delaybin]->SetDirectory(0);
	  TOBlayers[barrellayer][delaybin]->Sumw2();
	}
	TOBlayers[barrellayer][delaybin]->Fill(value);
      }
      else if(subdetid == 4){
	if(TIDlayers[TIDlayer][delaybin] == 0 or TIDlayers[TIDlayer][delaybin] == NULL) {
	  TIDlayers[TIDlayer][delaybin] = new TH1F(Form("TID_layer_%d_delay_%d",TIDlayer,delaybin),"",nBinsY,yMin,yMax);
	  TIDlayers[TIDlayer][delaybin]->SetDirectory(0);
	  TIDlayers[TIDlayer][delaybin]->Sumw2();
	}
	TIDlayers[TIDlayer][delaybin]->Fill(value);
      }
      else if(subdetid == 6  and clglobalZ > 0 and thickness > 400){
	if(TECPTlayers[TECPlayer][delaybin] == 0 or TECPTlayers[TECPlayer][delaybin] == NULL){
	  TECPTlayers[TECPlayer][delaybin] = new TH1F(Form("TECPT_layer_%d_delay_%d",TECPlayer,delaybin),"",nBinsY,yMin,yMax);
	  TECPTlayers[TECPlayer][delaybin]->SetDirectory(0);
	  TECPTlayers[TECPlayer][delaybin]->Sumw2();
	}
	TECPTlayers[TECPlayer][delaybin]->Fill(value);
      }
      else if(subdetid == 6  and clglobalZ > 0 and thickness < 400){
	if(TECPtlayers[TECPlayer][delaybin] == 0 or TECPtlayers[TECPlayer][delaybin] == NULL){
	  TECPtlayers[TECPlayer][delaybin] = new TH1F(Form("TECPt_layer_%d_delay_%d",TECPlayer,delaybin),"",nBinsY,yMin,yMax);
	  TECPtlayers[TECPlayer][delaybin]->SetDirectory(0);
	  TECPtlayers[TECPlayer][delaybin]->Sumw2();
	}
	TECPtlayers[TECPlayer][delaybin]->Fill(value);
      }
      else if(subdetid == 6  and clglobalZ < 0 and thickness > 400){
	if(TECMTlayers[TECMlayer][delaybin] == 0 or TECMTlayers[TECMlayer][delaybin] == NULL){
	  TECMTlayers[TECMlayer][delaybin] = new TH1F(Form("TECMT_layer_%d_delay_%d",TECMlayer,delaybin),"",nBinsY,yMin,yMax);
	  TECMTlayers[TECMlayer][delaybin]->SetDirectory(0);
	  TECMTlayers[TECMlayer][delaybin]->Sumw2();
	}
	TECMTlayers[TECMlayer][delaybin]->Fill(value);
      }
      else if(subdetid == 6  and clglobalZ < 0 and thickness < 400){
	if(TECMtlayers[TECMlayer][delaybin] == 0 or TECMtlayers[TECMlayer][delaybin] == NULL){
	  TECMtlayers[TECMlayer][delaybin] = new TH1F(Form("TECMt_layer_%d_delay_%d",TECMlayer,delaybin),"",nBinsY,yMin,yMax);
	  TECMtlayers[TECMlayer][delaybin]->SetDirectory(0);
	  TECMtlayers[TECMlayer][delaybin]->Sumw2();
	}
	TECMtlayers[TECMlayer][delaybin]->Fill(value);
      }
    }    
    std::cout<<std::endl;
    std::cout<<"Loop on events terminated"<<std::endl;
    file->Close();
  }

  std::cout<<"Build the mean and MPV distributions asaf of delay "<<endl;  
  unsigned int iBadChannelFit = 0;
  for(auto imap : TIBlayers){ // loop on the different layer
    if(TIBlayersMean[imap.first] == 0 or TIBlayersMean[imap.first] == NULL){
      TIBlayersMean[imap.first] = new TH1F(Form("TIB_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TIBlayersMean[imap.first]->SetDirectory(0);
      TIBlayersMean[imap.first]->Sumw2();
    }
    if(TIBlayersMPV[imap.first] == 0 or TIBlayersMPV[imap.first] == NULL){
      TIBlayersMPV[imap.first] = new TH1F(Form("TIB_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TIBlayersMPV[imap.first]->SetDirectory(0);
      TIBlayersMPV[imap.first]->Sumw2();
    }
    for(auto idelay : imap.second){
      // mean value
      TIBlayersMean[imap.first]->SetBinContent(TIBlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TIBlayersMean[imap.first]->SetBinError(TIBlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());      
      int status = makeClusterDistributionFit(idelay.second,TIBlayersMPV[imap.first],"TIB",idelay.first,observable,outputDIR);
      if(status != 0) iBadChannelFit++;            	
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TIB layers "<<iBadChannelFit<<" over "<<TIBlayers.size()*TIBlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TOBlayers){ // loop on the different layer
    if(TOBlayersMean[imap.first] == 0 or TOBlayersMean[imap.first] == NULL){
      TOBlayersMean[imap.first] = new TH1F(Form("TOB_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TOBlayersMean[imap.first]->SetDirectory(0);
      TOBlayersMean[imap.first]->Sumw2();
    }
    if(TOBlayersMPV[imap.first] == 0 or TOBlayersMPV[imap.first] == NULL){
      TOBlayersMPV[imap.first] = new TH1F(Form("TOB_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TOBlayersMPV[imap.first]->SetDirectory(0);
      TOBlayersMPV[imap.first]->Sumw2();
    }
    for(auto idelay : imap.second){
      // mean value
      TOBlayersMean[imap.first]->SetBinContent(TOBlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TOBlayersMean[imap.first]->SetBinError(TOBlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TOBlayersMPV[imap.first],"TOB",idelay.first,observable,outputDIR);
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TOB layers "<<iBadChannelFit<<" over "<<TOBlayers.size()*TOBlayers[1].size()<<std::endl;  

  iBadChannelFit = 0;
  for(auto imap : TIDlayers){ // loop on the different layer
    if(TIDlayersMean[imap.first] == 0 or TIDlayersMean[imap.first] == NULL){
      TIDlayersMean[imap.first] = new TH1F(Form("TID_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TIDlayersMean[imap.first]->SetDirectory(0);
      TIDlayersMean[imap.first]->Sumw2();
    }
    if(TIDlayersMPV[imap.first] == 0 or TIDlayersMPV[imap.first] == NULL){
      TIDlayersMPV[imap.first] = new TH1F(Form("TID_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TIDlayersMPV[imap.first]->SetDirectory(0);
      TIDlayersMPV[imap.first]->Sumw2();
    }
    for(auto idelay : imap.second){
      // mean value
      TIDlayersMean[imap.first]->SetBinContent(TIDlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TIDlayersMean[imap.first]->SetBinError(TIDlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TIDlayersMPV[imap.first],"TID",idelay.first,observable,outputDIR);
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TID layers "<<iBadChannelFit<<" over "<<TIDlayers.size()*TIDlayers[1].size()<<std::endl;  
  
  iBadChannelFit = 0;
  for(auto imap : TECPTlayers){ // loop on the different layer
    if(TECPTlayersMean[imap.first] == 0 or TECPTlayersMean[imap.first] == NULL){
      TECPTlayersMean[imap.first] = new TH1F(Form("TECPT_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECPTlayersMean[imap.first]->SetDirectory(0);
      TECPTlayersMean[imap.first]->Sumw2();
    }
    if(TECPTlayersMPV[imap.first] == 0 or TECPTlayersMPV[imap.first] == NULL){
      TECPTlayersMPV[imap.first] = new TH1F(Form("TECPT_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECPTlayersMPV[imap.first]->SetDirectory(0);
      TECPTlayersMPV[imap.first]->Sumw2();
    }
    for(auto idelay : imap.second){
      // mean value
      TECPTlayersMean[imap.first]->SetBinContent(TECPTlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECPTlayersMean[imap.first]->SetBinError(TECPTlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TECPTlayersMPV[imap.first],"TECPT",idelay.first,observable,outputDIR);
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPT layers "<<iBadChannelFit<<" over "<<TECPTlayers.size()*TECPTlayers[1].size()<<std::endl;  
  
  iBadChannelFit = 0;
  for(auto imap : TECPtlayers){ // loop on the different layer
    if(TECPtlayersMean[imap.first] == 0 or TECPtlayersMean[imap.first] == NULL){
      TECPtlayersMean[imap.first] = new TH1F(Form("TECPt_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECPtlayersMean[imap.first]->SetDirectory(0);
      TECPtlayersMean[imap.first]->Sumw2();
    }
    if(TECPtlayersMPV[imap.first] == 0 or TECPtlayersMPV[imap.first] == NULL){
      TECPtlayersMPV[imap.first] = new TH1F(Form("TECPt_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECPtlayersMPV[imap.first]->SetDirectory(0);
      TECPtlayersMPV[imap.first]->Sumw2();
    }
    for(auto idelay : imap.second){
      // mean value
      TECPtlayersMean[imap.first]->SetBinContent(TECPtlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECPtlayersMean[imap.first]->SetBinError(TECPtlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TECPtlayersMPV[imap.first],"TECPt",idelay.first,observable,outputDIR);
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPt layers "<<iBadChannelFit<<" over "<<TECPtlayers.size()*TECPtlayers[1].size()<<std::endl;  
  
  iBadChannelFit = 0;
  for(auto imap : TECMTlayers){ // loop on the different layer
    if(TECMTlayersMean[imap.first] == 0 or TECMTlayersMean[imap.first] == NULL){
      TECMTlayersMean[imap.first] = new TH1F(Form("TECMT_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECMTlayersMean[imap.first]->SetDirectory(0);
      TECMTlayersMean[imap.first]->Sumw2();
    }
    if(TECMTlayersMPV[imap.first] == 0 or TECMTlayersMPV[imap.first] == NULL){
      TECMTlayersMPV[imap.first] = new TH1F(Form("TECMT_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECMTlayersMPV[imap.first]->SetDirectory(0);
      TECMTlayersMPV[imap.first]->Sumw2();
    }
    for(auto idelay : imap.second){
      // mean value
      TECMTlayersMean[imap.first]->SetBinContent(TECMTlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECMTlayersMean[imap.first]->SetBinError(TECMTlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TECMTlayersMPV[imap.first],"TECPT",idelay.first,observable,outputDIR);
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMT layers "<<iBadChannelFit<<" over "<<TECMTlayers.size()*TECMTlayers[1].size()<<std::endl;  
  
  iBadChannelFit = 0;
  for(auto imap : TECMtlayers){ // loop on the different layer
    if(TECMtlayersMean[imap.first] == 0 or TECMtlayersMean[imap.first] == NULL){
      TECMtlayersMean[imap.first] = new TH1F(Form("TECMt_layer_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECMtlayersMean[imap.first]->SetDirectory(0);
      TECMtlayersMean[imap.first]->Sumw2();
    }
    if(TECMtlayersMPV[imap.first] == 0 or TECMtlayersMPV[imap.first] == NULL){
      TECMtlayersMPV[imap.first] = new TH1F(Form("TECMt_layer_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECMtlayersMPV[imap.first]->SetDirectory(0);
      TECMtlayersMPV[imap.first]->Sumw2();
    }
    for(auto idelay : imap.second){
      // mean value
      TECMtlayersMean[imap.first]->SetBinContent(TECMtlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECMtlayersMean[imap.first]->SetBinError(TECMtlayersMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TECMtlayersMPV[imap.first],"TECMt",idelay.first,observable,outputDIR);
      if(status != 0) iBadChannelFit++;      
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMt layers "<<iBadChannelFit<<" over "<<TECMtlayers.size()*TECMtlayers[1].size()<<std::endl;  

  // correct profiles and fit them with a gaussian
  unsigned int badFits = 0;
  for(auto iprof : TIBlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get())
      if(result->CovMatrixStatus() != 3 or result->Status() != 0){
	badFits++;
    }
  }
  std::cout<<"TIB Mean bad gaussian fits "<<badFits<<std::endl;
  
  badFits = 0;
  for(auto iprof : TOBlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
	badFits++;
    }
  }
  std::cout<<"TOB Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TIDlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TID Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECPTlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECPT Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECPtlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECPt Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECMTlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECMT Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECMtlayersMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECMt Mean bad gaussian fits "<<badFits<<std::endl;

  // correct profiles and fit them with a gaussian
  badFits = 0;
  for(auto iprof : TIBlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TIB MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TOBlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TOB MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TIDlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TID MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECPTlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECPT MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECPtlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECPt MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECMTlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECMT MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECMtlayersMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECMt MPV bad gaussian fits "<<badFits<<std::endl;
  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Layer Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;
}

// create the plots in R slices 
static std::map<unsigned int, std::map<float,TH1F* > > TIBrs;
static std::map<unsigned int, std::map<float,TH1F* > > TIDrs;
static std::map<unsigned int, std::map<float,TH1F* > > TOBrs;
static std::map<unsigned int, std::map<float,TH1F* > > TECPTrs;
static std::map<unsigned int, std::map<float,TH1F* > > TECPtrs;
static std::map<unsigned int, std::map<float,TH1F* > > TECMTrs;
static std::map<unsigned int, std::map<float,TH1F* > > TECMtrs;
static std::map<unsigned int, std::map<float,TH1F* > > TIB;
static std::map<unsigned int, std::map<float,TH1F* > > TID;
static std::map<unsigned int, std::map<float,TH1F* > > TOB;
static std::map<unsigned int, std::map<float,TH1F* > > TECPT;
static std::map<unsigned int, std::map<float,TH1F* > > TECPt;
static std::map<unsigned int, std::map<float,TH1F* > > TECMT;
static std::map<unsigned int, std::map<float,TH1F* > > TECMt;

static std::map<unsigned int, TH1F* > TIBrsMean;
static std::map<unsigned int, TH1F* > TIDrsMean;
static std::map<unsigned int, TH1F* > TOBrsMean;
static std::map<unsigned int, TH1F* > TECPTrsMean;
static std::map<unsigned int, TH1F* > TECPtrsMean;
static std::map<unsigned int, TH1F* > TECMTrsMean;
static std::map<unsigned int, TH1F* > TECMtrsMean;

static std::map<unsigned int, TH1F* > TIBMean;
static std::map<unsigned int, TH1F* > TIDMean;
static std::map<unsigned int, TH1F* > TOBMean;
static std::map<unsigned int, TH1F* > TECPTMean;
static std::map<unsigned int, TH1F* > TECPtMean;
static std::map<unsigned int, TH1F* > TECMTMean;
static std::map<unsigned int, TH1F* > TECMtMean;

static std::map<unsigned int, TH1F* > TIBrsMPV;
static std::map<unsigned int, TH1F* > TIDrsMPV;
static std::map<unsigned int, TH1F* > TOBrsMPV;
static std::map<unsigned int, TH1F* > TECPTrsMPV;
static std::map<unsigned int, TH1F* > TECPtrsMPV;
static std::map<unsigned int, TH1F* > TECMTrsMPV;
static std::map<unsigned int, TH1F* > TECMtrsMPV;
static std::map<unsigned int, TH1F* > TIBMPV;
static std::map<unsigned int, TH1F* > TIDMPV;
static std::map<unsigned int, TH1F* > TOBMPV;
static std::map<unsigned int, TH1F* > TECPTMPV;
static std::map<unsigned int, TH1F* > TECPtMPV;
static std::map<unsigned int, TH1F* > TECMTMPV;
static std::map<unsigned int, TH1F* > TECMtMPV;


//// Per ring analysis
void RPlots(const std::vector<std::string> & fileNames,
	    const std::vector<unsigned int> & nEvents,
	    TTree* corrections,
	    const std::string & observable,
	    const std::string & outputDIR){



  std::cout<<"######################################################"<<std::endl;
  std::cout << "Preparing Ring plot for the different subdetectors "<< std::endl; 
  std::cout<<"######################################################"<<std::endl;

  int correction;
  corrections->SetBranchStatus("*",kFALSE);
  corrections->SetBranchStatus("detid",kTRUE);
  corrections->SetBranchStatus("correction",kTRUE);
  corrections->SetBranchAddress("correction",&correction);

  // create vectors for the different Profiles
  float yMin = 0, yMax = 0;
  int   nBinsY = 0;
  vector<double> delayBins;
  setLimitsAndBinning("delay",delayBins);
  setLimitsAndBinning(observable,yMin,yMax,nBinsY);

  for(size_t ifile = 0; ifile < fileNames.size(); ifile++){
    std::cout<<"Looping on clusters of file "<<ifile<<" out of "<<fileNames.size()<<std::endl;
    auto file = TFile::Open(fileNames.at(ifile).c_str(),"READ");
    auto tree = (TTree*) file->FindObjectAny("clusters");
    
    // TTree Reader appear not to be working with addFriend and EventList
    unsigned int detid;
    float    clglobalX, clglobalY, clglobalZ, thickness, obs, delay;
    tree->SetBranchStatus("*",kFALSE);
    tree->SetBranchStatus("detid",kTRUE);
    tree->SetBranchStatus("clglobalX",kTRUE);
    tree->SetBranchStatus("clglobalY",kTRUE);
    tree->SetBranchStatus("clglobalZ",kTRUE);
    tree->SetBranchStatus("thickness",kTRUE);
    tree->SetBranchStatus("delay",kTRUE);
    tree->SetBranchStatus(observable.c_str(),kTRUE);
    tree->SetBranchAddress("detid",&detid);
    tree->SetBranchAddress("clglobalX",&clglobalX);
    tree->SetBranchAddress("clglobalY",&clglobalY);
    tree->SetBranchAddress("clglobalZ",&clglobalZ);
    tree->SetBranchAddress("thickness",&thickness);
    tree->SetBranchAddress("delay",&delay);
    tree->SetBranchAddress(observable.c_str(),&obs);
    
    unsigned int iEvent = 0;
    for( ; iEvent < nEvents.at(ifile)/reductionFactor; iEvent++){    
      tree->GetEntry(iEvent);
      // take the correction from the detid
      corrections->GetEntryWithIndex(detid);

      cout.flush();
      if(iEvent % 100000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(nEvents.at(ifile)/reductionFactor)<<" % ";
      
      unsigned int subdetid = int((detid-0x10000000)/0x2000000);
      float    R           = sqrt(clglobalX*clglobalX+clglobalY*clglobalY+clglobalZ*clglobalZ);
      float    value       = obs;
      
      int delaybin = std::round(delay)-std::round(correction);
      
      if(subdetid == 3){
	int ring = int((R-TIBRing.rMin)/((TIBRing.rMax-TIBRing.rMin)/TIBRing.nDivision))+1;   
	if(TIBrs[ring][delaybin] == 0 or TIBrs[ring][delaybin] == NULL){
	  TIBrs[ring][delaybin] = new TH1F(Form("TIB_ring_%d_delay_%d",ring,delaybin),"",nBinsY,yMin,yMax);      
	  TIBrs[ring][delaybin]->SetDirectory(0);
	}
	TIBrs[ring][delaybin]->Fill(value);
      }    
      else if(subdetid == 5){
	int ring = int((R-TOBRing.rMin)/((TOBRing.rMax-TOBRing.rMin)/TOBRing.nDivision))+1;
	if(TOBrs[ring][delaybin] == 0 or TOBrs[ring][delaybin] == NULL){
	  TOBrs[ring][delaybin] = new TH1F(Form("TOB_ring_%d_delay_%d",ring,delaybin),"",nBinsY,yMin,yMax);      
	  TOBrs[ring][delaybin]->SetDirectory(0);
	}
	TOBrs[ring][delaybin]->Fill(value);
	
      }
      else if(subdetid == 4){
	int ring = int((R-TOBRing.rMin)/((TOBRing.rMax-TOBRing.rMin)/TOBRing.nDivision))+1;
	if(TIDrs[ring][delaybin] == 0 or TIDrs[ring][delaybin] == NULL){      
	  TIDrs[ring][delaybin] = new TH1F(Form("TID_ring_%d_delay_%d",ring,delaybin),"",nBinsY,yMin,yMax);      
	  TIDrs[ring][delaybin]->SetDirectory(0);
	}
	TIDrs[ring][delaybin]->Fill(value);	
      }
      
      else if(subdetid == 6  and clglobalZ > 0 and thickness > 400){
	int ring = int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision))+1;
	if(TECPTrs[ring][delaybin] == 0 or TECPTrs[ring][delaybin] == NULL){
	  TECPTrs[ring][delaybin] = new TH1F(Form("TECPT_ring_%d_delay_%d",ring,delaybin),"",nBinsY,yMin,yMax);      
	  TECPTrs[ring][delaybin]->SetDirectory(0);
	}
	TECPTrs[ring][delaybin]->Fill(value);
	
      }
      else if(subdetid == 6  and clglobalZ > 0 and thickness < 400){
	int ring = int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision))+1;
	if(TECPtrs[ring][delaybin] == 0 or TECPtrs[ring][delaybin] == NULL){
	  TECPtrs[ring][delaybin] = new TH1F(Form("TECPt_ring_%d_delay_%d",ring,delaybin),"",nBinsY,yMin,yMax);      
	  TECPtrs[ring][delaybin]->SetDirectory(0);
	}
	TECPtrs[ring][delaybin]->Fill(value);
	
      }
      else if(subdetid == 6  and clglobalZ < 0 and thickness > 400){
	int ring = int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision))+1;
	if(TECMTrs[ring][delaybin] == 0 or TECMTrs[ring][delaybin] == NULL){
	  TECMTrs[ring][delaybin] = new TH1F(Form("TECMT_ring_%d_delay_%d",ring,delaybin),"",nBinsY,yMin,yMax);      
	  TECMTrs[ring][delaybin]->SetDirectory(0);
	}
	TECMTrs[ring][delaybin]->Fill(value);	
      }
      else if(subdetid == 6  and clglobalZ < 0 and thickness < 400){
	int ring = int((R-TECRing.rMin)/((TECRing.rMax-TECRing.rMin)/TECRing.nDivision))+1;
	if(TECMtrs[ring][delaybin] == 0 or TECMtrs[ring][delaybin] == NULL){
	  TECMtrs[ring][delaybin] = new TH1F(Form("TECMt_ring_%d_delay_%d",ring,delaybin),"",nBinsY,yMin,yMax);      
	  TECMtrs[ring][delaybin]->SetDirectory(0);
	}
	TECMtrs[ring][delaybin]->Fill(value);
      }
    }
    std::cout<<std::endl;
    std::cout<<"Loop on events terminated"<<std::endl;
    file->Close();
  }

  std::cout<<"Build the mean and MPV distributions asaf of delay "<<endl;
  unsigned int iBadChannelFit = 0;
  for(auto imap : TIBrs){ // loop on the different layer                                                                                                                     
    if(TIBrsMean[imap.first] == 0 or TIBrsMean[imap.first] == NULL){
      TIBrsMean[imap.first] = new TH1F(Form("TIB_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TIBrsMean[imap.first]->SetDirectory(0);
    }
    if(TIBrsMPV[imap.first] == 0 or TIBrsMPV[imap.first] == NULL){
      TIBrsMPV[imap.first] = new TH1F(Form("TIB_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TIBrsMPV[imap.first]->SetDirectory(0);
    }

    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TIBrsMean[imap.first]->SetBinContent(TIBrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TIBrsMean[imap.first]->SetBinError(TIBrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TIBrsMPV[imap.first],"TIB",idelay.first,observable,outputDIR,"partitions");
      if(status != 0) iBadChannelFit++;
    }
  }

  std::cout<<"Bad channel fit from Landau+Gaus fit in TIB rings "<<iBadChannelFit<<" over "<<TIBrs.size()*TIBrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TOBrs){ // loop on the different layer                                                                                                                    
    if(TOBrsMean[imap.first] == 0 or TOBrsMean[imap.first] == NULL){
      TOBrsMean[imap.first] = new TH1F(Form("TOB_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TOBrsMean[imap.first]->SetDirectory(0);
    }
    if(TOBrsMPV[imap.first] == 0 or TOBrsMPV[imap.first] == NULL){
      TOBrsMPV[imap.first] = new TH1F(Form("TOB_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TOBrsMPV[imap.first]->SetDirectory(0);
    }
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TOBrsMean[imap.first]->SetBinContent(TOBrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TOBrsMean[imap.first]->SetBinError(TOBrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TOBrsMPV[imap.first],"TOB",idelay.first,observable,outputDIR,"partitions");
      if(status != 0) iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TOB rings "<<iBadChannelFit<<" over "<<TOBrs.size()*TOBrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TIDrs){ // loop on the different layer                                                                                                                     
    if(TIDrsMean[imap.first] == 0 or TIDrsMean[imap.first] == NULL){
      TIDrsMean[imap.first] = new TH1F(Form("TID_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TIDrsMean[imap.first]->SetDirectory(0);
    }
    if(TIDrsMPV[imap.first] == 0 or TIDrsMPV[imap.first] == NULL){
      TIDrsMPV[imap.first] = new TH1F(Form("TID_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TIDrsMPV[imap.first]->SetDirectory(0);
    }
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TIDrsMean[imap.first]->SetBinContent(TIDrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TIDrsMean[imap.first]->SetBinError(TIDrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TIDrsMPV[imap.first],"TID",idelay.first,observable,outputDIR);
      if(status != 0) iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TID rings "<<iBadChannelFit<<" over "<<TIDrs.size()*TIDrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TECPTrs){ // loop on the different layer                                                                                                                     
    if(TECPTrsMean[imap.first] == 0 or TECPTrsMean[imap.first] == NULL){
      TECPTrsMean[imap.first] = new TH1F(Form("TECPT_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECPTrsMean[imap.first]->SetDirectory(0);
    }
    if(TECPTrsMPV[imap.first] == 0 or TECPTrsMPV[imap.first] == NULL){
      TECPTrsMPV[imap.first] = new TH1F(Form("TECPT_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECPTrsMPV[imap.first]->SetDirectory(0);
    }
    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TECPTrsMean[imap.first]->SetBinContent(TECPTrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECPTrsMean[imap.first]->SetBinError(TECPTrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TECPTrsMPV[imap.first],"TECPT",idelay.first,observable,outputDIR,"partitions");
      if(status != 0) iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPT rings "<<iBadChannelFit<<" over "<<TECPTrs.size()*TECPTrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TECPtrs){ // loop on the different layer                                                                                                                     
    if(TECPtrsMean[imap.first] == 0 or TECPtrsMean[imap.first] == NULL){
      TECPtrsMean[imap.first] = new TH1F(Form("TECPt_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECPtrsMean[imap.first]->SetDirectory(0);
    }
    if(TECPtrsMPV[imap.first] == 0 or TECPtrsMPV[imap.first] == NULL){
      TECPtrsMPV[imap.first] = new TH1F(Form("TECPt_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECPtrsMPV[imap.first]->SetDirectory(0);
    }
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TECPtrsMean[imap.first]->SetBinContent(TECPtrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECPtrsMean[imap.first]->SetBinError(TECPtrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TECPtrsMPV[imap.first],"TECPt",idelay.first,observable,outputDIR,"partitions");
      if(status != 0) iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECPt rings "<<iBadChannelFit<<" over "<<TECPtrs.size()*TECPtrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TECMTrs){ // loop on the different layer                                                                                                                     
    if(TECMTrsMean[imap.first] == 0 or TECMTrsMean[imap.first] == NULL){
      TECMTrsMean[imap.first] = new TH1F(Form("TECMT_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECMTrsMean[imap.first]->SetDirectory(0);
    }
    if(TECMTrsMPV[imap.first] == 0 or TECMTrsMPV[imap.first] == NULL){
      TECMTrsMPV[imap.first] = new TH1F(Form("TECMT_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECMTrsMPV[imap.first]->SetDirectory(0);
    }    
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TECMTrsMean[imap.first]->SetBinContent(TECMTrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECMTrsMean[imap.first]->SetBinError(TECMTrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TECMTrsMPV[imap.first],"TECMT",idelay.first,observable,outputDIR,"partitions");
      if(status != 0) iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMT rings "<<iBadChannelFit<<" over "<<TECMTrs.size()*TECMTrs[1].size()<<std::endl;

  iBadChannelFit = 0;
  for(auto imap : TECMtrs){ // loop on the different layer                                                                                                                     
    if(TECMtrsMean[imap.first] == 0 or TECMtrsMean[imap.first] == NULL){
      TECMtrsMean[imap.first] = new TH1F(Form("TECMt_ring_%d_mean",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECMtrsMean[imap.first]->SetDirectory(0);
    }
    if(TECMtrsMPV[imap.first] == 0 or TECMtrsMPV[imap.first] == NULL){
      TECMtrsMPV[imap.first] = new TH1F(Form("TECMt_ring_%d_mpv",imap.first),"",delayBins.size()-1,&delayBins[0]);
      TECMtrsMPV[imap.first]->SetDirectory(0);
    }
    for(auto idelay : imap.second){
      // mean value                                                                                                                                                             
      TECMtrsMean[imap.first]->SetBinContent(TECMtrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMean());
      TECMtrsMean[imap.first]->SetBinError(TECMtrsMean[imap.first]->FindBin(idelay.first),idelay.second->GetMeanError());
      int status = makeClusterDistributionFit(idelay.second,TECMtrsMPV[imap.first],"TECMt",idelay.first,observable,outputDIR,"partitions");
      if(status != 0) iBadChannelFit++;
    }
  }
  std::cout<<"Bad channel fit from Landau+Gaus fit in TECMt rings "<<iBadChannelFit<<" over "<<TECMtrs.size()*TECMtrs[1].size()<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Loop on events terminated"<<std::endl;

  // correct profiles and fit them with a gaussian
  unsigned int badFits = 0;
  for(auto iprof : TIBrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TIB Mean bad gaussian fits "<<badFits<<std::endl;
  badFits = 0;
  for(auto iprof : TOBrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TOB Mean bad gaussian fits "<<badFits<<std::endl;
  badFits = 0;
  for(auto iprof : TIDrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TID Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECPTrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECPT Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECPtrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECPt Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECMTrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECMT Mean bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECMtrsMean){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECMt Mean bad gaussian fits "<<badFits<<std::endl;

  // correct profiles and fit them with a gaussian
  badFits = 0;
  for(auto iprof : TIBrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TIB MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TOBrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TOB MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"Analyze TID profiles mpv layers"<<std::endl;
  badFits = 0;
  for(auto iprof : TIDrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TID MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECPTrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECPT MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECPtrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECPt MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECMTrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECMT MPV bad gaussian fits "<<badFits<<std::endl;

  badFits = 0;
  for(auto iprof : TECMtrsMPV){
    TFitResultPtr result = fitHistogram(iprof.second,isGaussian,"Q");
    if(result.Get()){
      if(result->CovMatrixStatus() != 3 or result->Status() != 0)
        badFits++;
    }
  }
  std::cout<<"TECMt MPV bad gaussian fits "<<badFits<<std::endl;

  std::cout<<"#################################"<<std::endl;
  std::cout<<"##### End of Ring Analysis #####"<<std::endl;
  std::cout<<"#################################"<<std::endl;
}


/// main function that run the analysis
void delayValidation(string inputDIR,  // inputDIR with all files, use *.root to select them
		     string fileCorr = "../../data/nocorrection.root",  // possible file with correction
		     string nameToGrep   = "",
		     string observable   = "maxChargeCorrected",   // observable to be considered: maxCharge, S/N ..etc
		     bool plotPartitions = true, // best delay setting in each partition 
		     bool plotLayer      = true, // best delay setting per layer
		     bool plotSlices     = false, // best delay setting per ring
		     string outputDIR    = "prompt" // output directory name
		     ){

  // prepare style and load macros
  setTDRStyle();
  // not dump stat and fit info
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gROOT->SetBatch(kTRUE);

  system(("mkdir -p "+outputDIR).c_str());

  // open the file and prepare the cluster tree, adding the other trees as frined --> memory consuming
  std::cout<<"#############################"<<std::endl;
  std::cout<<"###### delayValidation ######"<<std::endl;
  std::cout<<"#############################"<<std::endl;

  std::cout<<"### Make input file list"<<std::endl;
  // file names and events
  std::vector<string> fileNames;
  system(("find "+inputDIR+" -name \"*"+nameToGrep+"*.root\" > file.temp").c_str());
  std::ifstream infile;
  string line;
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
  std::vector<unsigned int> nEvents;
  for(auto ifile : fileNames){
    auto file = TFile::Open(ifile.c_str(),"READ");
    auto tree = (TTree*) file->FindObjectAny("clusters");
    nEvents.push_back(tree->GetEntries());
    cout<<"file "<<ifile<<" nevents "<<nEvents.back()<<endl;
  }

  TFile* _file1 = TFile::Open(fileCorr.c_str());
  TTree* delayCorrections = (TTree*) _file1->FindObjectAny("delayCorrections");

  // create the plots per layer
  if(plotLayer){
   
    // run per layer analysis
    LayerPlots(fileNames,nEvents,delayCorrections,observable,outputDIR);

    // canvas for layer one mean observable result
    TCanvas* c1_mean = prepareCanvas("TIB_layers_mean",observable);
    plotAll(c1_mean,TIBlayersMean);
    c1_mean->Print(Form("%s/TIB_layers_mean.root",outputDIR.c_str()));
    // canvas for layer one MPV observable result 
    TCanvas* c1_mpv = prepareCanvas("TIB_layers_mpv",observable);
    plotAll(c1_mpv,TIBlayersMPV);
    c1_mpv->Print(Form("%s/TIB_layers_mpv.root",outputDIR.c_str()));
    
    /// Layer TID
    TCanvas* c2_mean = prepareCanvas("TID_layers_mean",observable);
    plotAll(c2_mean,TIDlayersMean);
    c2_mean->Print(Form("%s/TID_layers_mean.root",outputDIR.c_str()));
    TCanvas* c2_mpv = prepareCanvas("TID_layers_mpv",observable);
    plotAll(c2_mpv,TIDlayersMPV);
    c2_mpv->Print(Form("%s/TID_layers_mpv.root",outputDIR.c_str()));
    
    /// Layer TOB
    TCanvas* c3_mean = prepareCanvas("TOB_layers_mean",observable);
    plotAll(c3_mean,TOBlayersMean);
    c3_mean->Print(Form("%s/TOB_layers_mean.root",outputDIR.c_str()));
    TCanvas* c3_mpv = prepareCanvas("TOB_layers_mpv",observable);
    plotAll(c3_mpv,TOBlayersMPV);
    c3_mpv->Print(Form("%s/TOB_layers_mpv.root",outputDIR.c_str()));

    /// Layer TECP thin sensors
    TCanvas* c4_mean = prepareCanvas("TECPt_layers_mean",observable);
    plotAll(c4_mean,TECPtlayersMean);
    c4_mean->Print(Form("%s/TECPt_layers_mean.root",outputDIR.c_str()));
    TCanvas* c4_mpv = prepareCanvas("TECPt_layers_mpv",observable);
    plotAll(c4_mpv,TECPtlayersMPV);
    c4_mpv->Print(Form("%s/TECPt_layers_mpv.root",outputDIR.c_str()));
    
    /// Layer TECP thick sensors
    TCanvas* c5_mean = prepareCanvas("TECPT_layers_mean",observable);
    plotAll(c5_mean,TECPTlayersMean);
    c5_mean->Print(Form("%s/TECPT_layers_mean.root",outputDIR.c_str()));
    TCanvas* c5_mpv = prepareCanvas("TECPT_layers_mpv",observable);
    plotAll(c5_mpv,TECPTlayersMPV);
    c5_mpv->Print(Form("%s/TECPT_layers_mpv.root",outputDIR.c_str()));

    /// Layer TECM thin sensors
    TCanvas* c6_mean = prepareCanvas("TECMt_layers_mean",observable);
    plotAll(c6_mean,TECMtlayersMean);
    c6_mean->Print(Form("%s/TECMt_layers_mean.root",outputDIR.c_str()));
    TCanvas* c6_mpv = prepareCanvas("TECMt_layers_mpv",observable);
    plotAll(c6_mpv,TECMtlayersMPV);
    c6_mpv->Print(Form("%s/TECMt_layers_mpv.root",outputDIR.c_str()));

    /// Layer TECM thick sensors
    TCanvas* c7_mean = prepareCanvas("TECMT_layers_mean",observable);
    plotAll(c7_mean,TECMTlayersMean);
    c7_mean->Print(Form("%s/TECMT_layers_mean.root",outputDIR.c_str()));
    TCanvas* c7_mpv = prepareCanvas("TECMT_layers_mpv",observable);
    plotAll(c7_mpv,TECMTlayersMPV);
    c7_mpv->Print(Form("%s/TECMT_layers_mpv.root",outputDIR.c_str()));    

    std::vector<TH1F* > alllayersMean;
    for(auto tib : TIBlayersMean)
      alllayersMean.push_back(tib.second);
    for(auto tob : TOBlayersMean)
      alllayersMean.push_back(tob.second);
    for(auto tid : TIDlayersMean)
      alllayersMean.push_back(tid.second);
    for(auto tec : TECPTlayersMean)
      alllayersMean.push_back(tec.second);
    for(auto tec : TECPtlayersMean)
      alllayersMean.push_back(tec.second);
    for(auto tec : TECMTlayersMean)
      alllayersMean.push_back(tec.second);
    for(auto tec : TECMtlayersMean)
      alllayersMean.push_back(tec.second);  

    std::vector<TH1F* > alllayersMPV;
    for(auto tib : TIBlayersMPV)
      alllayersMPV.push_back(tib.second);
    for(auto tob : TOBlayersMPV)
      alllayersMPV.push_back(tob.second);
    for(auto tid : TIDlayersMPV)
      alllayersMPV.push_back(tid.second);
    for(auto tec : TECPTlayersMPV)
      alllayersMPV.push_back(tec.second);
    for(auto tec : TECPtlayersMPV)
      alllayersMPV.push_back(tec.second);
    for(auto tec : TECMTlayersMPV)
      alllayersMPV.push_back(tec.second);
    for(auto tec : TECMtlayersMPV)
      alllayersMPV.push_back(tec.second);

    // store all the different plots --> plot maximum value for each layer    
    TCanvas* c8_mean (new TCanvas("c_layers_mean","",800,650));  
    plotMaxima(c8_mean,alllayersMean,outputDIR,"layers_mean");
    TCanvas* c8_mpv (new TCanvas("c_layers_mpv","",800,650));  
    plotMaxima(c8_mpv,alllayersMPV,outputDIR,"layers_mpv");

    TFile* allLayerMeanFile = new TFile((outputDIR+"/outputAllLayerMean.root").c_str(),"RECREATE");
    allLayerMeanFile->cd();
    for(auto ihist: alllayersMean)
      ihist->Write();
    allLayerMeanFile->Close();
    
    TFile* allLayerMPVFile = new TFile((outputDIR+"/outputAllLayerMPV.root").c_str(),"RECREATE");
    allLayerMPVFile->cd();
    for(auto ihist: alllayersMPV)
      ihist->Write();
    allLayerMPVFile->Close();

    TIBlayers.clear();
    TOBlayers.clear();
    TIDlayers.clear();
    TECPTlayers.clear();
    TECPtlayers.clear();
    TECMTlayers.clear();
    TECMtlayers.clear();

    TIBlayersMean.clear();
    TOBlayersMean.clear();
    TIDlayersMean.clear();
    TECPTlayersMean.clear();
    TECPtlayersMean.clear();
    TECMTlayersMean.clear();
    TECMtlayersMean.clear();
    alllayersMean.clear();
    
    TIBlayersMPV.clear();
    TOBlayersMPV.clear();
    TIDlayersMPV.clear();
    TECPTlayersMPV.clear();
    TECPtlayersMPV.clear();
    TECMTlayersMPV.clear();
    TECMtlayersMPV.clear();
    alllayersMPV.clear();
  }

  // Per rings
  if(plotSlices){

    // make the hisograms
    RPlots(fileNames,nEvents,delayCorrections,observable,outputDIR);
    
    // Per ring in TIB
    TCanvas* c1b_mean = prepareCanvas("TIB_distance_mean",observable);
    plotAll(c1b_mean,TIBrsMean);
    c1b_mean->Print(Form("%s/TIB_distance_mean.root",outputDIR.c_str()));
    TCanvas* c1b_mpv = prepareCanvas("TIB_distance_mpv",observable);
    plotAll(c1b_mpv,TIBrsMPV);
    c1b_mpv->Print(Form("%s/TIB_distance_mpv.root",outputDIR.c_str()));

    // Per ring in TID
    TCanvas* c2b_mean = prepareCanvas("TID_distance_mean",observable);
    plotAll(c2b_mean,TIDrsMean);
    c2b_mean->Print(Form("%s/TID_distance_mean.root",outputDIR.c_str()));
    TCanvas* c2b_mpv = prepareCanvas("TID_distance_mpv",observable);
    plotAll(c2b_mpv,TIDrsMPV);
    c2b_mpv->Print(Form("%s/TID_distance_mpv.root",outputDIR.c_str()));

    // Per ring in TOB
    TCanvas* c3b_mean = prepareCanvas("TOB_distance_mean",observable);
    plotAll(c3b_mean,TOBrsMean);
    c3b_mean->Print(Form("%s/TOB_distance_mean.root",outputDIR.c_str()));
    TCanvas* c3b_mpv = prepareCanvas("TOB_distance_mpv",observable);
    plotAll(c3b_mpv,TOBrsMPV);
    c3b_mpv->Print(Form("%s/TOB_distance_mpv.root",outputDIR.c_str()));

    // Per ring in TECM
    TCanvas* c4b_mean = prepareCanvas("TECMT_distance_mean",observable);
    plotAll(c4b_mean,TECMTrsMean);
    c4b_mean->Print(Form("%s/TECMT_distance_mean.root",outputDIR.c_str()));
    TCanvas* c4b_mpv = prepareCanvas("TECMT_distance_mpv",observable);
    plotAll(c4b_mpv,TECMTrsMPV);
    c4b_mpv->Print(Form("%s/TECMT_distance_mpv.root",outputDIR.c_str()));
    TCanvas* c5b_mean = prepareCanvas("TECMt_distance_mean",observable);
    plotAll(c5b_mean,TECMtrsMean);
    c5b_mean->Print(Form("%s/TECMt_distance_mean.root",outputDIR.c_str()));
    TCanvas* c5b_mpv = prepareCanvas("TECMt_distance_mpv",observable);
    plotAll(c5b_mpv,TECMtrsMPV);
    c5b_mpv->Print(Form("%s/TECMt_distance_mpv.root",outputDIR.c_str()));

    // Per ring in TECM
    TCanvas* c6b_mean = prepareCanvas("TECPT_distance_mean",observable);
    plotAll(c6b_mean,TECPTrsMean);
    c6b_mean->Print(Form("%s/TECPT_distance_mean.root",outputDIR.c_str()));
    TCanvas* c6b_mpv = prepareCanvas("TECPT_distance_mpv",observable);
    plotAll(c6b_mpv,TECPTrsMPV);
    c6b_mpv->Print(Form("%s/TECPT_distance_mpv.root",outputDIR.c_str()));
    TCanvas* c7b_mean = prepareCanvas("TECPt_distance_mean",observable);
    plotAll(c7b_mean,TECPtrsMean);
    c7b_mean->Print(Form("%s/TECPt_distance_mean.root",outputDIR.c_str()));
    TCanvas* c7b_mpv = prepareCanvas("TECPt_distance_mpv",observable);
    plotAll(c7b_mpv,TECPtrsMPV);
    c7b_mpv->Print(Form("%s/TECPt_distance_mpv.root",outputDIR.c_str()));

    std::vector<TH1F* > allrsMean;
    for(auto tib : TIBrsMean)
      allrsMean.push_back(tib.second);
    for(auto tob : TOBrsMean)
      allrsMean.push_back(tob.second);
    for(auto tid : TIDrsMean)
      allrsMean.push_back(tid.second);
    for(auto tec : TECPTrsMean)
      allrsMean.push_back(tec.second);
    for(auto tec : TECPtrsMean)
      allrsMean.push_back(tec.second);
    for(auto tec : TECMTrsMean)
      allrsMean.push_back(tec.second);
    for(auto tec : TECMtrsMean)
      allrsMean.push_back(tec.second);

    std::vector<TH1F* > allrsMPV;
    for(auto tib : TIBrsMPV)
      allrsMPV.push_back(tib.second);
    for(auto tob : TOBrsMPV)
      allrsMPV.push_back(tob.second);
    for(auto tid : TIDrsMPV)
      allrsMPV.push_back(tid.second);
    for(auto tec : TECPTrsMPV)
      allrsMPV.push_back(tec.second);
    for(auto tec : TECPtrsMPV)
      allrsMPV.push_back(tec.second);
    for(auto tec : TECMTrsMPV)
      allrsMPV.push_back(tec.second);
    for(auto tec : TECMtrsMPV)
      allrsMPV.push_back(tec.second);

    TCanvas* c8b_mean (new TCanvas("c_rings_mean","",800,650));  
    plotMaxima(c8b_mean,allrsMean,outputDIR,"rings");
    TCanvas* c8b_mpv (new TCanvas("c_rings_mpv","",800,650));  
    plotMaxima(c8b_mpv,allrsMPV,outputDIR,"rings");
    
    TIBrs.clear();
    TIDrs.clear();
    TOBrs.clear();
    TECPTrs.clear();
    TECPtrs.clear();
    TECMTrs.clear();
    TECMtrs.clear();

    TIBrsMean.clear();
    TIDrsMean.clear();
    TOBrsMean.clear();
    TECPTrsMean.clear();
    TECPtrsMean.clear();
    TECMTrsMean.clear();
    TECMtrsMean.clear();
    allrsMean.clear();

    TIBrsMPV.clear();
    TIDrsMPV.clear();
    TOBrsMPV.clear();
    TECPTrsMPV.clear();
    TECPtrsMPV.clear();
    TECMTrsMPV.clear();
    TECMtrsMPV.clear();
    allrsMPV.clear();
    
  }

  // Plot per partition 
  if(plotPartitions){

    //change ring definition to collapse all of them
    TIBRing.nDivision = 1;
    TIDRing.nDivision = 1;
    TOBRing.nDivision = 1;
    TECRing.nDivision = 1;
    
    // cumulate all rings
    RPlots(fileNames,nEvents,delayCorrections,observable,outputDIR);

    // create the plots per partition
    std::vector<TH1F* > allPartitionMean;
    for(auto tib : TIBrsMean)
      allPartitionMean.push_back(tib.second);
    for(auto tob : TOBrsMean)
      allPartitionMean.push_back(tob.second);
    for(auto tid : TIDrsMean)
      allPartitionMean.push_back(tid.second);
    for(auto tec : TECPTrsMean)
      allPartitionMean.push_back(tec.second);
    for(auto tec : TECPtrsMean)
      allPartitionMean.push_back(tec.second);
    for(auto tec : TECMTrsMean)
      allPartitionMean.push_back(tec.second);
    for(auto tec : TECMtrsMean)
      allPartitionMean.push_back(tec.second);

    std::vector<TH1F* > allPartitionMPV;
    for(auto tib : TIBrsMPV)
      allPartitionMPV.push_back(tib.second);
    for(auto tob : TOBrsMPV)
      allPartitionMPV.push_back(tob.second);
    for(auto tid : TIDrsMPV)
      allPartitionMPV.push_back(tid.second);
    for(auto tec : TECPTrsMPV)
      allPartitionMPV.push_back(tec.second);
    for(auto tec : TECPtrsMPV)
      allPartitionMPV.push_back(tec.second);
    for(auto tec : TECMTrsMPV)
      allPartitionMPV.push_back(tec.second);
    for(auto tec : TECMtrsMPV)
      allPartitionMPV.push_back(tec.second);


    TCanvas* c1_mean = prepareCanvas("Partitions_mean",observable);
    plotAll(c1_mean,allPartitionMean,"ring1");
    c1_mean->Print(Form("%s/Partitions_mean.root",outputDIR.c_str()));

    TCanvas* c1_mpv = prepareCanvas("Partitions_mpv",observable);
    plotAll(c1_mpv,allPartitionMPV,"ring1");
    c1_mpv->Print(Form("%s/Partitions_mpv.root",outputDIR.c_str()));

    TFile* outputAllPartitionMean = new TFile((outputDIR+"/outputAllPartitionMean.root").c_str(),"RECREATE");
    outputAllPartitionMean->cd();
    for(auto ihist : allPartitionMean)
      ihist->Write();
    outputAllPartitionMean->Close();

    TFile* outputAllPartitionMPV = new TFile((outputDIR+"/outputAllPartitionMPV.root").c_str(),"RECREATE");
    outputAllPartitionMPV->cd();
    for(auto ihist : allPartitionMPV)
      ihist->Write();
    outputAllPartitionMPV->Close();
    
    TIBrs.clear();
    TIDrs.clear();
    TOBrs.clear();
    TECPTrs.clear();
    TECPtrs.clear();
    TECMTrs.clear();
    TECMtrs.clear();

    TIBrsMean.clear();
    TIDrsMean.clear();
    TOBrsMean.clear();
    TECPTrsMean.clear();
    TECPtrsMean.clear();
    TECMTrsMean.clear();
    TECMtrsMean.clear();

    TIBrsMPV.clear();
    TIDrsMPV.clear();
    TOBrsMPV.clear();
    TECPTrsMPV.clear();
    TECPtrsMPV.clear();
    TECMTrsMPV.clear();
    TECMtrsMPV.clear();

    allPartitionMean.clear();
    allPartitionMPV.clear();
  } 

  outputFitFile->Close();
}

