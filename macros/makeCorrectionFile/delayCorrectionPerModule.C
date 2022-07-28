#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"

#include "../delayUtils.h"

// take as input the output root file produced by delayValidationPerModule (TTree with floating point correction for each detId).
// It creates a new file with a TTree with only Detid, fedChannel, delay in step of 1.04 ns
static int reductionFactor = 100;

void  delayCorrectionPerModule(string fileName, string outputDIR, string outputName, bool saveFits = false, float delayCutForPlotting = 4, string observable = "maxChargeCorrected"){

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  // compute corrections
  std::cout<<"###############################"<<std::endl;
  std::cout<<"#### computing corrections ####"<<std::endl;
  std::cout<<"###############################"<<std::endl;

  // open and create TTreeReader for the input tree
  std::shared_ptr<TFile> inputFile (TFile::Open(fileName.c_str()));
  std::shared_ptr<TTree> inputTree ((TTree*) inputFile->FindObjectAny("delayCorrection"));

  system(("mkdir -p "+outputDIR).c_str());
  
  TTreeReader reader(inputTree.get());
  TTreeReaderValue<unsigned int> Detid_i    (reader,"Detid");
  TTreeReaderValue<unsigned int> dcuId_i    (reader,"dcuId");
  TTreeReaderValue<unsigned int> fedCh_i    (reader,"fedCh");
  TTreeReaderValue<unsigned int> fedId_i    (reader,"fedId");
  TTreeReaderValue<unsigned int> fecRing_i  (reader,"fecRing");
  TTreeReaderValue<unsigned int> ccuAdd_i   (reader,"ccuAdd");
  TTreeReaderValue<unsigned int> ccuChan_i  (reader,"ccuChan");
  TTreeReaderValue<int>   notFound_i (reader,"notFound");
  TTreeReaderValue<int>   nFilledBinRejected_i (reader,"nFilledBinRejected");
  TTreeReaderValue<int>   fitRejected_i (reader,"fitRejected");
  TTreeReaderValue<int>   peakOutBoundaryRejected_i (reader,"peakOutBoundaryRejected");
  TTreeReaderValue<int>   amplitudeRejected_i (reader,"amplitudeRejected");
  TTreeReaderValue<int>   significanceRejected_i (reader,"significanceRejected");
  TTreeReaderValue<int>   sigmaRejected_i (reader,"sigmaRejected");
  TTreeReaderValue<float> delayCorr_i (reader,"delayCorr");
  // to reconstruct the Gaussian fit vs delay for each module
  TTreeReaderValue<float> measuredMeanAmplitude_i (reader,"measuredMeanAmplitude");
  TTreeReaderValue<float> measuredSigma_i (reader,"measuredSigma");
  TTreeReaderValue<float> measuredDelay_i (reader,"measuredDelay");
  TTreeReaderValue<vector<float> >  amplitude (reader,"amplitude");
  TTreeReaderValue<vector<float> >  amplitudeUnc (reader,"amplitudeUnc");

  // output file and output tree structure
  std::shared_ptr<TFile> outputFile (new TFile((outputDIR+"/"+outputName+".root").c_str(),"RECREATE"));
  outputFile->cd();
  std::shared_ptr<TTree> outputTree (new TTree("delayCorrection","delayCorrection"));
  unsigned int Detid;
  unsigned int dcuId;
  unsigned int fedId;
  unsigned int fedCh;
  float    delayCorr;
  outputTree->Branch("Detid",&Detid,"Detid/i");
  outputTree->Branch("dcuId",&dcuId,"dcuId/i");
  outputTree->Branch("fedId",&fedId,"fedId/i");
  outputTree->Branch("fedCh",&fedCh,"fedCh/i");
  outputTree->Branch("delayCorr",&delayCorr,"delayCorr/F");

  // Bad devices
  std::map<unsigned int,float> notFoundMap;
  std::map<unsigned int,float> filledBinsRejectedMap;
  std::map<unsigned int,float> fitRejectedMap;
  std::map<unsigned int,float> peakOutBoundaryRejectedMap;
  std::map<unsigned int,float> amplitudeRejectedMap;
  std::map<unsigned int,float> significanceRejectedMap;
  std::map<unsigned int,float> sigmaRejectedMap;
  std::map<unsigned int,float> rejectedMap;
  std::map<unsigned int,float> rawDelayMap;
  std::map<unsigned int,float> delayMap;
  std::map<unsigned int,std::string> largeDelayCoordinatesMap;
  std::map<unsigned int,std::string> rejectedCoordinatesMap;
  std::vector<unsigned int> detIdMap;

  std::cout<<"### Start loop "<<std::endl;
  // start loop on the input tree
  vector<double> limits;
  setLimitsAndBinning("delay",limits);
  TF1  fitfunc ("fitfunc","[0]*TMath::Gaus(x,[1],[2])",limits.front(),limits.back());
  TH1F profile ("profile","",limits.size()-1,&limits[0]);
  TCanvas canvas ("canvas","",600,600);
  canvas.cd();
  unsigned int ientry = 0;
  while(reader.Next()){

    // move from APV-chip based to dcuId/detid based
    if(std::find(detIdMap.begin(),detIdMap.end(),*Detid_i) != detIdMap.end()) continue;
    detIdMap.push_back(*Detid_i);

    ientry++;
    Detid = *Detid_i;
    fedCh = *fedCh_i;
    delayCorr = std::round(*delayCorr_i*24./25.)*25./24.;// delay in unitis of 24/25        
    if(delayCorr == -0) delayCorr = 0;
    outputTree->Fill();

    notFoundMap[Detid] = 0;
    rejectedMap[Detid] = 0;
    filledBinsRejectedMap[Detid] = 0;
    fitRejectedMap[Detid] = 0;
    peakOutBoundaryRejectedMap[Detid] = 0;
    amplitudeRejectedMap[Detid] = 0;
    significanceRejectedMap[Detid] = 0;
    significanceRejectedMap[Detid] = 0;
    sigmaRejectedMap[Detid] = 0;
    
    if(*notFound_i){
      notFoundMap[Detid] = 1;
      rejectedCoordinatesMap[Detid] = to_string(*fedId_i)+" "+to_string(*fedCh_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
      continue;
    }
    else if (*nFilledBinRejected_i){
      rejectedMap[Detid] = 1;
      rejectedCoordinatesMap[Detid] = to_string(*fedId_i)+" "+to_string(*fedCh_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
      filledBinsRejectedMap[Detid] = 1;
      continue;
    }
    else if (*fitRejected_i){
      rejectedMap[Detid] = 1;
      rejectedCoordinatesMap[Detid] = to_string(*fedId_i)+" "+to_string(*fedCh_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
      fitRejectedMap[Detid] = 1;
      continue;
    }
    else if(*peakOutBoundaryRejected_i){
      rejectedMap[Detid] = 1;
      rejectedCoordinatesMap[Detid] = to_string(*fedId_i)+" "+to_string(*fedCh_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
      peakOutBoundaryRejectedMap[Detid] = 1;
      continue;
    }
    else if(*amplitudeRejected_i){
      rejectedMap[Detid] = 1;
      rejectedCoordinatesMap[Detid] = to_string(*fedId_i)+" "+to_string(*fedCh_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
      amplitudeRejectedMap[Detid] = 1;
      continue;
    }
    else if(*significanceRejected_i){
      rejectedMap[Detid] = 1;
      rejectedCoordinatesMap[Detid] = to_string(*fedId_i)+" "+to_string(*fedCh_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
      significanceRejectedMap[Detid] = 1;
      continue;
    }
    else if(*sigmaRejected_i){
      rejectedMap[Detid] = 1;
      rejectedCoordinatesMap[Detid] = to_string(*fedId_i)+" "+to_string(*fedCh_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
      sigmaRejectedMap[Detid] = 1;
      continue;
    }
    
    rawDelayMap[Detid] = *delayCorr_i;
    delayMap[Detid]    = delayCorr;

    // save the coordinates
    if(fabs(*delayCorr_i) > delayCutForPlotting)
      largeDelayCoordinatesMap[Detid] = to_string(*fedId_i)+" "+to_string(*fedCh_i)+" "+to_string(*fecRing_i)+" "+to_string(*ccuAdd_i)+" "+to_string(*ccuChan_i);
    
    // save good fits
    if(saveFits and ientry % reductionFactor == 0){
      fitfunc.SetParameter(0,*measuredMeanAmplitude_i);
      fitfunc.SetParameter(1,*measuredDelay_i);
      fitfunc.SetParameter(2,*measuredSigma_i);
      profile.Reset();
      for(int iBin = 1; iBin < profile.GetNbinsX(); iBin++){
	if(amplitude->at(iBin) == 0) continue;
	profile.SetBinContent(iBin,amplitude->at(iBin));
	profile.SetBinError(iBin,amplitudeUnc->at(iBin));
      }
      profile.GetXaxis()->SetTitle("delay [ns]");
      if(TString(observable).Contains("maxCharge"))
	profile.GetXaxis()->SetTitle("Leading strip charge (ADC)");
      else if(observable == "clCorrectedSignalOverNoise")
	profile.GetXaxis()->SetTitle("Cluster S/N");
      else if(observable == "clCorrectedCharge")
	profile.GetXaxis()->SetTitle("Cluster charge (ADC)");
      profile.SetMarkerColor(kBlack);
      profile.SetMarkerSize(1);
      profile.SetMarkerStyle(20);
      profile.Draw("PE");
      fitfunc.SetLineColor(kRed);
      fitfunc.SetLineWidth(2);	  
      fitfunc.Draw("Lsame");
      profile.GetYaxis()->SetRangeUser(fitfunc.GetMinimum(limits.front(),limits.back())*0.75,fitfunc.GetMaximum(limits.front(),limits.back())*1.25);
      CMS_lumi(&canvas,"");
      
      canvas.SaveAs((outputDIR+"/plot_detid_"+to_string(Detid)+".png").c_str(),"png");
      canvas.SaveAs((outputDIR+"/plot_detid_"+to_string(Detid)+".pdf").c_str(),"pdf");
    }
  }

  std::cout<<"### Loop finished "<<std::endl;  

  std::cout<<"### Save tree "<<std::endl;  

  // write output
  outputFile->cd();
  outputTree->BuildIndex("Detid");  
  outputTree->Write();

  std::cout<<"### Make output text file"<<std::endl;

  ofstream rawDelayFile ((outputDIR+"/rawDelayCorrection.txt").c_str());
  for(auto imap : rawDelayMap){
    rawDelayFile << imap.first << "   "<<imap.second<<"\n";
  }
  rawDelayFile.close();

  ofstream delayFile ((outputDIR+"/delayCorrection.txt").c_str());
  for(auto imap : delayMap){
    delayFile << imap.first << "   "<<imap.second<<"\n";
  }
  delayFile.close();

  ofstream largeDelayChannel((outputDIR+"/largeDelayChannel.txt").c_str());
  for(auto imap : largeDelayCoordinatesMap)
    largeDelayChannel << imap.first << "   "<<imap.second<<"\n";
  largeDelayChannel.close();

  ofstream rejectedChannel((outputDIR+"/rejectedChannel.txt").c_str());
  for(auto imap : rejectedCoordinatesMap)
    rejectedChannel << imap.first << "   "<<imap.second<<"\n";
  rejectedChannel.close();

  ofstream notFoundFile((outputDIR+"/notFound.txt").c_str());
  for(auto imap : notFoundMap)
    notFoundFile << imap.first << "   "<<imap.second<<"\n";
  notFoundFile.close();

  ofstream filledBinsRejectedFile((outputDIR+"/filledBinsRejected.txt").c_str());
  for(auto imap : filledBinsRejectedMap)
    filledBinsRejectedFile << imap.first << "   "<<imap.second<<"\n";
  filledBinsRejectedFile.close();

  ofstream fitRejectedFile((outputDIR+"/fitRejected.txt").c_str());
  for(auto imap : fitRejectedMap)
    fitRejectedFile << imap.first << "   "<<imap.second<<"\n";
  fitRejectedFile.close();

  ofstream peakOutBoundaryRejectedFile((outputDIR+"/peakOutBoundaryRejected.txt").c_str());
  for(auto imap : peakOutBoundaryRejectedMap)
    peakOutBoundaryRejectedFile << imap.first << "   "<<imap.second<<"\n";
  peakOutBoundaryRejectedFile.close();

  ofstream amplitudeRejectedFile((outputDIR+"/amplitudeRejected.txt").c_str());
  for(auto imap : amplitudeRejectedMap)
    amplitudeRejectedFile << imap.first << "   "<<imap.second<<"\n";
  amplitudeRejectedFile.close();

  ofstream significanceRejectedFile((outputDIR+"/significanceRejected.txt").c_str());
  for(auto imap : significanceRejectedMap)
    significanceRejectedFile << imap.first << "   "<<imap.second<<"\n";
  significanceRejectedFile.close();

  ofstream sigmaRejectedFile((outputDIR+"/sigmaRejected.txt").c_str());
  for(auto imap : sigmaRejectedMap)
    sigmaRejectedFile << imap.first << "   "<<imap.second<<"\n";
  sigmaRejectedFile.close();

  ofstream rejectedFile((outputDIR+"/rejected.txt").c_str());
  for(auto imap : rejectedMap)
    rejectedFile << imap.first << "   "<<imap.second<<"\n";
  rejectedFile.close();

  // produce tracker maps from text files
  std::cout<<"### Make tracker maps "<<std::endl;
  system(("trackerMapPlot "+outputDIR+"/rawDelayCorrection.txt NULL Delay-Raw 0 -8 8 "+outputDIR+"/rawDelayCorrection").c_str());
  system(("trackerMapPlot "+outputDIR+"/delayCorrection.txt NULL Delay 0 -8 8 "+outputDIR+"/delayCorrection").c_str());
  system(("trackerMapPlot "+outputDIR+"/notFound.txt NULL NotFound 0 0 1 "+outputDIR+"/notFound").c_str());
  system(("trackerMapPlot "+outputDIR+"/filledBinsRejected.txt NULL nBinsRejected 0 0 1 "+outputDIR+"/filledBinsRejected").c_str());
  system(("trackerMapPlot "+outputDIR+"/fitRejected.txt NULL fitRejected 0 0 1 "+outputDIR+"/fitRejected").c_str());
  system(("trackerMapPlot "+outputDIR+"/peakOutBoundaryRejected.txt peakOutRejected Delay 0 0 1 "+outputDIR+"/peakOutBoundaryRejected").c_str());
  system(("trackerMapPlot "+outputDIR+"/amplitudeRejected.txt NULL ampRejected 0 0 1 "+outputDIR+"/amplitudeRejected").c_str());
  system(("trackerMapPlot "+outputDIR+"/significanceRejected.txt NULL signifRejected 0 0 1 "+outputDIR+"/significanceRejected").c_str());
  system(("trackerMapPlot "+outputDIR+"/sigmaRejected.txt NULL sigmaRejected 0 0 1 "+outputDIR+"/sigmaRejected").c_str());
  system(("trackerMapPlot "+outputDIR+"/rejected.txt NULL Rejected 0 0 1 "+outputDIR+"/rejected").c_str());



  return;
}


