#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <iostream>
#include <TFractionFitter.h>
#include <TCanvas.h>
#include <TPad.h>



void FindChi(){
  TFile *f1 = new TFile("PromptRecoPhotonCorPfMet.root", "READ");
  int nBins  = 15 ;
  float chi2 = 0;
  float fitchi2 = 0;
  int NDF = 0;
  TH1F *h_met_rereco=(TH1F*)f1->Get("h_met_rereco");
  //// rebinning
  
  h_met_rereco->Rebin(4);
  

  //// adding two histogram
  //TH1F *eeffMET = new TH1F("eeffMET", "ee and ff MET added together", 200, 0., 200.);
  

  ///// normalization
  Float_t t1=h_met_rereco->Integral(0,50);
  h_met_rereco->Scale(1./t1);
  h_met_rereco->SetLineColor(kGreen+3);

  TFile *F3 = new TFile("UnnormalizedAsymmErrorCorFake.root", "READ");
  TGraphErrors *g_fake_cor=(TGraphErrors*)F3->Get("unnormalizedrerecoffcorerror");
  
  TH1F *h_fake_cor = new TH1F("h_fake_cor","Fake MET after reweighting",50,0.0,50.0);
  
  float sum = 0;
  for(int i=0;i<50;++i){
    float y = g_fake_cor->Eval(i*4);
    sum = sum + y;
    printf("y for %d:%f\n",i,y);
    float y1 = g_fake_cor->GetErrorY(i);
    printf("error of fake for %d:%f\n",i,y1);
    h_fake_cor->SetBinContent(i*4,y);
    h_fake_cor->SetBinError(i*4,y1);
  }

  TFile *F4 = new TFile("UnnormalizedAsymmErrorCorElectron.root", "READ");
  TGraphErrors *g_electron_cor=(TGraphErrors*)F4->Get("unnormalizedrerecoeecorerror");
  
  
  TH1F *h_electron_cor = new TH1F("h_electron_cor","Electron MET after reweighting",80,0.0,80.0);
  
  //float sum = 0;
  for(int i=0;i<50;++i){
    float y = g_electron_cor->Eval(i*4);
    
    //sum = sum + y;
    printf("y for %d:%f\n",i,y);
    float y1 = g_electron_cor->GetErrorY(i);
    
    h_electron_cor->SetBinContent(i*4,y);
    h_electron_cor->SetBinError(i*4,y1);
    
  }


  /*Float_t t2=h_met_ffSample_RhoReweighted->Integral(0,50);
  h_met_ffSample_RhoReweighted->Scale(1./t2);
  Float_t t3=h_met_candidate->Integral(0,50);
  h_met_candidate->Scale(1./t3);
  
  

  eeffMET->Add(h_met_eeSample_PurityReweighted);
  eeffMET->Add(h_met_ffSample_PurityReweighted);
  eeffMET->Rebin(4);
  Float_t t6=eeffMET->Integral(0,50);
  eeffMET->Scale(1./t6);*/

  /// setting color
  h_met_rereco->SetLineColor(kGreen+3);
  h_electron_cor->SetLineColor(4);
  h_fake_cor->SetLineColor(2);
  //eeffMET->SetLineColor(2);
  
  
  //h_met_ffSample_RhoReweighted->Draw();
  
  //h_met_eeSample_diEMPtRhoReweighted->Draw();
  //h_met_ffSample_RhoReweighted->Draw("sames");
  //eeffMET->Draw("sames");

  for(int i=1; i<=nBins; ++i){
    float x   = h_electron_cor->GetBinContent(i);
    //printf("x:%f for bin %d\n",x,i);
    float dx  = h_electron_cor->GetBinError(i);
    //printf("dx:%f for bin %d\n",dx,i);
    float y   = h_met_rereco->GetBinContent(i);
    //printf("y:%f for bin %d\n",y,i);
    float dy  = h_met_rereco->GetBinError(i);
    //printf("dy:%f for bin %d\n",dy,i);
    if(sqrt(dx*dx+dy*dy)!=0){
      float chi = (x-y)/(sqrt(dx*dx+dy*dy));
      printf("%f,\n",chi);
      chi2     += chi*chi;
    }
  }
  printf("Chi2: %f for Bins upto: %d\n",chi2,nBins);
  printf("Chi2/dof %f for Bins upto: %d\n",chi2/nBins,nBins);


  TObjArray *component = new TObjArray(2);


  component->Add(h_electron_cor);
  component->Add(h_fake_cor);
   
   TFractionFitter* fit = new TFractionFitter(h_met_rereco, component); // initialise
   fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
   fit->SetRangeX(0,55);                    // use only the first 15 bins in the fit
   Int_t status = fit->Fit();              // perform the fit
   fitchi2=fit->GetChisquare();
   NDF=fit->GetNDF();
   std::cout << "fit status: " << status << std::endl;
   std::cout << "chi2:       " << fitchi2 << std::endl;
   std::cout << "chi2/ndf:   " << fitchi2/NDF << std::endl;
   TH1F *CandidateOverFit = new TH1F("CandidateOverFit", "Candidate diEMPt to fit diEMPt ratio", 200, 0.0, 200.0);
   if (status == 0) {                       // check on fit status
     TH1F* result = (TH1F*) fit->GetPlot();
     h_met_rereco->SetLineColor(4);
     h_met_rereco->Draw("Ep");
     result->SetLineColor(2);
     result->Draw("same");
     //TH1F *CandidateOverFit = new TH1F("CandidateOverFit", "Candidate MET to fit MET ratio", 200, 0.0, 200.0);
     //CandidateOverFit->Sumw2();
     CandidateOverFit->Divide(h_met_rereco, result); 
     //CandidateOverFit->Draw();
   }
  
  

}
