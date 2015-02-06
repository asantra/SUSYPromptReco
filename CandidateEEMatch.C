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


void CandidateEEMatch(){
  TFile *F2 = new TFile("PromptRecoPhotonCorPfMet.root","READ");
  TH1F *h_met_rereco=(TH1F*)F2->Get("h_met_rereco");
  h_met_rereco->Sumw2();
  h_met_rereco->Rebin(4);
  
  Float_t x1 = h_met_rereco->Integral(0,50);
  printf("integral200:%f\n", x1);
  //h_met_rereco->Scale(1./x1);
  
  TFile *F3 = new TFile("Normalized40AsymmErrorCorElectron.root", "READ");
  TGraphErrors *g_electron_cor=(TGraphErrors*)F3->Get("normalized40rerecoeecorerror");// normalized to MET 40 of data
  TFile *F4 = new TFile("Normalized40AsymmErrorCorFake.root", "READ");
  TGraphErrors *g_fake_cor=(TGraphErrors*)F4->Get("normalized40rerecoffcorerror");
  
  TH1F *h_electron_cor = new TH1F("h_electron_cor","Electron MET after reweighting",80,0.0,80.0);
  TH1F *h_fake_cor = new TH1F("h_fake_cor","Fake MET after reweighting",80,0.0,80.0);
  TH1F *h_mixture_cor = new TH1F("h_mixture_cor","Electron+Fake MET after purity reweighting",80,0.0,80.0);
  Double_t Err;
  float sum   = 0, summ40 = 0, suma40 = 0, suma60 = 0, toterr = 0;
  float see   = h_met_rereco->Integral(0,10);
  //float see60 = h_met_rereco->Integral(15,50);
  float see60 = h_met_rereco->IntegralAndError(15,50,Err);
  const float frac1 = 0.466175 ;
  for(int i=0;i<50;++i){
    float y = g_electron_cor->Eval(i*4);
    float z = g_fake_cor->Eval(i*4);
    //printf("y for %d:%f\n",i,y);
    /**/
    float r = frac1*y+(1-frac1)*z;
    if(i<=10)summ40 = summ40 + r;
    if(i>=15)sum = sum + r;
  }
  for(int i=0;i<50;++i){
    float y = g_electron_cor->Eval(i*4);
    float z = g_fake_cor->Eval(i*4);

    float r = frac1*y+(1-frac1)*z;
    r = r*see/summ40;
    if(i<=10)suma40 = suma40 + r;
    if(i>=15)suma60 = suma60 + r;
    float y1 = g_electron_cor->GetErrorY(i);
    float z1 = g_fake_cor->GetErrorY(i);
    float err = sqrt(frac1*frac1*y1*y1+(1-frac1)*(1-frac1)*z1*z1);
    if(i>=15)toterr = toterr + err*err*err*err;
    h_electron_cor->SetBinContent(i*4,y);
    h_electron_cor->SetBinError(i*4,y1);
    h_fake_cor->SetBinContent(i*4,z);
    h_fake_cor->SetBinError(i*4,z1);


    h_mixture_cor->SetBinContent(i*4,r);
    h_mixture_cor->SetBinError(i*4,err);
  }
  printf("unreweighted mixture integral after 60:%f\n",sum);
  printf("unreweighted mixture integral upto 40 :%f\n",summ40);
  printf("reweighted mixture integral upto 40   :%f\n",suma40);
  printf("reweighted mixture integral above 60  :%f\n",suma60);
  printf("reweighted mixture integral error 60  :%f\n",sqrt(toterr));
  printf("candidate integral upto 40            :%f\n",see);
  printf("candidate integral after 60           :%f\n",see60);
  printf("candidate integral error after 60     :%f\n",Err);
  


  TCanvas *Q = new TCanvas("Q","Comparison",900,900); //1200,900
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad1->Draw();
   pad1->cd();
   
   gStyle->SetStatTextColor(kGreen+3);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2); 
   h_met_rereco->SetLineColor(kGreen+3);
   h_met_rereco->GetXaxis()->SetRangeUser(0,80);
   h_met_rereco->Draw();
   //Fakerhoreweighted_prompt->Draw();
   gStyle->SetOptStat(111111);
   
   

   gStyle->SetStatTextColor(2);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_mixture_cor->SetLineColor(2);
   h_mixture_cor->Draw("sames");
   gStyle->SetOptStat(111111);
   
   
   
   
   Q->cd();
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad2->Draw();
   pad2->cd();
   TH1F *h2 = new TH1F("h2","candidate over mixture",200,0,200);
   for(int i =0 ; i<50; ++i){
     /*float y = g_electron_cor->Eval(i*4);
     float ery= g_electron_cor->GetErrorY(i);
     float y2= g_fake_cor->Eval(i*4);
     float ery2= g_fake_cor->GetErrorY(i);
     if(y2!=0)h2->SetBinContent(i,y/y2);
     float erz=0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     printf("error%d:%f\n",i,erz);
     h2->SetBinError(i,erz);*/

     float y = h_met_rereco->GetBinContent(i);
     float ery= h_met_rereco->GetBinError(i);
     float y2= h_mixture_cor->GetBinContent(i*4);
     float ery2= h_mixture_cor->GetBinError(i*4);
     if(y2!=0)h2->SetBinContent(i*4,y/y2);
     float erz=0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     //printf("error%d:%f\n",i,erz);
     h2->SetBinError(i*4,erz);

   }
   
   //h2->Sumw2();
   
   //h2->SetStats(0);
   

   for(int i = 0; i< 50 ; ++i){
     float so = h2->GetBinContent(i);
     //printf("so:%f\n",so);
   }
   //gStyle->SetStatTextColor(kGreen+3);
   
   h2->GetXaxis()->SetRangeUser(0,80);
   h2->GetYaxis()->SetRangeUser(0.2,2.0);
   //h2->SetMarkerStyle(21);
   h2->Draw("ep");
   
   
   
   
   //g2->Draw("AP");
   
   //Q->SaveAs("CandidateReweighted40Electron0_46Fake0_54MixtureCorPfMetReReco.eps");
   //Q->SaveAs("CandidateReweighted40Electron0_46Fake0_54MixtureCorPfMetReReco.pdf");

   TCanvas *R = new TCanvas("R","Ratio",600,450); //1200,900
   R->cd();
   h2->Draw();

  /*for(int i =0 ;i<50; ++i){
    Float_t y1 = g_electron_cor->Eval(i);
    Float_t y2 = g_electron_cor->GetErrorY(i);
    h_electron_cor->SetBinContent(i,y1);
    h_electron_cor->SetBinError(i,y2);
  }*/

  //h_electron_cor->Draw();
  return;
  
 

}
