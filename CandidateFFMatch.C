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


void CandidateFFMatch(){
  TFile *F2 = new TFile("PromptRecoPhotonCorPfMet.root","READ");
  TH1F *h_met_rereco=(TH1F*)F2->Get("h_met_rereco");
  h_met_rereco->Sumw2();
  h_met_rereco->Rebin(4);
  
  Float_t x1 = h_met_rereco->Integral(0,50);
  h_met_rereco->Scale(1./x1);
  
  TFile *F3 = new TFile("AsymmErrorCorFake.root", "READ");
  TGraphErrors *g_fake_cor=(TGraphErrors*)F3->Get("rerecoffcorerror");
  
  TH1F *h_fake_cor = new TH1F("h_fake_cor","Fake MET after reweighting",50,0.0,50.0);
  
  float sum = 0;
  for(int i=0;i<80;++i){
    float y = g_fake_cor->Eval(i*4);
    sum = sum + y;
    printf("y for %d:%f\n",i,y);
    float y1 = g_fake_cor->GetErrorY(i);
    printf("error of fake for %d:%f\n",i,y1);
    h_fake_cor->SetBinContent(i*4,y);
    h_fake_cor->SetBinError(i*4,y1);
  }
  printf("sum:%f\n",sum);
  float see = h_fake_cor->Integral(0,200);
  printf("see:%f\n",see);
  


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
   //Fakerhoreweighted_rereco->Draw();
   gStyle->SetOptStat(111111);
   
   

   gStyle->SetStatTextColor(2);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_fake_cor->SetLineColor(2);
   h_fake_cor->Draw("sames");
   gStyle->SetOptStat(111111);
   
   
   
   
   Q->cd();
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad2->Draw();
   pad2->cd();
   TH1F *h2 = new TH1F("h2","fake over candidate",50,0,50);
   for(int i =0 ; i<50; ++i){
     float y = g_fake_cor->Eval(i*4);
     float ery= g_fake_cor->GetErrorY(i);
     float y2= h_met_rereco->GetBinContent(i);
     float ery2= h_met_rereco->GetBinError(i);
     if(y2!=0)h2->SetBinContent(i,y/y2);
     float erz=0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     printf("error%d:%f\n",i,erz);
     h2->SetBinError(i,erz);

   }
   
   //h2->Sumw2();
   
   //h2->SetStats(0);
   

   for(int i = 0; i< 80 ; ++i){
     float so = h2->GetBinContent(i);
     printf("so:%f\n",so);
   }
   //gStyle->SetStatTextColor(kGreen+3);
   
   h2->GetXaxis()->SetRangeUser(0,20);
   h2->GetYaxis()->SetRangeUser(0.2,2.0);
   //h2->SetMarkerStyle(21);
   h2->Draw("ep");
   
   
   
   
   //g2->Draw("AP");
   
   Q->SaveAs("CandidateReweightedFakeCorPfMetReReco.eps");
   Q->SaveAs("CandidateReweightedFakeCorPfMetReReco.pdf");

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
