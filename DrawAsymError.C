#include <TH1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <iostream>


void DrawAsymError(){

  TFile *F22 = new TFile("PromptRecoPhotonCorPfMet.root","READ");
  TH1F *h_met_rereco = (TH1F*)F22->Get("h_met_rereco");
  float integral40 = h_met_rereco->Integral(0,40);
  printf("integral40: %f\n",integral40);

  TFile *F2 = new TFile("PromptRecoFakeCorPfMet.root","READ");
  
  TH1F *RhoRatioWeightedff[1000];
  TH1F *RhoRatioWeightedffPrompt[1000];
  TH1F *h_met_rereco_diEMPtrhoreweighted = (TH1F*)F2->Get("h_met_rereco_diEMPtrhoreweighted");
  TH1F *h_met_prompt_diEMPtrhoreweighted = (TH1F*)F2->Get("h_met_prompt_diEMPtrhoreweighted"); 

  h_met_rereco_diEMPtrhoreweighted->Rebin(4);
  h_met_prompt_diEMPtrhoreweighted->Rebin(4);
  /*Float_t x1 = h_met_rereco_diEMPtrhoreweighted->Integral(0,50);
  h_met_rereco_diEMPtrhoreweighted->Scale(1./x1);

  Float_t x2 = h_met_prompt_diEMPtrhoreweighted->Integral(0,50);
  h_met_prompt_diEMPtrhoreweighted->Scale(1./x2);*/

  for(int g=0; g<1000; ++g){
    char *nameGraph        = new char[50];
    char *nameGraph2        = new char[50];
    sprintf(nameGraph,"RhoReweightedffReReco%d",g+1);
    RhoRatioWeightedff[g] = (TH1F*)F2->Get(nameGraph);
    RhoRatioWeightedff[g]->Rebin(4);
    /*Float_t y1 = RhoRatioWeightedff[g]->Integral(0,50);
    RhoRatioWeightedff[g]->Scale(1./y1);*/

    sprintf(nameGraph2,"RhoReweightedffPrompt%d",g+1);
    RhoRatioWeightedffPrompt[g] = (TH1F*)F2->Get(nameGraph2);
    RhoRatioWeightedffPrompt[g]->Rebin(4);
    
  }

  


  float xbin[200][1000]                 = {{0},{0}};
  float xbinprompt[200][1000]           = {{0},{0}};
  Int_t N                              = 50; 
  float Bin[200]                        = {0};
  float xerror1[200]                    = {0};
  float xerror2[200]                    = {0};
  float metvalue[200]                   = {0};
  float metvalueprompt[200]             = {0};
  float normalizedmetvalue[200]         = {0};
  float normalizedmetvalueprompt[200]   = {0};
  float metvalueerror[200]              = {0};
  float metvalueerrorprompt[200]        = {0};
  float normalizedmetvalueerror[200]    = {0};
  float normalizedmetvalueerrorprompt[200]= {0};
  float totalerror[200]                 = {0};
  float totalerrorprompt[200]           = {0};
  float normalizedtotalerror[200]       = {0};
  float normalizedtotalerrorprompt[200] = {0}; 

  float central[200]                    = {0};
  float centralprompt[200]              = {0};
  float errorup[200]                    = {0};
  float errorupprompt[200]              = {0};
  float errordown[200]                  = {0};
  float errordownprompt[200]            = {0};
  float normalizederrorup[200]          = {0};
  float normalizederrorupprompt[200]    = {0};
  float normalizederrordown[200]        = {0};
  float normalizederrordownprompt[200]  = {0};
  float integral                        = 0;
  float integralprompt                  = 0;
  
  

  for(int j=0 ; j<=50 ; ++j){
    metvalue[j] = h_met_rereco_diEMPtrhoreweighted->GetBinContent(j);
    metvalueerror[j] = h_met_rereco_diEMPtrhoreweighted->GetBinError(j);
    if(j<=10)integral = integral + metvalue[j] ;// getting value upto 40
    metvalueprompt[j] = h_met_prompt_diEMPtrhoreweighted->GetBinContent(j);
    metvalueerrorprompt[j] = h_met_prompt_diEMPtrhoreweighted->GetBinError(j);
    integralprompt = integralprompt + metvalueprompt[j] ;

    Bin[j] = j*4; 
    float min(10000), max(0), minprompt(10000), maxprompt(0);
    float sum = 0;
    float sumprompt = 0;
    for(int k=0 ; k<1000;++k){
      xbin[j][k]       = RhoRatioWeightedff[k]->GetBinContent(j);
      sum = sum + xbin[j][k];
      if(max < xbin[j][k])max=xbin[j][k];
      if(min > xbin[j][k])min=xbin[j][k];
  
      xbinprompt[j][k] = RhoRatioWeightedffPrompt[k]->GetBinContent(j);
      sumprompt = sumprompt+xbinprompt[j][k];
      if(maxprompt < xbinprompt[j][k])maxprompt=xbinprompt[j][k];
      if(minprompt > xbinprompt[j][k])minprompt=xbinprompt[j][k];
      
    }
    
    central[j]   = sum/1000;
    //printf("central[%d]: %f\n",j,central[j]);
    errorup[j]   = 0.68*(max-central[j]);
    //printf("errorup[%d]: %f \n",j,errorup[j]);
    errordown[j] = 0.68*(central[j]-min); 
    //printf("errordown[%d]: %f \n",j,errordown[j]);


    centralprompt[j]   =  sumprompt/1000;
    errorupprompt[j]   = (maxprompt-centralprompt[j])*0.68;
    errordownprompt[j] = (centralprompt[j]-minprompt)*0.68; 

    
    
  }

  printf("intergal: %f\n",integral);
  for(int j=0 ; j<=50 ; ++j){
    normalizedmetvalue[j] = metvalue[j]*integral40/integral; // normalized to MET 40 of the data
    normalizedmetvalueerror[j] = metvalueerror[j]*integral40/integral;
    //printf("normalizedmetvalue[%d]:%f \n", j,normalizedmetvalue[j]);
    normalizederrorup[j]  = errorup[j]*integral40/integral;
    //printf("normalizederrorup[%d]:%f \n", j,normalizederrorup[j]);
    normalizederrordown[j]= errordown[j]*integral40/integral;
    //printf("normalizederrordown[%d]:%f \n", j,normalizederrordown[j]);

    normalizedmetvalueprompt[j] = metvalueprompt[j]/integralprompt;
    normalizedmetvalueerrorprompt[j] = metvalueerrorprompt[j]/integralprompt;
    normalizederrorupprompt[j]  = errorupprompt[j]/integralprompt;
    normalizederrordownprompt[j]= errordownprompt[j]/integralprompt;

    totalerror[j] = sqrt((metvalueerror[j]*metvalueerror[j])+(errorup[j]+errordown[j])*(errorup[j]+errordown[j]));

    totalerrorprompt[j] = sqrt((metvalueerrorprompt[j]*metvalueerrorprompt[j])+(errorupprompt[j]+errordownprompt[j])*(errorupprompt[j]+errordownprompt[j]));


    normalizedtotalerror[j] = sqrt((normalizedmetvalueerror[j]*normalizedmetvalueerror[j])+(normalizederrorup[j]+normalizederrordown[j])*(normalizederrorup[j]+normalizederrordown[j]));

    normalizedtotalerrorprompt[j] = sqrt((normalizedmetvalueerrorprompt[j]*normalizedmetvalueerrorprompt[j])+(normalizederrorupprompt[j]+normalizederrordownprompt[j])*(normalizederrorupprompt[j]+normalizederrordownprompt[j]));
    
  }
  
  
  TGraphErrors *g1, *g2;

  g1 = new TGraphErrors(N, Bin, normalizedmetvalue, xerror1, normalizedtotalerror);
  g2 = new TGraphErrors(N, Bin, metvalueprompt, xerror1, totalerrorprompt);

  //float q = g1->Integral(0,50);
  //g1->Scale(1./q);
  //printf("q: %f",q);
  //g1->SetMarkerStyle(kFullDotSmall);
  //g1->SetMarkerColor(kRed);
  //g2->SetMarkerStyle(kFullDotSmall);
  //g2->SetMarkerColor(kBlue);
  //g1->GetXaxis()->SetLimits(0, 80);
  //g2->GetXaxis()->SetLimits(0, 80);

  g1->Draw("A*");
  //g2->Draw();
  TFile* fout = new TFile("Normalized40AsymmErrorCorFake.root","RECREATE");
  fout->cd();
  g1->Write("normalized40rerecoffcorerror");
  g2->Write("unnormalizedpromptffrawerror");

  
  

  //g2->Draw("sames,P");

  
  /*TCanvas *Q = new TCanvas("Q","Comparison",600,450);
  Q->cd();*/

  /*gStyle->SetStatTextColor(4);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  h_met_ffSample->Draw();
  gStyle->SetOptStat(111111);*/


  /*gStyle->SetStatTextColor(2);
  gStyle->SetStatY(0.5);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  RhoRatioWeightedee[0]->Draw();
  
  for(int j=1;j<1000;++j){
    if(j%20==0)printf("processed:%d\n",j);

    
    RhoRatioWeightedee[j]->Draw("sames");
  }


  Q->SaveAs("RhoErrorPropagatedeeMETReRecoCor.eps");
  Q->SaveAs("RhoErrorPropagatedeeMETReRecoCor.pdf");


  TCanvas *R = new TCanvas("R","Comparison",600,450);
  R->cd();

  /*gStyle->SetStatTextColor(4);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  h_met_ffSample->Draw();
  gStyle->SetOptStat(111111);*/


 /* gStyle->SetStatTextColor(2);
  gStyle->SetStatY(0.5);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  RhoRatioWeightedff[0]->Draw();
  
  for(int j=1;j<1000;++j){
    if(j%20==0)printf("processed:%d\n",j);

    
    RhoRatioWeightedff[j]->Draw("sames");
  }


  R->SaveAs("ErrorPropagatedffMETReRecoCor.eps");
  R->SaveAs("ErrorPropagatedffMETReRecoCor.pdf");*/
  

}
