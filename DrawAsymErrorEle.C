#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <iostream>


void DrawAsymErrorEle(){
  
  TFile *F22 = new TFile("PromptRecoPhotonCorPfMet.root","READ");
  TH1F *h_met_rereco = (TH1F*)F22->Get("h_met_rereco");
  float integral40 = h_met_rereco->Integral(0,40);
  printf("integral40: %f\n",integral40);


  TFile *F2 = new TFile("PromptRecoElectronCorPfMet.root","READ");
  
  TH1F *RhoRatioWeightedee[1000];
  TH1F *RhoRatioWeightedeeprompt[1000];
  TH1F *DiEMPtRatioWeightedee[1000];
  TH1F *DiEMPtRatioWeightedeeprompt[1000];

  TH1F *h_met_rereco_diEMPtrhodiEMPtreweighted = (TH1F*)F2->Get("h_met_rereco_diEMPtrhodiEMPtreweighted");
  TH1F *h_met_prompt_diEMPtrhodiEMPtreweighted = (TH1F*)F2->Get("h_met_prompt_diEMPtrhodiEMPtreweighted"); 
  
  h_met_rereco_diEMPtrhodiEMPtreweighted->Rebin(4);
  h_met_prompt_diEMPtrhodiEMPtreweighted->Rebin(4);
  
  for(int g=0; g<1000; ++g){
    char *nameGraph         = new char[50];
    char *nameGraph2        = new char[50];
    char *nameGraph3        = new char[50];
    char *nameGraph4        = new char[50];
    
    sprintf(nameGraph,"RhoReweightedeeReReco%d",g+1);
    RhoRatioWeightedee[g] = (TH1F*)F2->Get(nameGraph);
    RhoRatioWeightedee[g]->Rebin(4);

    sprintf(nameGraph2,"RhoReweightedeePrompt%d",g+1);
    RhoRatioWeightedeeprompt[g] = (TH1F*)F2->Get(nameGraph2);
    RhoRatioWeightedeeprompt[g]->Rebin(4);

    sprintf(nameGraph3,"DiEMPtReweightedeeReReco%d",g+1);
    DiEMPtRatioWeightedee[g] = (TH1F*)F2->Get(nameGraph3);
    DiEMPtRatioWeightedee[g]->Rebin(4);

    sprintf(nameGraph4,"DiEMPtReweightedeePrompt%d",g+1);
    DiEMPtRatioWeightedeeprompt[g] = (TH1F*)F2->Get(nameGraph4);
    DiEMPtRatioWeightedeeprompt[g]->Rebin(4);
    
    
  }

  


  float xbindiempt[200][1000]                 = {{0},{0}};
  float xbinrho[200][1000]                    = {{0},{0}};
  float xbindiemptprompt[200][1000]           = {{0},{0}};
  float xbinrhoprompt[200][1000]              = {{0},{0}};
  int N                                       = 50; 
  float Bin[200]                              = {0};
  float xerror1[200]                          = {0};
  float xerror2[200]                          = {0};
  float metvalue[200]                         = {0};
  float metvalueerror[200]                    = {0};
  float metvalueprompt[200]                   = {0};
  float metvalueerrorprompt[200]              = {0};
  float metrhovalue[200]                      = {0};
  float metrhovalueprompt[200]                = {0};
  float normalizedmetvalue[200]               = {0};
  float normalizedmetvalueerror[200]          = {0};
  float normalizedmetvalueprompt[200]         = {0};
  float normalizedmetvalueerrorprompt[200]    = {0};
  float normalizedmetrhovalue[200]            = {0};
  float normalizedmetrhovalueprompt[200]      = {0};

  float centraldiempt[200]                    = {0};
  float centraldiemptprompt[200]              = {0};
  float centralrho[200]                       = {0};
  float centralrhoprompt[200]                 = {0};
  float errorupdiempt[200]                    = {0};
  float errorupdiemptprompt[200]              = {0};
  float erroruprho[200]                       = {0};
  float erroruprhoprompt[200]                 = {0};
  float errordowndiempt[200]                  = {0};
  float errordowndiemptprompt[200]            = {0};
  float errordownrho[200]                     = {0};
  float errordownrhoprompt[200]               = {0};
  float normalizederrorupdiempt[200]          = {0};
  float normalizederrorupdiemptprompt[200]    = {0};
  float normalizederroruprho[200]             = {0};
  float normalizederroruprhoprompt[200]       = {0};
  float normalizederrordowndiempt[200]        = {0};
  float normalizederrordowndiemptprompt[200]  = {0};
  float normalizederrordownrho[200]           = {0};
  float normalizederrordownrhoprompt[200]     = {0};
  float integral                              = 0;
  float integralprompt                        = 0;



  float totalerror[200]                       = {0};
  float totalerrorprompt[200]                 = {0};
  float normalizedtotalerror[200]             = {0};
  float normalizedtotalerrorprompt[200]       = {0};

  
  

  for(int j=0 ; j<=50 ; ++j){
    Bin[j] = j*4; 
    float mindiempt(1000000), minrho(1000000), mindiemptprompt(1000000), minrhoprompt(1000000); 
    float maxdiempt(0), maxrho(0), maxdiemptprompt(0), maxrhoprompt(0);
    float sumdiempt             = 0;
    float sumdiemptprompt       = 0;
    float sumrho                = 0;
    float sumrhoprompt          = 0;
    

    metvalue[j] = h_met_rereco_diEMPtrhodiEMPtreweighted->GetBinContent(j);
    metvalueerror[j] = h_met_rereco_diEMPtrhodiEMPtreweighted->GetBinError(j);
    metvalueprompt[j] = h_met_prompt_diEMPtrhodiEMPtreweighted->GetBinContent(j);
    metvalueerrorprompt[j] = h_met_prompt_diEMPtrhodiEMPtreweighted->GetBinError(j);
    if(j<=10)integral = integral + metvalue[j]; // upto 40 GeV
    integralprompt = integralprompt + metvalueprompt[j];
  
    for(int k=0 ; k<1000;++k){
      xbinrho[j][k]       = RhoRatioWeightedee[k]->GetBinContent(j);
      sumrho = sumrho+ xbinrho[j][k];
      if(maxrho < xbinrho[j][k])maxrho=xbinrho[j][k];
      if(minrho > xbinrho[j][k])minrho=xbinrho[j][k];
  
      xbinrhoprompt[j][k] = RhoRatioWeightedeeprompt[k]->GetBinContent(j);
      sumrhoprompt = sumrhoprompt+xbinrhoprompt[j][k];
      if(maxrhoprompt < xbinrhoprompt[j][k])maxrhoprompt=xbinrhoprompt[j][k];
      if(minrhoprompt > xbinrhoprompt[j][k])minrhoprompt=xbinrhoprompt[j][k];


      xbindiempt[j][k]       = DiEMPtRatioWeightedee[k]->GetBinContent(j);
      sumdiempt = sumdiempt+ xbindiempt[j][k];
      if(maxdiempt < xbindiempt[j][k])maxdiempt=xbindiempt[j][k];
      if(mindiempt > xbindiempt[j][k])mindiempt=xbindiempt[j][k];
  
      xbindiemptprompt[j][k] = DiEMPtRatioWeightedeeprompt[k]->GetBinContent(j);
      sumdiemptprompt = sumdiemptprompt+xbindiemptprompt[j][k];
      if(maxdiemptprompt < xbindiemptprompt[j][k])maxdiemptprompt=xbindiemptprompt[j][k];
      if(mindiemptprompt > xbindiemptprompt[j][k])mindiemptprompt=xbindiemptprompt[j][k];
      
    }
    
    centralrho[j]   = sumrho/1000.;
    
    //printf("centralrho[%d]: %f\n",j,metvalue[j]);
    erroruprho[j]   = 0.68*(maxrho-centralrho[j]);
    //printf("erroruprho[%d]: %f \n",j,erroruprho[j]);
    errordownrho[j] = 0.68*(centralrho[j]-minrho); 
    //printf("errordownrho[%d]: %f \n",j,errordownrho[j]);

    centraldiempt[j]   = sumdiempt/1000.;
    errorupdiempt[j]   = 0.68*(maxdiempt-centraldiempt[j]);
    errordowndiempt[j] = 0.68*(centraldiempt[j]-mindiempt); 
    

    centralrhoprompt[j]   = sumrhoprompt/1000.;
    erroruprhoprompt[j]   = (maxrhoprompt-centralrhoprompt[j])*0.68;
    errordownrhoprompt[j] = (centralrhoprompt[j]-minrhoprompt)*0.68; 

    centraldiemptprompt[j]   =  sumdiemptprompt/1000.;
    errorupdiemptprompt[j]   = (maxdiemptprompt-centraldiemptprompt[j])*0.68;
    errordowndiemptprompt[j] = (centraldiemptprompt[j]-mindiemptprompt)*0.68;

    
    //printf("metvalue[%d]:%f\n",j,metvalue[j]);
    //printf("metvalueoprompt[%d]:%f\n",j,metvalueprompt[j]);
  
    totalerror[j] = sqrt((errorupdiempt[j]+errordowndiempt[j])*(errorupdiempt[j]+errordowndiempt[j])+(erroruprho[j]+errordownrho[j])*(erroruprho[j]+errordownrho[j]));

    totalerrorprompt[j] = sqrt((errorupdiemptprompt[j]+errordowndiemptprompt[j])*(errorupdiemptprompt[j]+errordowndiemptprompt[j])+(erroruprhoprompt[j]+errordownrhoprompt[j])*(erroruprhoprompt[j]+errordownrhoprompt[j]));
    
  }
  float check = 0;
  float check2 = 0;
  for(int j=0 ; j<=50 ; ++j){
    normalizedmetvalue[j] = metvalue[j]*integral40/integral;//normalized to MET 40 of data
    check = check+normalizedmetvalue[j];
    normalizedmetvalueerror[j] = metvalueerror[j]*integral40/integral;
    normalizedmetvalueprompt[j] = metvalueprompt[j]/integralprompt;
    check2 = check2 + normalizedmetvalueprompt[j];
    normalizedmetvalueerrorprompt[j] = metvalueerrorprompt[j]/integralprompt; 
    //printf("normalizedmetvalue[%d]: %f\n",j, normalizedmetvalue[j]);
    normalizederroruprho[j]   = erroruprho[j]*integral40/integral;
    //printf("normalizederrorup[%d]: %f \n",j,normalizederroruprho[j]);
    normalizederrordownrho[j] = errordownrho[j]*integral40/integral; 
    //printf("normalizederrordown[%d]: %f \n",j,normalizederrordownrho[j]);

    normalizederrorupdiempt[j]   = errorupdiempt[j]*integral40/integral;
    normalizederrordowndiempt[j] = errordowndiempt[j]*integral40/integral; 
    
    normalizederroruprhoprompt[j]   = erroruprhoprompt[j]/integralprompt;
    normalizederrordownrhoprompt[j] = errordownrhoprompt[j]/integralprompt; 

    normalizederrorupdiemptprompt[j]   = errorupdiemptprompt[j]/integralprompt;
    normalizederrordowndiemptprompt[j] = errordowndiemptprompt[j]/integralprompt;

    normalizedtotalerror[j] = sqrt((normalizederrorupdiempt[j]+normalizederrordowndiempt[j])*(normalizederrorupdiempt[j]+normalizederrordowndiempt[j])+(normalizederroruprho[j]+normalizederrordownrho[j])*(normalizederroruprho[j]+normalizederrordownrho[j])+(normalizedmetvalueerror[j]*normalizedmetvalueerror[j]));

    normalizedtotalerrorprompt[j] = sqrt((normalizederrorupdiemptprompt[j]+normalizederrordowndiemptprompt[j])*(normalizederrorupdiemptprompt[j]+normalizederrordowndiemptprompt[j])+(normalizederroruprhoprompt[j]+normalizederrordownrhoprompt[j])*(normalizederroruprhoprompt[j]+normalizederrordownrhoprompt[j])+(normalizedmetvalueerrorprompt[j]*normalizedmetvalueerrorprompt[j]));

  }
  printf("normalized check: %f \n",check);
  printf("normalized check2: %f \n",check2);
  
  TGraphErrors *g1, *g2;

  g1 = new TGraphErrors(N, Bin, normalizedmetvalue, xerror1, normalizedtotalerror);
  g2 = new TGraphErrors(N, Bin, metvalueprompt, xerror1, totalerrorprompt);

  g1->SetMarkerColor(kRed);
  //g2->SetMarkerColor(kOrange);
  //g1->GetXaxis()->SetLimits(0, 80);

  g1->Draw("A*");

  TFile* fout2 = new TFile("Normalized40AsymmErrorCorElectron.root","RECREATE");
  fout2->cd();
  g1->Write("normalized40rerecoeecorerror");
  g2->Write("unnormalizedprompteerawerror");

  //g2->Draw("sames,P");

  
  
  

}
