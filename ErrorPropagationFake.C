// ----------------------------------------------------------------------------------- //
/*
  ROOT macro for illustrating error propagation using random Gaussian numbers.

  The example is based on first doing the error propagation analytically, and then
  verifying it by running a so-called Monte-Carlo (MC) program.

  For more information on error propagation, see:
    R. J. Barlow: page 48-61 
    P. R. Bevington: page 36-48

  Author: Troels C. Petersen (NBI)
  Email:  petersen@nbi.dk
  Date:   8th of September 2011
*/
// ----------------------------------------------------------------------------------- //

#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TGraphErrors.h>


//---------------------------------------------------------------------------------- 
void ErrorPropagationFake()
//---------------------------------------------------------------------------------- 
{
  


  // Set parameters:
  TRandom3 r;
  TRandom3 r1;
  
  
  
  

  // Make histograms:

  TFile *Fnew = new TFile("FakeRatio.root","READ");
  TH1F *FakeDiEMPtReweightedRho= (TH1F*)Fnew->Get("FakeDiEMPtReweightedRho");
  TH1F *FakeDiEMPtReweightedRhoPrompt= (TH1F*)Fnew->Get("FakeDiEMPtReweightedRhoPrompt");
  
  TFile* fout = new TFile("RandomRatioFF.root","RECREATE");
  
  
  
  TH1F *RatioFakediEMPtReweightedRho[1000]; // for ff rho
  TH1F *RatioFakediEMPtReweightedRhoPrompt[1000];

  
  char *histnameFakediEMPtReweightedRho       = new char[50];
  char *histnameFakediEMPtReweightedRhoPrompt = new char[50];
  char *histtitleFakediEMPtReweightedRho      = new char[50];
  char *histtitleFakediEMPtReweightedRhoPrompt= new char[50];



//----------------------------------------------------------------------------------
// Loop over process:
//----------------------------------------------------------------------------------



  

  double mu1[40]          =  {0};
  double sig1[40]         =  {0};

  double mu2[40]          =  {0};
  double sig2[40]         =  {0};

  

  float rho12 =  0.0; 

  for(int q=0 ; q < 40;++q){
    mu1[q]  = FakeDiEMPtReweightedRho->GetBinContent(q);
    sig1[q] = FakeDiEMPtReweightedRho->GetBinError(q);
    mu2[q]  = FakeDiEMPtReweightedRhoPrompt->GetBinContent(q); 
    sig2[q] = FakeDiEMPtReweightedRhoPrompt->GetBinError(q);  
  }


  if ((rho12 < -1.0) || (rho12 > 1.0)) {
    printf("  ERROR: Correlation factor not in interval [-1,1], as it is %6.2f \n", rho12);
    return;
  }



//----------------------------------------------------------------------------------
// Loop over process:
//----------------------------------------------------------------------------------



  for(int k=0; k < 1000; ++k){
    if(k%50 == 0)printf("processed rho:%d \n", k);//
    sprintf(histnameFakediEMPtReweightedRho, "FakeRhoRandom%d",k+1);
    sprintf(histtitleFakediEMPtReweightedRho,"random rho ratio for ff rereco %d",k+1);
    RatioFakediEMPtReweightedRho[k]=new TH1F(histnameFakediEMPtReweightedRho, histtitleFakediEMPtReweightedRho,40,0.0,40.0);

    sprintf(histnameFakediEMPtReweightedRhoPrompt, "FakeRhoRandomPrompt%d",k+1);
    sprintf(histtitleFakediEMPtReweightedRhoPrompt,"random rho ratio for ff prompt %d",k+1);
    RatioFakediEMPtReweightedRhoPrompt[k]=new TH1F(histnameFakediEMPtReweightedRhoPrompt,histtitleFakediEMPtReweightedRhoPrompt,40,0.0,40.0);
    



    for(int j=0; j < 40; ++j){  
      float w = r.Gaus(mu1[j], sig1[j]);
      RatioFakediEMPtReweightedRho[k]->SetBinContent(j,w);
      float f = r1.Gaus(mu2[j], sig2[j]);
      RatioFakediEMPtReweightedRhoPrompt[k]->SetBinContent(j,f);
    } 
    RatioFakediEMPtReweightedRho[k]->Write();
    RatioFakediEMPtReweightedRhoPrompt[k]->Write();
  }

  
  
}


