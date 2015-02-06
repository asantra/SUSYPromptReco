
{  
   TFile f8("PromptRecoFakeCorPfMet.root");
   h_rho_rereco->Sumw2();
   Float_t y1=h_rho_rereco->Integral(0,40);
   h_rho_rereco->Scale(1.0/y1);
   TH1F *FakeRho_rereco=(TH1F*)h_rho_rereco->Clone("FakeRho_rereco");

   h_rho_prompt->Sumw2();
   Float_t y12=h_rho_prompt->Integral(0,40);
   h_rho_prompt->Scale(1.0/y12);
   TH1F *FakeRho_prompt=(TH1F*)h_rho_prompt->Clone("FakeRho_prompt");

   h_rho_rereco_diEMPtreweighted->Sumw2();
   Float_t y41=h_rho_rereco_diEMPtreweighted->Integral(0,40);
   h_rho_rereco_diEMPtreweighted->Scale(1.0/y41);
   TH1F *FakediEMPtReweightedRho_rereco=(TH1F*)h_rho_rereco_diEMPtreweighted->Clone("FakediEMPtReweightedRho_rereco");

   h_rho_prompt_diEMPtreweighted->Sumw2();
   Float_t y42=h_rho_prompt_diEMPtreweighted->Integral(0,40);
   h_rho_prompt_diEMPtreweighted->Scale(1.0/y42);
   TH1F *FakediEMPtReweightedRho_prompt=(TH1F*)h_rho_prompt_diEMPtreweighted->Clone("FakediEMPtReweightedRho_prompt");

   h_diEMPt_rereco->Sumw2();
   Float_t y2=h_diEMPt_rereco->Integral(0,250);
   h_diEMPt_rereco->Scale(1.0/y2);
   TH1F *FakeDiEMPt_rereco=(TH1F*)h_diEMPt_rereco->Clone("FakeDiEMPt_rereco");

   h_diEMPt_prompt->Sumw2();
   Float_t y22=h_diEMPt_prompt->Integral(0,250);
   h_diEMPt_prompt->Scale(1.0/y22);
   TH1F *FakeDiEMPt_prompt=(TH1F*)h_diEMPt_prompt->Clone("FakeDiEMPt_prompt");
   
  
   TFile f10("PromptRecoPhotonCorPfMet.root");
   h_rho_rereco->Sumw2();
   Float_t x2=h_rho_rereco->Integral(0,40);
   h_rho_rereco->Scale(1.0/x2);
   TH1F *PhotonRho_rereco=(TH1F*)h_rho_rereco->Clone("PhotonRho_rereco");

   h_rho_prompt->Sumw2();
   Float_t x22=h_rho_prompt->Integral(0,40);
   h_rho_prompt->Scale(1.0/x22);
   TH1F *PhotonRho_prompt=(TH1F*)h_rho_prompt->Clone("PhotonRho_prompt");

   h_diEMPt_rereco->Sumw2();
   Float_t x3=h_diEMPt_rereco->Integral(0,250);
   h_diEMPt_rereco->Scale(1.0/x3);
   TH1F *PhotonDiEMPt_rereco=(TH1F*)h_diEMPt_rereco->Clone("PhotonDiEMPt_rereco");

   h_diEMPt_prompt->Sumw2();
   Float_t x32=h_diEMPt_prompt->Integral(0,250);
   h_diEMPt_prompt->Scale(1.0/x32);
   TH1F *PhotonDiEMPt_prompt=(TH1F*)h_diEMPt_prompt->Clone("PhotonDiEMPt_prompt");

   TH1F *FakeDiEMPt = new TH1F("FakeDiEMPt","photon over fake DiEMPt",250,0,250);
   FakeDiEMPt->Sumw2();
   FakeDiEMPt->Divide(PhotonDiEMPt_rereco,FakeDiEMPt_rereco);
   

   TH1F *FakeDiEMPtPrompt = new TH1F("FakeDiEMPtPrompt","photon over fake DiEMPt",250,0,250);
   FakeDiEMPtPrompt->Sumw2();
   FakeDiEMPtPrompt->Divide(PhotonDiEMPt_prompt,FakeDiEMPt_prompt);
   

   TH1F *FakeRho = new TH1F("FakeRho","photon over electron rho",40,0,40);
   FakeRho->Sumw2();
   FakeRho->Divide(PhotonRho_rereco,FakeRho_rereco);

   TH1F *FakeDiEMPtReweightedRho = new TH1F("FakeDiEMPtReweightedRho","photon over fake rho",40,0,40);
   FakeDiEMPtReweightedRho->Sumw2();
   FakeDiEMPtReweightedRho->Divide(PhotonRho_rereco,FakediEMPtReweightedRho_rereco);
   

   TH1F *FakeRhoPrompt = new TH1F("FakeRhoPrompt","photon over fake rho",40,0,40);
   FakeRhoPrompt->Sumw2();
   FakeRhoPrompt->Divide(PhotonRho_prompt,FakeRho_prompt);
   

   TH1F *FakeDiEMPtReweightedRhoPrompt = new TH1F("FakeDiEMPtReweightedRhoPrompt","photon over fake rho",40,0,40);
   FakeDiEMPtReweightedRhoPrompt->Sumw2();
   FakeDiEMPtReweightedRhoPrompt->Divide(PhotonRho_prompt,FakediEMPtReweightedRho_prompt);
   

   TFile* fout = new TFile("FakeRatio.root","RECREATE");
   fout->cd();
   FakeRho->Write();
   FakeDiEMPt->Write();
   FakeRhoPrompt->Write();
   FakeDiEMPtPrompt->Write();
   FakeDiEMPtReweightedRho->Write();
   FakeDiEMPtReweightedRhoPrompt->Write();
   
   

}
