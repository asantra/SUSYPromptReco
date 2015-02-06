
{  

  //////// electron reweighting  ////////////
   TFile f8("PromptRecoElectronCorPfMet.root");
   h_rho_rereco->Sumw2();
   Float_t y1=h_rho_rereco->Integral(0,40);
   h_rho_rereco->Scale(1.0/y1);
   TH1F *ElectronRho_rereco=(TH1F*)h_rho_rereco->Clone("ElectronRho_rereco");

   h_rho_prompt->Sumw2();
   Float_t y12=h_rho_prompt->Integral(0,40);
   h_rho_prompt->Scale(1.0/y12);
   TH1F *ElectronRho_prompt=(TH1F*)h_rho_prompt->Clone("ElectronRho_prompt");

   h_rho_rereco_diEMPtreweighted->Sumw2();
   Float_t y13=h_rho_rereco_diEMPtreweighted->Integral(0,40);
   h_rho_rereco_diEMPtreweighted->Scale(1.0/y13);
   TH1F *ElectrondiEMPtreweightedRho_rereco=(TH1F*)h_rho_rereco_diEMPtreweighted->Clone("ElectrondiEMPtreweightedRho_rereco");

   h_rho_prompt_diEMPtreweighted->Sumw2();
   Float_t y14=h_rho_prompt_diEMPtreweighted->Integral(0,40);
   h_rho_prompt_diEMPtreweighted->Scale(1.0/y14);
   TH1F *ElectrondiEMPtreweightedRho_prompt=(TH1F*)h_rho_prompt_diEMPtreweighted->Clone("ElectrondiEMPtreweightedRho_prompt");

   h_rho_rereco_ffdiEMPtreweighted->Sumw2();
   Float_t y15=h_rho_rereco_ffdiEMPtreweighted->Integral(0,40);
   h_rho_rereco_ffdiEMPtreweighted->Scale(1.0/y15);
   TH1F *ElectronffdiEMPtreweightedRho_rereco=(TH1F*)h_rho_rereco_ffdiEMPtreweighted->Clone("ElectronffdiEMPtreweightedRho_rereco");

   h_rho_prompt_ffdiEMPtreweighted->Sumw2();
   Float_t y16=h_rho_prompt_ffdiEMPtreweighted->Integral(0,40);
   h_rho_prompt_ffdiEMPtreweighted->Scale(1.0/y16);
   TH1F *ElectronffdiEMPtreweightedRho_prompt=(TH1F*)h_rho_prompt_ffdiEMPtreweighted->Clone("ElectronffdiEMPtreweightedRho_prompt");


   h_diEMPt_rereco->Sumw2();
   Float_t y2=h_diEMPt_rereco->Integral(0,250);
   h_diEMPt_rereco->Scale(1.0/y2);
   TH1F *ElectronDiEMPt_rereco=(TH1F*)h_diEMPt_rereco->Clone("ElectronDiEMPt_rereco");

   h_diEMPt_prompt->Sumw2();
   Float_t y22=h_diEMPt_prompt->Integral(0,250);
   h_diEMPt_prompt->Scale(1.0/y22);
   TH1F *ElectronDiEMPt_prompt=(TH1F*)h_diEMPt_prompt->Clone("ElectronDiEMPt_prompt");





 //////// fake reweighting //////////////

   TFile f9("PromptRecoFakeCorPfMet.root");
   h_rho_rereco->Sumw2();
   Float_t z1=h_rho_rereco->Integral(0,40);
   h_rho_rereco->Scale(1.0/z1);
   TH1F *FakeRho_rereco=(TH1F*)h_rho_rereco->Clone("FakeRho_rereco");

   h_rho_prompt->Sumw2();
   Float_t z12=h_rho_prompt->Integral(0,40);
   h_rho_prompt->Scale(1.0/z12);
   TH1F *FakeRho_prompt=(TH1F*)h_rho_prompt->Clone("FakeRho_prompt");

   h_diEMPt_rereco->Sumw2();
   Float_t z2=h_diEMPt_rereco->Integral(0,250);
   h_diEMPt_rereco->Scale(1.0/z2);
   TH1F *FakeDiEMPt_rereco=(TH1F*)h_diEMPt_rereco->Clone("FakeDiEMPt_rereco");

   h_diEMPt_prompt->Sumw2();
   Float_t z22=h_diEMPt_prompt->Integral(0,250);
   h_diEMPt_prompt->Scale(1.0/z22);
   TH1F *FakeDiEMPt_prompt=(TH1F*)h_diEMPt_prompt->Clone("FakeDiEMPt_prompt");


  ///////// photon reweighting ///////////
   
  
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



   ///// division of histograms //////////

   TH1F *ElectronDiEMPt = new TH1F("ElectronDiEMPt","photon over electron DiEMPt",250,0,250);
   ElectronDiEMPt->Sumw2();
   ElectronDiEMPt->Divide(PhotonDiEMPt_rereco,ElectronDiEMPt_rereco);
   //ElectronDiEMPt->Draw();

   TH1F *ElectronDiEMPtPrompt = new TH1F("ElectronDiEMPtPrompt","photon over electron DiEMPt",250,0,250);
   ElectronDiEMPtPrompt->Sumw2();
   ElectronDiEMPtPrompt->Divide(PhotonDiEMPt_prompt,ElectronDiEMPt_prompt);
   //ElectronDiEMPtPrompt->Draw();

   TH1F *ElectronRho = new TH1F("ElectronRho","photon over electron rho",40,0,40);
   ElectronRho->Sumw2();
   ElectronRho->Divide(PhotonRho_rereco,ElectronRho_rereco);
   //ElectronRho->Draw("sames");

   TH1F *ElectronRhoPrompt = new TH1F("ElectronRhoPrompt","photon over electron rho",40,0,40);
   ElectronRhoPrompt->Sumw2();
   ElectronRhoPrompt->Divide(PhotonRho_prompt,ElectronRho_prompt);

   TH1F *ElectrondiEMPtReweightedRho = new TH1F("ElectrondiEMPtReweightedRho","photon over diEMPt reweighted to candidate electron rho",40,0,40);
   ElectrondiEMPtReweightedRho->Sumw2();
   ElectrondiEMPtReweightedRho->Divide(PhotonRho_rereco,ElectrondiEMPtreweightedRho_rereco);
   //ElectronRho->Draw("sames");

   TH1F *ElectrondiEMPtReweightedRhoPrompt = new TH1F("ElectrondiEMPtReweightedRhoPrompt","photon over diEMPt Reweighted to candidate electron rho",40,0,40);
   ElectrondiEMPtReweightedRhoPrompt->Sumw2();
   ElectrondiEMPtReweightedRhoPrompt->Divide(PhotonRho_prompt,ElectrondiEMPtreweightedRho_prompt);
   //ElectronRhoPrompt->Draw("sames");


   TH1F *FakeElectronDiEMPt = new TH1F("FakeElectronDiEMPt","fake over electron DiEMPt",250,0,250);
   FakeElectronDiEMPt->Sumw2();
   FakeElectronDiEMPt->Divide(FakeDiEMPt_rereco,ElectronDiEMPt_rereco);
   //ElectronDiEMPt->Draw();

   TH1F *FakeElectronDiEMPtPrompt = new TH1F("FakeElectronDiEMPtPrompt","fake over electron DiEMPt",250,0,250);
   FakeElectronDiEMPtPrompt->Sumw2();
   FakeElectronDiEMPtPrompt->Divide(FakeDiEMPt_prompt,ElectronDiEMPt_prompt);
   //FakeElectronDiEMPtPrompt->Draw();

   TH1F *FakeElectronRho = new TH1F("FakeElectronRho","fake over electron rho",40,0,40);
   FakeElectronRho->Sumw2();
   FakeElectronRho->Divide(FakeRho_rereco,ElectronRho_rereco);
   //ElectronRho->Draw("sames");

   TH1F *FakeElectronRhoPrompt = new TH1F("FakeElectronRhoPrompt","fake over electron rho",40,0,40);
   FakeElectronRhoPrompt->Sumw2();
   FakeElectronRhoPrompt->Divide(FakeRho_prompt,ElectronRho_prompt);



   TH1F *FakeElectrondiEMPtReweightedRho = new TH1F("FakeElectrondiEMPtReweightedRho","fake over diEMPt Reweighted to ff electron rho",40,0,40);
   FakeElectrondiEMPtReweightedRho->Sumw2();
   FakeElectrondiEMPtReweightedRho->Divide(FakeRho_rereco,ElectronffdiEMPtreweightedRho_rereco);
   //ElectronRho->Draw("sames");

   TH1F *FakeElectrondiEMPtReweightedRhoPrompt = new TH1F("FakeElectrondiEMPtReweightedRhoPrompt","fake over diEMPt Reweighted to ff electron rho",40,0,40);
   FakeElectrondiEMPtReweightedRhoPrompt->Sumw2();
   FakeElectrondiEMPtReweightedRhoPrompt->Divide(FakeRho_prompt,ElectronffdiEMPtreweightedRho_prompt);
   //FakeElectronRhoPrompt->Draw("sames");



   TFile* fout = new TFile("ElectronRatio.root","RECREATE");
   fout->cd();
   ElectronRho->Write();
   ElectronDiEMPt->Write();
   ElectronRhoPrompt->Write();
   ElectronDiEMPtPrompt->Write();
   ElectrondiEMPtReweightedRho->Write();
   ElectrondiEMPtReweightedRhoPrompt->Write();
   FakeElectronRho->Write();
   FakeElectronDiEMPt->Write();
   FakeElectronRhoPrompt->Write();
   FakeElectronDiEMPtPrompt->Write();
   FakeElectrondiEMPtReweightedRho->Write();
   FakeElectrondiEMPtReweightedRhoPrompt->Write();
   
   
   

}
