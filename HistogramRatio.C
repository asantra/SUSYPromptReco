 {
   TFile f1("PromptRecoFakeRawPfMet.root");
   h_met_rereco->Sumw2();
   h_met_rereco->Rebin(4);
   Float_t x1=h_met_rereco->Integral(0,50);
   h_met_rereco->Scale(1.0/x1);
   h_met_rereco->SetLineColor(2);
   h_met_rereco->GetXaxis()->SetRangeUser(0,80);
   TH1F *Fake_rereco=(TH1F*)h_met_rereco->Clone("Fake_rereco");
  
   TFile f3("PromptRecoPhotonRawPfMet.root");
   h_met_rereco->Sumw2();
   h_met_rereco->Rebin(4);
   Float_t x11=h_met_rereco->Integral(0,50);
   h_met_rereco->Scale(1.0/x11);
   h_met_rereco->SetLineColor(kGreen+3);
   h_met_rereco->GetXaxis()->SetRangeUser(0,80);
   TH1F *Photon_rereco=(TH1F*)h_met_rereco->Clone("Photon_rereco");

   //Fake_rereco->GetXaxis()->SetRange(0,100); //font in pixels
   /*Fake_rereco->GetXaxis()->SetLabelSize(16); //in pixels
   Fake_rereco->GetYaxis()->SetLabelFont(63); //font in pixels
   Fake_rereco->GetYaxis()->SetLabelSize(16); //in pixels*/
   
   TFile f2("PromptRecoElectronRawPfMet.root");
   h_met_rereco->Sumw2();
   h_met_rereco->Rebin(4);
   Float_t x2=h_met_rereco->Integral(0,50);
   h_met_rereco->Scale(1.0/x2);
   h_met_rereco->SetLineColor(4);
   h_met_rereco->GetXaxis()->SetRangeUser(0,80);
   TH1F *Electron_rereco=(TH1F*)h_met_rereco->Clone("Electron_rereco");
   //Electron_rereco->GetXaxis()->SetRange(0,100);

   
   TCanvas *Q = new TCanvas("Q","Comparison",900,900); //1200,900
   TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad1->Draw();
   pad1->cd();
   
   gStyle->SetStatTextColor(2);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2); 
   Fake_rereco->Draw();
   //Fakerhoreweighted_rereco->Draw();
   gStyle->SetOptStat(111111);
   
   

   gStyle->SetStatTextColor(kGreen+3);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   Photon_rereco->Draw("sames");
   gStyle->SetOptStat(111111);
   
   
   
   
   Q->cd();
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad2->Draw();
   pad2->cd();
   
   //Fake_rereco->SetStats(0);
   //Electron_rereco->SetStats(0);
   TH1F *h2 = new TH1F("h2","photon over fake",50,0,50);
   h2->Sumw2();
   
   //h2->SetStats(0);
   h2->Divide(Photon_rereco, Fake_rereco);
   //gStyle->SetStatTextColor(kGreen+3);
   
   h2->GetXaxis()->SetRangeUser(0,20);
   h2->GetYaxis()->SetRangeUser(0.2,2.0);
   //h2->SetMarkerStyle(21);
   h2->Draw("ep");
   
   Q->SaveAs("CandidateFakeReRecoRawPfMetRatio.eps");

   TCanvas *R = new TCanvas("R","Ratio",600,450); //1200,900
   R->cd();
   h2->Draw();


   //

 
}
       
