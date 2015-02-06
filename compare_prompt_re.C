
{  
   

   TFile f8("PromptRecoElectronCorPfMet.root");
   h_met_prompt_ffdiEMPtrhoreweighted->Sumw2();
   h_met_prompt_ffdiEMPtrhoreweighted->Rebin(4);
   Float_t x2=h_met_prompt_ffdiEMPtrhoreweighted->Integral(0,50);
   h_met_prompt_ffdiEMPtrhoreweighted->Scale(1.0/x2);
   h_met_prompt_ffdiEMPtrhoreweighted->SetLineColor(4);
   h_met_prompt_ffdiEMPtrhoreweighted->GetXaxis()->SetRangeUser(0,80);
   TH1F *Electron_prompt_ffdiEMPtrhoreweighted=(TH1F*)h_met_prompt_ffdiEMPtrhoreweighted->Clone("Electron_prompt_ffdiEMPtrhoreweighted");
   
  
   TFile f10("PromptRecoFakeCorPfMet.root");
   h_met_prompt->Sumw2();
   h_met_prompt->Rebin(4);
   Float_t x1=h_met_prompt->Integral(0,50);
   h_met_prompt->Scale(1.0/x1);
   h_met_prompt->SetLineColor(2);
   h_met_prompt->GetXaxis()->SetRangeUser(0,80);
   TH1F *Fake_prompt=(TH1F*)h_met_prompt->Clone("Fake_prompt");
   

   TCanvas *Q = new TCanvas("Q","Comparison",600,450);//1200,900
   
   Q->cd();
   
   TH1F *h2 = new TH1F("h2","electron over fake",50,0,50);
   h2->Sumw2();
   h2->Divide(Electron_prompt_ffdiEMPtrhoreweighted,Fake_prompt);
   //Electron_prompt->SetStats(0);
   h2->GetXaxis()->SetRange(0,20);
   
   h2->GetYaxis()->SetRangeUser(0.2,2.0);
   h2->Draw();
   
   
   Q->Modified(); 


}
