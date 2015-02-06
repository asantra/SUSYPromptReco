
{  
   TFile f2("PromptRecoFakePercent.root");
   Float_t x1=h_trail_photon_pt_percentdifference1->Integral();
   h_trail_photon_pt_percentdifference1->Sumw2();
   h_trail_photon_pt_percentdifference1->Scale(1.0/x1);
   h_trail_photon_pt_percentdifference1->SetLineColor(2);
   TH1F *Fake=(TH1F*)h_trail_photon_pt_percentdifference1->Clone("Fake");
   //Fake->GetXaxis()->SetRange(-0.5,0.5);
   //h_trail_photon_pt_percentdifference1->GetXaxis()->SetRange(0,75);
   std::cout << "Fake:" << Fake->Integral() << std::endl;
   
  
   TFile f3("PromptRecoElectronPercentAll.root");
   Float_t x2=h_trail_photon_pt_percentdifference1->Integral();
   h_trail_photon_pt_percentdifference1->Sumw2();
   h_trail_photon_pt_percentdifference1->Scale(1.0/x2);
   h_trail_photon_pt_percentdifference1->SetLineColor(4);
   TH1F *Electron=(TH1F*)h_trail_photon_pt_percentdifference1->Clone("Electron");
   std::cout << "Electron:" << Electron->Integral() << std::endl;
   Electron->GetXaxis()->SetRange(-0.5,0.5);
   //h_trail_photon_pt_rereco->GetXaxis()->SetRange(0,75);
   

   TCanvas *Q = new TCanvas("Q","Comparison",1200,900);
   //Q->Divide(2,1);
   /*Q->Modified(); Q->Update();
   TPaveStats *stats =(TPaveStats*)Q->GetPrimitive("stats");
   stats->SetName("h_trail_photon_pt_percentdifference1");
   stats->SetX1NDC(.4);
   stats->SetX2NDC(.6);
   stats->SetTextColor(4);

   Q->Update();
   TPaveStats *stats2 = (TPaveStats*)Q->GetPrimitive("stats");
   stats2->SetName("h_trail_photon_pt_rereco");
   stats2->SetX1NDC(.1);
   stats2->SetX2NDC(.3);
   stats2->SetTextColor(2);*/
   
   /*TPaveStats *st1 = (TPaveStats*)h_trail_photon_pt_rereco->GetListOfFunctions()->FindObject("stats"); 
   TPaveStats *st2 = (TPaveStats*)h_trail_photon_pt_percentdifference1->GetListOfFunctions()->FindObject("stats"); 
   st1->SetX1NDC(.5); st1->SetX2NDC(.7); 
   st2->SetX1NDC(.2); st2->SetX2NDC(.4); */
   
   gStyle->SetStatTextColor(2);
   Fake->Draw();
   gStyle->SetOptStat(111111);


   gStyle->SetStatTextColor(4);
   Electron->Draw("sames");
   gStyle->SetOptStat(111111);

   
   
   
    //to for the generation of the 'stat" boxes 
   // to draw stat boxes in different positions
   
   Q->Modified(); 



   ///// not needed for now  /////
   /**/
   
   
   //h_nvtx_candidate->Sumw2();
   /*TFile f3("Data_PromptReco_2012_PileUp.root");
   Float_t x2=pileup->Integral();
   cout << x2 << endl;
   TH1F *pileup=(TH1F*)pileup->Clone("pileup");
   
   pileup->Scale(1.0/x2);
   pileup->SetLineColor(2);
   //h_nvtx_weight->Sumw2();
   ReweightedMCPileUp->Draw();
   pileup->Draw("sames");*/

}
