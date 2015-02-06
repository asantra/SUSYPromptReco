#define PromptReco_cxx
#include "PromptReco.h"
#include <iostream>



void
photonEffectiveAreas(double _eta, double* _effA)
{
  double& effACH(_effA[0]);
  double& effANH(_effA[1]);
  double& effAPh(_effA[2]);
  double& effAW(_effA[3]);

  // Source: CutBasedPhotonID2012 twiki
  if(_eta < 1.){
    effACH = 0.012;
    effANH = 0.03;
    effAPh = 0.148;
    effAW  = 0.075;
  }
  else if(_eta < 1.479){
    effACH = 0.010;
    effANH = 0.057;
    effAPh = 0.13;
    effAW  = 0.0617;
  }
}

float dRCalc(float etaLead, float phiLead, float etaTrail, float phiTrail){
    
  float dphi = fabs(phiLead - phiTrail);
  if (dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  float deta = fabs(etaLead - etaTrail);
  float dR = sqrt(deta*deta + dphi*dphi);
  return dR;
    
}

float dPhiCalc(float phiLead, float phiTrail){
  float dphi = fabs(phiLead - phiTrail);
  if(dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  return dphi;
}

float findDiEMPt(float ELead, float EtaLead, float PhiLead, float ETrail, float EtaTrail, float PhiTrail){
  float theta1 = 2*atan(exp(-EtaLead));
  float theta2 = 2*atan(exp(-EtaTrail));
  float PX1 = ELead*sin(theta1)*cos(PhiLead);
  float PY1 = ELead*sin(theta1)*sin(PhiLead);
  float PX2 = ETrail*sin(theta2)*cos(PhiTrail);
  float PY2 = ETrail*sin(theta2)*sin(PhiTrail);
  float DiEMPt = sqrt((PX1+PX2)*(PX1+PX2)+(PY1+PY2)*(PY1+PY2));
  return DiEMPt;

}

float findMass(float ELead, float EtaLead, float PhiLead, float ETrail, float EtaTrail, float PhiTrail){
  float theta1 = 2*atan(exp(-EtaLead));
  float theta2 = 2*atan(exp(-EtaTrail));
  float PX1 = ELead*sin(theta1)*cos(PhiLead);
  float PY1 = ELead*sin(theta1)*sin(PhiLead);
  float PX2 = ETrail*sin(theta2)*cos(PhiTrail);
  float PY2 = ETrail*sin(theta2)*sin(PhiTrail);
  float PZ1 = ELead*cos(theta1);
  float PZ2 = ETrail*cos(theta2);
  float Mass = sqrt((ELead+ETrail)*(ELead+ETrail)-((PX1+PX2)*(PX1+PX2)+(PY1+PY2)*(PY1+PY2)+(PZ1+PZ2)*(PZ1+PZ2)));
  return Mass;

}

float findPt(float Energy, float Eta, float Phi){
  float theta = 2*atan(exp(-Eta));
  float PX1 = Energy*sin(theta)*cos(Phi);
  float PY1 = Energy*sin(theta)*sin(Phi);
  float PT  = sqrt(PX1*PX1+PY1*PY1);
  return PT;
}

void PromptReco::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L PromptReco.C
//      Root > PromptReco t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   TFile *f1 = new TFile("NotImportant.root","RECREATE");
   /* Defining histograms */

   ///// met histograms  //////
   /*TH1F* h_met_ffSample(new TH1F("h_met_ffSample","ff #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_eeSample(new TH1F("h_met_eeSample","ee #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_all(new TH1F("h_met_all","prompt all #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_prompt(new TH1F("h_met_prompt","prompt #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_metX_prompt(new TH1F("h_metX_prompt","prompt #slash{E}_{T} X component, Asymmetric Pt ;#slash{E}_{TX} (GeV);Events / GeV", 400, -200., 200.));
   TH1F* h_metY_prompt(new TH1F("h_metY_prompt","prompt #slash{E}_{T} Y component, Asymmetric Pt ;#slash{E}_{TY} (GeV);Events / GeV", 400, -200., 200.));
   TH1F* h_lead_photon_pt_prompt(new TH1F("h_lead_photon_pt_prompt","prompt lead photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_trail_photon_pt_prompt(new TH1F("h_trail_photon_pt_prompt","prompt trail photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_lead_photon_phi_prompt(new TH1F("h_lead_photon_phi_prompt","prompt lead photon #Phi, Asymmetric Pt ;#Phi ;Events", 80, -4., 4.));
   TH1F* h_trail_photon_phi_prompt(new TH1F("h_trail_photon_phi_prompt","prompt trail photon #Phi, Asymmetric Pt ;#Phi ;Events ", 80, -4., 4.));
   TH1F* h_lead_photon_eta_prompt(new TH1F("h_lead_photon_eta_prompt","prompt lead photon #eta, Asymmetric Pt ;#eta ;Events", 40, -2., 2.));
   TH1F* h_trail_photon_eta_prompt(new TH1F("h_trail_photon_eta_prompt","prompt trail photon #eta, Asymmetric Pt ;#eta ;Events ", 40, -2., 2.));
   TH1F* h_diEMPt_prompt(new TH1F("h_diEMPt_prompt","prompt diEMPt, Asymmetric Pt ;DiEM #P_T (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_reduced_met_prompt(new TH1F("h_reduced_met_prompt","prompt MET, Asymmetric Pt ;reduced MET (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_jetpt_prompt(new TH1F("h_jetpt_prompt","prompt jet #P_T, Asymmetric Pt ; #P_T (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_jetphi_prompt(new TH1F("h_jetphi_prompt","prompt jet #Phi, Asymmetric Pt ; #Phi;Events ", 80, -4., 4.));
   TH1F* h_jeteta_prompt(new TH1F("h_jeteta_prompt","prompt jet #eta, Asymmetric Pt ; #eta;Events ", 40, -2., 2.));*/
   
  
  /* TH1F* h_met_candidate_NoCut(new TH1F("h_met_candidate_NoCut","candidate No Cut #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_PhoIsoCut(new TH1F("h_met_candidate_PhoIsoCut","candidate with PhoIsoCut #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_PhoNeutralIso(new TH1F("h_met_candidate_PhoNeutralIso","candidate with Photon and Neutral iso #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_PhoNeutralChargedIso(new TH1F("h_met_candidate_PhoNeutralChargedIso","candidate with Photon , Charged and Neutral iso #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_PhoNeutralChargedWorstIso(new TH1F("h_met_candidate_PhoNeutralChargedWorstIso","candidate with Photon, Charged , Worst and Neutral iso #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));*/

   /////// mixing of two mets  ///////
   /*TH1F* h_met_eeSample_PurityReweighted(new TH1F("h_met_eeSample_PurityReweighted","ee Sample Purity Reweighted #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_ffSample_PurityReweighted(new TH1F("h_met_ffSample_PurityReweighted","ff Sample Purity Reweighted #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));*/
   

   ///// reduced met histograms //////
   
   /*TH1F* h_reduced_met_ffSample(new TH1F("h_reduced_met_ffSample","ff reduced #slash{E}_{T}, Asymmetric Pt ;reduced #slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_reduced_met_eeSample(new TH1F("h_reduced_met_eeSample","ee reduced #slash{E}_{T}, Asymmetric Pt ;reduced #slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));


   TH1F* h_met_eeSample_RhoReweighted(new TH1F("h_met_eeSample_RhoReweighted","ee #slash{E}_{T}, rho reweighted to candidate ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_ffSample_RhoReweighted(new TH1F("h_met_ffSample_RhoReweighted","ff #slash{E}_{T}, rho reweighted to candidate;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_eeSample_diEMPtRhoReweighted(new TH1F("h_met_eeSample_diEMPtRhoReweighted","ee #slash{E}_{T}, rho and diEMPt reweighted to candidate;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));*/
   /*TH1F* h_dR_candidate(new TH1F("h_dR_candidate","candidate #Delta R, Asymmetric Pt;#Delta R ;Events ", 100, 0., 5.));
   TH1F* h_dPhi_candidate(new TH1F("h_dPhi_candidate","candidate #Delta#phi, Asymmetric Pt ;#Delta#phi ;Events ", 80, 0., 4.));
   TH1F* h_dPhi_candidate_dr06(new TH1F("h_dPhi_candidate_dr06","candidate #Delta#phi,#Delta R > 0.6 , Asymmetric Pt;#Delta#phi ;Events ", 80, 0., 4.));*/
   // histograms for purity
   ///// candidate //////
  /* TH1F* h_et_all_passing_shower_cut_candidate(new TH1F("h_et_all_passing_shower_cut_candidate","all passing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate","lead passing trail failing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate","lead failing trail passing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate(new TH1F("h_et_all_failing_shower_cut_candidate","all failing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// fake ///////
   TH1F* h_et_all_passing_shower_cut_fake(new TH1F("h_et_all_passing_shower_cut_fake","all passing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_fake(new TH1F("h_et_lead_passing_trail_failing_shower_cut_fake","lead passing trail failing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_fake(new TH1F("h_et_lead_failing_trail_passing_shower_cut_fake","lead failing trail passing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_fake(new TH1F("h_et_all_failing_shower_cut_fake","all failing shower cut, Asymmetric Pt;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// electron //////
   TH1F* h_et_all_passing_shower_cut_electron(new TH1F("h_et_all_passing_shower_cut_electron","all passing shower cut, Asymmetric Pt;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_electron(new TH1F("h_et_lead_passing_trail_failing_shower_cut_electron","lead passing trail failing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_electron(new TH1F("h_et_lead_failing_trail_passing_shower_cut_electron","lead failing trail passing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_electron(new TH1F("h_et_all_failing_shower_cut_electron","all failing shower cut, Asymmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));

   ////// MHT histograms ////////
   TH1F* h_mht_candidate(new TH1F("h_mht_candidate","candidate MHT, Asymmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_NoCut(new TH1F("h_mht_candidate_NoCut","candidate No Cut MHT, Asymmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_PhoIsoCut(new TH1F("h_mht_candidate_PhoIsoCut","candidate with PhoIsoCut MHT, Asymmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_PhoNeutralIso(new TH1F("h_mht_candidate_PhoNeutralIso","candidate with Photon and Neutral iso MHT, Asymmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_PhoNeutralChargedIso(new TH1F("h_mht_candidate_PhoNeutralChargedIso","candidate with Photon , Charged and Neutral iso MHT, Asymmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_PhoNeutralChargedWorstIso(new TH1F("h_mht_candidate_PhoNeutralChargedWorstIso","candidate with Photon, Charged , Worst and Neutral iso MHT, Asymmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_fake(new TH1F("h_mht_fake","fake MHT, Asymmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_electron(new TH1F("h_mht_electron","electron MHT, Asymmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));

   /////// MHT Delta Phi histograms ///////
   TH1F* h_mht_deltaphi_candidate(new TH1F("h_mht_deltaphi_candidate","#Delta#Phi between MHT and nearest photon, Asymmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_NoCut(new TH1F("h_mht_deltaphi_candidate_NoCut","#Delta#Phi MHT and nearest canddiate with no cut, Asymmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_PhoIsoCut(new TH1F("h_mht_deltaphi_candidate_PhoIsoCut","#Delta#Phi MHT and candidate with PhoIsoCut, Asymmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_PhoNeutralIso(new TH1F("h_mht_deltaphi_candidate_PhoNeutralIso","#Delta#Phi MHT and candidate with Photon and Neutral iso, Asymmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_PhoNeutralChargedIso(new TH1F("h_mht_deltaphi_candidate_PhoNeutralChargedIso","#Delta#Phi MHT and candidate with Photon , Charged and Neutral iso, Asymmetric Pt ; #Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_PhoNeutralChargedWorstIso(new TH1F("h_mht_deltaphi_candidate_PhoNeutralChargedWorstIso","#Delta#Phi MHT and candidate with Photon, Charged , Worst and Neutral iso, Asymmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));


   // 2D histograms
   /*TH2F* h_leadeta_trailphi_candidate(new TH2F("h_leadeta_trailphi_candidate_HLT","lead #eta vs trail #phi, Asymmetric Pt ;lead #eta ;trail #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_leadeta_leadphi_candidate(new TH2F("h_leadeta_leadphi_candidate","lead #eta vs lead #phi, Asymmetric Pt ;lead #eta ;lead #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_traileta_leadphi_candidate(new TH2F("h_traileta_leadphi_candidate","trail #eta vs lead #phi, Asymmetric Pt;trail #eta ;lead #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_traileta_trailphi_candidate(new TH2F("h_traileta_trailphi_candidate","trail #eta vs trail #phi, Asymmetric Pt;trail #eta ;trail #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_leadeta_traileta_candidate(new TH2F("h_leadeta_traileta_candidate","lead #eta vs trail #eta, Asymmetric Pt;lead #eta ;trail #eta ", 80, -2., 2.,80,-2.,2.));
   TH2F* h_leadphi_trailphi_candidate(new TH2F("h_leadphi_trailphi_candidate","lead #phi vs trail #phi, Asymmetric Pt;lead #phi ;trail #phi ", 160, -4., 4.,160,-4.,4.));
   TH2F* h_leadpt_trailpt_candidate(new TH2F("h_leadpt_trailpt_candidate","lead Pt vs trail Pt , Asymmetric Pt;lead Pt;trail Pt ", 300, 0., 300.,300,0.,300.));*/
   // rho histograms
   /*TH1F* h_rho_ffSample(new TH1F("h_rho_ffSample","ff Rho25, Asymmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_eeSample(new TH1F("h_rho_eeSample","ee Rho25, Asymmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_candidate(new TH1F("h_rho_candidate","candidate Rho25, Asymmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));*/
   // nVtx histograms
   //TH1F* h_nVtx_ffSample(new TH1F("h_nVtx_ffSample","ff nVtx ;nVtx;Events ", 50, 0., 50.));
   //TH1F* h_nVtx_eeSample(new TH1F("h_nVtx_eeSample","ee nVtx ;nVtx;Events ", 50, 0., 50.));
   //TH1F* h_nVtx_candidate(new TH1F("h_nVtx_candidate","candidate nVtx ;nVtx;Events ", 50, 0., 50.));

   // diEMPt histograms
   /*TH1F* h_diEMPt_eeSample(new TH1F("h_diEMPt_eeSample","ee diEMPt, Asymmetric Pt Cut ;diEMPt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_diEMPt_candidate(new TH1F("h_diEMPt_candidate","candidate diEMPt, Asymmetric Pt Cut;diEMPt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_diEMPt_ffSample(new TH1F("h_diEMPt_ffSample","ff diEMPt, Asymmetric Pt Cut ;diEMPt (GeV);Events / GeV", 200, 0., 200.));*/

   // Saving the error weights //
   /*h_met_eeSample->Sumw2();
   h_met_ffSample->Sumw2();
   h_met_candidate->Sumw2();
   /*h_rho_eeSample->Sumw2();
   h_rho_ffSample->Sumw2();
   h_rho_candidate->Sumw2();
   h_diEMPt_eeSample->Sumw2();
   h_diEMPt_candidate->Sumw2();
   h_diEMPt_ffSample->Sumw2();*/
   //h_met_ffSample_RhoReweighted->Sumw2();
   //h_met_eeSample_diEMPtRhoReweighted->Sumw2();
   /*h_met_eeSample_PurityReweighted->Sumw2();
   h_met_ffSample_PurityReweighted->Sumw2();
   

   // getting the reweighting //

   TFile *Fnew = new TFile("OldHLT.root","READ");
   TH1F *h_rho_eeSample = (TH1F*)Fnew->Get("h_rho_eeSample");
   TH1F *h_rho_ffSample = (TH1F*)Fnew->Get("h_rho_ffSample");
   TH1F *h_rho_candidate = (TH1F*)Fnew->Get("h_rho_candidate");
   h_rho_eeSample->Scale(1./h_rho_eeSample->Integral(0,50));
   h_rho_ffSample->Scale(1./h_rho_ffSample->Integral(0,50));
   h_rho_candidate->Scale(1./h_rho_candidate->Integral(0,50));
   TH1F *candidateOvereeRho = new TH1F("candidateOvereeRho", "candidate  to ee Rho ratio", 50, 0., 50.);
   candidateOvereeRho->Sumw2(); 
   candidateOvereeRho->Divide(h_rho_candidate,h_rho_eeSample);
   TH1F *candidateOverffRho = new TH1F("candidateOverffRho", "candidate  to ff Rho ratio", 50, 0., 50.);
   candidateOverffRho->Sumw2();
   candidateOverffRho->Divide(h_rho_candidate,h_rho_ffSample);

   //TFile *Fnew2 = new TFile("Weighted.root","READ");
   TH1F *h_diEMPt_eeSample = (TH1F*)Fnew->Get("h_diEMPt_eeSample");
   TH1F *h_diEMPt_candidate = (TH1F*)Fnew->Get("h_diEMPt_candidate");
   h_diEMPt_eeSample->Scale(1./h_diEMPt_eeSample->Integral(0,200));
   h_diEMPt_candidate->Scale(1./h_diEMPt_candidate->Integral(0,200));
   TH1F *candidateOvereediEMPt = new TH1F("candidateOvereediEMPt", "candidate  to ee diEMPt ratio", 200, 0., 200.);
   candidateOvereediEMPt->Sumw2();
   candidateOvereediEMPt->Divide(h_diEMPt_candidate,h_diEMPt_eeSample);*/
   
    
   
   

   if (fChain == 0) return;
   
   
   
   Long64_t nentries = fChain->GetEntriesFast();
   //printf("GetEntriesFast: %lld",nentries);
   Long64_t checkentry= fChain->GetEntries();
   //printf("GetEntries: %lld",checkentry);
   
   
   /////// writing events to a file ////////
   
   
   FILE *PromptAll1     = fopen("file_All_Prompt1.txt", "w");
   FILE *PromptAll2     = fopen("file_All_Prompt2.txt", "w");
   FILE *PromptAll3     = fopen("file_All_Prompt3.txt", "w");
   FILE *PromptAll4     = fopen("file_All_Prompt4.txt", "w");
   FILE *PromptAll5     = fopen("file_All_Prompt5.txt", "w");
   FILE *PromptAll6     = fopen("file_All_Prompt6.txt", "w");
   FILE *PromptAll7     = fopen("file_All_Prompt7.txt", "w");
   FILE *PromptAll8     = fopen("file_All_Prompt8.txt", "w");
   FILE *PromptAll9     = fopen("file_All_Prompt9.txt", "w");
   FILE *PromptAll10    = fopen("file_All_Prompt10.txt", "w");
   FILE *PromptAll11    = fopen("file_All_Prompt11.txt", "w");

   FILE *PromptEle1     = fopen("file_ele_Prompt1.txt", "w");
   FILE *PromptEle2     = fopen("file_ele_Prompt2.txt", "w");
   FILE *PromptEle3     = fopen("file_ele_Prompt3.txt", "w");
   FILE *PromptEle4     = fopen("file_ele_Prompt4.txt", "w");
   FILE *PromptEle5     = fopen("file_ele_Prompt5.txt", "w");
   FILE *PromptEle6     = fopen("file_ele_Prompt6.txt", "w");
   FILE *PromptEle7     = fopen("file_ele_Prompt7.txt", "w");
   FILE *PromptEle8     = fopen("file_ele_Prompt8.txt", "w");
   FILE *PromptEle9     = fopen("file_ele_Prompt9.txt", "w");
   FILE *PromptEle10    = fopen("file_ele_Prompt10.txt", "w");
   FILE *PromptEle11    = fopen("file_ele_Prompt11.txt", "w");
   FILE *PromptEle12    = fopen("file_ele_Prompt12.txt", "w");

   ////// Reading from a file /////////
   /*FILE *myfile;
   myfile=fopen("file_ReReco.txt", "r");
   char line[1200];
   int z=0;
   Long64_t ReRecoEventNo[308440]={0};
   Long64_t ReRecoEntry[308440]={0};
   Int_t ReRecoRunNo[308440]={0};
   while(fgets(line, sizeof line, myfile) != NULL){
     fscanf(myfile,"%d\t%lld\t%lld", &ReRecoRunNo[z], &ReRecoEventNo[z], &ReRecoEntry[z]);
     //printf("%d\t%lld\t%lld\n", ReRecoRunNo[z], ReRecoEventNo[z], ReRecoEntry[z]);
     z++;
   }*/


   
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if(jentry%100000==0)printf("Entry processed: %lld\n", jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if(runNo_prompt <= 191235)fprintf(PromptAll1,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 191235 && runNo_prompt <= 192014)fprintf(PromptAll2,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 192014 && runNo_prompt <= 193572)fprintf(PromptAll3,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 193572 && runNo_prompt <= 193768)fprintf(PromptAll4,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 193768 && runNo_prompt <= 193963)fprintf(PromptAll5,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 193963 && runNo_prompt <= 194158)fprintf(PromptAll6,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 194158 && runNo_prompt <= 194353)fprintf(PromptAll7,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 194353 && runNo_prompt <= 194548)fprintf(PromptAll8,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 194548 && runNo_prompt <= 194743)fprintf(PromptAll9,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 194743 && runNo_prompt <= 194938)fprintf(PromptAll10,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
      if(runNo_prompt > 194938 && runNo_prompt <= 195130)fprintf(PromptAll11,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);

      /* /////// Photon, electron and fake Selection /////// */

      int   realphoton = 0;
      int   candidate  = 0;
      int   fakephoton = 0;
      int   fake       = 0;
      int   elephoton  = 0;
      int   electron   = 0;
      int   noCut      = 0;
      int   OnlyPho    = 0;
      int   PhoNeutral = 0;
      int   PhoNeuCharg= 0;
      //int   PhoNChargeW= 0;
      for(unsigned int i=0;i < photon_neutraliso_prompt->size();++i){
        float absEta = fabs(photon_eta_prompt->at(i));
        double effA[4];
        photonEffectiveAreas(absEta, effA);
        float theta = 2*atan(exp(-photon_eta_prompt->at(i)));
        float pt = fabs(photon_e_prompt->at(i)*sin(theta)); 

        bool neutral     = photon_neutraliso_prompt->at(i) - rho_prompt * effA[1] - 0.04 * pt < 3.5;
        bool charged     = photon_chargeiso_prompt->at(i) - rho_prompt * effA[0] < 2.6;
        bool chargedlim  = photon_chargeiso_prompt->at(i) - rho_prompt * effA[0] < 15.;
        bool photoniso   = photon_photoniso_prompt->at(i) - rho_prompt * effA[2] - 0.005 * pt < 1.3;
        //printf("absEta: %f and effACh: %f effPh:%f and effN:%f\n", absEta,effA[0],effA[2],effA[1]);
        //bool worstiso    = photon_worstiso->at(i) -rho*effA[3] < 2.6;
        bool showercut   = photon_showershape_prompt->at(i) < 0.011 ; 
        bool shower      = photon_showershape_prompt->at(i) < 0.014 ;
        bool pixelcut    = photon_pixelseed_prompt->at(i) == 0;
        bool SymPtCut    = pt > 40;
       
        if(pixelcut  )noCut++;
        if(photoniso && pixelcut )OnlyPho++;
        if(neutral && photoniso && pixelcut  )PhoNeutral++;
        if(neutral && photoniso && charged && pixelcut )PhoNeuCharg++;
        //if(neutral && photoniso && charged && pixelcut )PhoNChargeW++;
  
        if(neutral && charged && photoniso && showercut && pixelcut)   realphoton++;
        if(neutral && !charged && photoniso && showercut && pixelcut && chargedlim )  fakephoton++;
        if(neutral && charged && photoniso && showercut && !pixelcut )  elephoton++;

        if(neutral && charged && photoniso && shower && pixelcut  )  candidate++;
        if(neutral && !charged && photoniso && shower && pixelcut && chargedlim )  fake++;
        if(neutral && charged && photoniso && shower && !pixelcut )  electron++;
       
        
      }
      if(realphoton >= 2){
        //printf("Event number: %lld\n",eventNo);
        //h_met_all->Fill(met_et_prompt);
        float DR=dRCalc(photon_eta_prompt->at(0),photon_phi_prompt->at(0),photon_eta_prompt->at(1),photon_phi_prompt->at(1));
        //h_dR_candidate->Fill(DR);
        float DPHI=dPhiCalc(photon_phi_prompt->at(0),photon_phi_prompt->at(1));
        //h_dPhi_candidate->Fill(DPHI);
        if(DR > 0.6){
          /////////  Writing to a file /////////////
          
          /*if(eventNo_prompt==278137497 && runNo_prompt==195757){
            std::cout << "runno1: " << runNo_prompt << " and eventnoprompt1: " << eventNo_prompt << std::endl;
            std::cout << "photonesize: " << photon_e_prompt->size() << " and photoneprompt: " << photon_e_prompt->at(0) << std::endl;
            break;
          }*/
          //int f=0;
          /*int counter = 0;
          for(int f=0;f<308440;++f){
            if(ReRecoEventNo[f]==eventNo_prompt && ReRecoRunNo[f]==runNo_prompt){
              //printf("eventNo: %lld , runNo: %lld, ReRecoEventNo: %lld, ReRecoRunNo: %lld\n",eventNo,runNo,ReRecoEventNo[f],ReRecoRunNo[f]);
              counter++;
              break;
            }
          }
          if(counter==1){
            float DIEMPT=findDiEMPt(photon_e_prompt->at(0),photon_eta_prompt->at(0),photon_phi_prompt->at(0),photon_e_prompt->at(1),photon_eta_prompt->at(1),photon_phi_prompt->at(1));
            float PxTotal = 0;
            float PyTotal = 0;
            for(size_t k = 0; k < photon_e_prompt->size(); ++k){
              float theta = 2*atan(exp(-photon_eta_prompt->at(k)));
              float PX    = photon_e_prompt->at(k)*sin(theta)*cos(photon_phi_prompt->at(k)); 
              float PY    = photon_e_prompt->at(k)*sin(theta)*sin(photon_phi_prompt->at(k));
              PxTotal    += PX;
              PyTotal    += PY;
            }
            float redX    = met_X_prompt+PxTotal;
            float redY    = met_Y_prompt+PyTotal;
            float redMET  = sqrt(redX*redX+redY*redY);
            h_lead_photon_pt_prompt->Fill(photon_e_prompt->at(0)*2*atan(exp(-photon_eta_prompt->at(0))));
            h_trail_photon_pt_prompt->Fill(photon_e_prompt->at(1)*2*atan(exp(-photon_eta_prompt->at(1))));
            h_lead_photon_phi_prompt->Fill(photon_phi_prompt->at(0));
            h_trail_photon_phi_prompt->Fill(photon_phi_prompt->at(1));
            h_lead_photon_eta_prompt->Fill(photon_eta_prompt->at(0));
            h_trail_photon_eta_prompt->Fill(photon_eta_prompt->at(1));
            h_diEMPt_prompt->Fill(DIEMPT);
            h_met_prompt->Fill(met_et_prompt);
            h_metX_prompt->Fill(met_X_prompt);
            h_metY_prompt->Fill(met_Y_prompt);
            h_reduced_met_prompt->Fill(redMET);
            for(size_t s=0; s<jet_pt_prompt->size(); ++s){
              float dr=0;
              int count = 0;
              for(size_t j=0; j<photon_e_prompt->size(); ++j){
                dr=dRCalc(jet_eta_prompt->at(s),jet_phi_prompt->at(s),photon_eta_prompt->at(j),photon_phi_prompt->at(j));
                if(dr<0.6)count++;
              }
              if(count == 0){
                h_jetpt_prompt->Fill(jet_pt_prompt->at(s));
                h_jeteta_prompt->Fill(jet_eta_prompt->at(s));
                h_jetphi_prompt->Fill(jet_phi_prompt->at(s));
              }
            }
          }*/
        }
      }


      /*if(eventNo==705190660){
        printf("MET:                      %f \n",met_et);
        printf("MET_X:                    %f \n",met_X);
        printf("MET_Y:                    %f \n",met_Y);
        for(size_t i=0; i<photon_e->size(); ++i){
          printf("%u th photon pt :         %f \n",i,photon_e->at(i)*2*atan(exp(-photon_eta->at(i))));
          printf("%u th photon eta:         %f \n",i,photon_eta->at(i));
          printf("%u th photon phi:         %f \n",i,photon_phi->at(i));
        }
        float PxTotal = 0;
        float PyTotal = 0;
        for(size_t k = 0; k < photon_e->size(); ++k){
          float theta = 2*atan(exp(-photon_eta->at(k)));
          float PX    = photon_e->at(k)*sin(theta)*cos(photon_phi->at(k)); 
          float PY    = photon_e->at(k)*sin(theta)*sin(photon_phi->at(k));
          PxTotal    += PX;
          PyTotal    += PY;
        }
        float redX    = met_X+PxTotal;
        float redY    = met_Y+PyTotal;
        float redMET  = sqrt(redX*redX+redY*redY);
        float DIEMPT=findDiEMPt(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
        printf("diEMPt of photons:        %f \n",DIEMPT);
        printf("reduced MET:              %f \n", redMET);
        
        for(size_t j=0; j<jet_pt->size(); ++j){
          printf("corrected %u th jet pt:   %f \n",j,jet_pt->at(j));
          printf("corrected %u th jet eta:  %f \n",j,jet_eta->at(j));
          printf("corrected %u th jet phi:  %f \n",j,jet_phi->at(j));

        }
      }*/
      if(fakephoton >= 2){
        float DR=dRCalc(photon_eta_prompt->at(0),photon_phi_prompt->at(0),photon_eta_prompt->at(1),photon_phi_prompt->at(1));
        if(DR > 0.6){
          
          /*float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j = 0; j < photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          float redX    = met_X+PxTotal;
          float redY    = met_Y+PyTotal;
          float redMET  = sqrt(redX*redX+redY*redY);
          //h_reduced_met_ffSample->Fill(redMET);
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            //h_mht_fake->Fill(MHT);
          }
          float DIEMPT=findDiEMPt(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
          //h_diEMPt_ffSample->Fill(DIEMPT);
          //h_rho_ffSample->Fill(rho);
          /*int Bin=candidateOverffRho->FindBin(rho);
          float Weight=candidateOverffRho->GetBinContent(Bin);
          h_met_ffSample_RhoReweighted->Fill(met_et, Weight);
          h_met_ffSample_PurityReweighted->Fill(met_et, Weight*0.48459);*/
          
        }
      }
      if(elephoton >= 2){
        float DR=dRCalc(photon_eta_prompt->at(0),photon_phi_prompt->at(0),photon_eta_prompt->at(1),photon_phi_prompt->at(1));
        float MASS=findMass(photon_e_prompt->at(0),photon_eta_prompt->at(0),photon_phi_prompt->at(0),photon_e_prompt->at(1),photon_eta_prompt->at(1),photon_phi_prompt->at(1));
        if(DR > 0.6 && MASS>=75 && MASS<=105 ){

          if(runNo_prompt <= 193572)fprintf(PromptEle1,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 193572 && runNo_prompt <= 195130)fprintf(PromptEle2,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 195130 && runNo_prompt <= 196688)fprintf(PromptEle3,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 196688 && runNo_prompt <= 198246)fprintf(PromptEle4,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 198246 && runNo_prompt <= 199804)fprintf(PromptEle5,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 199804 && runNo_prompt <= 201362)fprintf(PromptEle6,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 201362 && runNo_prompt <= 202920)fprintf(PromptEle7,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 202920 && runNo_prompt <= 204478)fprintf(PromptEle8,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 204478 && runNo_prompt <= 206036)fprintf(PromptEle9,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 206036 && runNo_prompt <= 206815)fprintf(PromptEle10,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 206815 && runNo_prompt <= 207594)fprintf(PromptEle11,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          if(runNo_prompt > 207594 && runNo_prompt <= 209152)fprintf(PromptEle12,"%d\t%lld\t%lld\n",runNo_prompt,eventNo_prompt,jentry);
          
          //printf("elephoton Entry: %lld\n", jentry);
          //h_met_eeSample->Fill(met_et);
          
          /*float DIEMPT=findDiEMPt(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          float redX    = met_X+PxTotal;
          float redY    = met_Y+PyTotal;
          float redMET  = sqrt(redX*redX+redY*redY);
          //h_reduced_met_eeSample->Fill(redMET);
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            } 
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            //h_mht_electron->Fill(MHT);
          }
          /*int Bin=candidateOvereeRho->FindBin(rho);
          float Weight=candidateOvereeRho->GetBinContent(Bin);
          h_met_eeSample_RhoReweighted->Fill(met_et, Weight);
          int Bin2=candidateOvereediEMPt->FindBin(DIEMPT);
          float Weight2=candidateOvereediEMPt->GetBinContent(Bin2);
          h_met_eeSample_diEMPtRhoReweighted->Fill(met_et, Weight*Weight2);
          h_met_eeSample_PurityReweighted->Fill(met_et, 0.51541*Weight*Weight2);*/
          //h_diEMPt_eeSample->Fill(DIEMPT);
          //h_rho_eeSample->Fill(rho);
        }
      }
      // calculating purity from matrix method
      /*if(candidate >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR>0.6){
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)<0.011){
            //h_et_all_passing_shower_cut_candidate->Fill(photon_e->at(0)); h_met_candidate->Fill(met_et); 
            float PxTotal = 0;
            float PyTotal = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              float theta = 2*atan(exp(-photon_eta->at(j)));
              float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
              float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
              PxTotal    += PX;
              //printf("PxTotal: %f and PX: %f for entry: %lld\n", PxTotal, PX, jentry);
              PyTotal    += PY;
              //printf("PyTotal: %f and PY: %f for entry: %lld\n", PyTotal, PY, jentry);
            }
            for(size_t k=0; k<jet_pt->size(); ++k){
              float dr=0;
              int counter = 0;
              for(size_t j=0; j<photon_e->size(); ++j){
                dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
                if(dr<0.6)counter++;
              }
              if(counter == 0){
                float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
                float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
                PxTotal    += PXj;
                //printf("PxTotal: %f and PXj: %f for entry: %lld\n", PxTotal, PXj, jentry);
                PyTotal    += PYj;
                //printf("PyTotal: %f and PYj: %f for entry: %lld\n", PyTotal, PYj, jentry);
              }
            }
            if(PxTotal!=0 && PyTotal!=0){
              float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
              //printf("final PxTotal: %f PyTotal: %f for entry: %lld and MHT: %f\n",PxTotal, PyTotal , jentry, MHT);
              //h_mht_candidate->Fill(MHT); 
              float MHTPhi=atan2(-PyTotal,-PxTotal);
              float deltaphi[10]={0};
              for(unsigned int i=0; i< photon_e->size();++i){
                deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
                if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
              }
              float mindphi=deltaphi[0];
              for(unsigned int i=0;i<photon_e->size();++i){
                if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
              }
              //h_mht_deltaphi_candidate->Fill(mindphi);
               
            }
            float DIEMPT=findDiEMPt(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
            //h_diEMPt_candidate->Fill(DIEMPT);
            //h_rho_candidate->Fill(rho); 
          }
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)>0.011)h_et_lead_passing_trail_failing_shower_cut_candidate->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)<0.011)h_et_lead_failing_trail_passing_shower_cut_candidate->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)>0.011)h_et_all_failing_shower_cut_candidate->Fill(photon_e->at(0));
        }
      }

      if(fake >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR>0.6){
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)<0.011)h_et_all_passing_shower_cut_fake->Fill(photon_e->at(0));
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)>0.011)h_et_lead_passing_trail_failing_shower_cut_fake->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)<0.011)h_et_lead_failing_trail_passing_shower_cut_fake->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)>0.011)h_et_all_failing_shower_cut_fake->Fill(photon_e->at(0));
        }
      }

      if(electron >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        float MASS=findMass(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
        //printf("MASS:%f for entry:%lld\n",MASS,jentry);
        if(DR>0.6 && MASS>=75 && MASS<=105){
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)<0.011)h_et_all_passing_shower_cut_electron->Fill(photon_e->at(0));
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)>0.011)h_et_lead_passing_trail_failing_shower_cut_electron->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)<0.011)h_et_lead_failing_trail_passing_shower_cut_electron->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)>0.011)h_et_all_failing_shower_cut_electron->Fill(photon_e->at(0));
        }
      }

      if(noCut >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            //if(counter > 0)printf("counter: %d for entry:%lld\n ",counter, jentry);
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_NoCut->Fill(MHT);
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(unsigned int i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(unsigned int i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_NoCut->Fill(mindphi);  
          }
          h_met_candidate_NoCut->Fill(met_et);
          
        }
      }
      if(OnlyPho >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_PhoIsoCut->Fill(MHT);
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(unsigned int i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(unsigned int i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_PhoIsoCut->Fill(mindphi);  
          }
          h_met_candidate_PhoIsoCut->Fill(met_et);
        }
      }
      if(PhoNeutral >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_PhoNeutralIso->Fill(MHT); 
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(unsigned int i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(unsigned int i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_PhoNeutralIso->Fill(mindphi); 
          }
          h_met_candidate_PhoNeutralIso->Fill(met_et);
        }
      }
      if(PhoNeuCharg >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj; 
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_PhoNeutralChargedIso->Fill(MHT);
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(unsigned int i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(unsigned int i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_PhoNeutralChargedIso->Fill(mindphi);  
          }
          h_met_candidate_PhoNeutralChargedIso->Fill(met_et);
        }
      }
      if(PhoNChargeW >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            } 
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_PhoNeutralChargedWorstIso->Fill(MHT);
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(unsigned int i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(unsigned int i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_PhoNeutralChargedWorstIso->Fill(mindphi);  
          }
          h_met_candidate_PhoNeutralChargedWorstIso->Fill(met_et);
        }
      }*/

      // if (Cut(ientry) < 0) continue;*/
     
      
   }
   
   //printf("total same events:%d\n",counter);
    
   
   fclose(PromptAll1);
   fclose(PromptAll2);
   fclose(PromptAll3);
   fclose(PromptAll4);
   fclose(PromptAll5);
   fclose(PromptAll6);
   fclose(PromptAll7);
   fclose(PromptAll8);
   fclose(PromptAll9);
   fclose(PromptAll10);
   fclose(PromptAll11);
   
   fclose(PromptEle1);
   fclose(PromptEle2);
   fclose(PromptEle3);
   fclose(PromptEle4);
   fclose(PromptEle5);
   fclose(PromptEle6);
   fclose(PromptEle7);
   fclose(PromptEle8);
   fclose(PromptEle9);
   fclose(PromptEle10);
   fclose(PromptEle11);
   fclose(PromptEle12);
   //fclose(myfile);
   f1->cd();
   f1->Write();
   
}
                      
