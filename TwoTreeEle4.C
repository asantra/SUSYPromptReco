#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.
#include <vector>

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

void TwoTreeEle4(){

  ////// get trees  ///////
  TChain b1a("tree");
  TChain b1b("tree");
  TChain b2a("NNtree");
  TChain b2b("NNtree");
  b1a.Add("PhotonAll_PfMet_ReReco.root"); // raw PfMet from rereco
  b1b.Add("PhotonAll_PfMet_ReReco.root");
  b2a.Add("PhotonAll_Prompt.root");
  b2b.Add("PhotonAll_Prompt.root");

  ///// variable declaration for tree 1 ///////
  UInt_t          runno;
  Int_t           lumino;
  ULong64_t       eventno;
  Int_t           verticesn;
  Int_t           photonsn;
  Float_t         metet;
  Float_t         metX;
  Float_t         metY;
  Float_t         rho;
  vector<float>   *photone = 0;
  vector<float>   *photoneta = 0;
  vector<float>   *photonphi = 0;
  vector<float>   *photonshowershape = 0;
  vector<int>     *photonpixelseed = 0;
  vector<float>   *photonchargeiso = 0;
  vector<float>   *photonneutraliso = 0;
  vector<float>   *photonphotoniso = 0;
  vector<float>   *jetpt = 0;
  //vector<float>   *jetptuncorrected = 0;
  vector<float>   *jeteta = 0;
  vector<float>   *jetphi = 0;
  //vector<float>   *jetetauncorrected = 0;
  //vector<float>   *jetphiuncorrected = 0;
  b1a.SetBranchAddress("runNo",&runno);
  b1a.SetBranchAddress("eventNo",&eventno);
  b1b.SetBranchAddress("photon_e", &photone);
  b1b.SetBranchAddress("lumiNo",&lumino);
  b1b.SetBranchAddress("eventNo",&eventno);
  b1b.SetBranchAddress("vertices_n",&verticesn);
  b1b.SetBranchAddress("photons_n",&photonsn);
  b1b.SetBranchAddress("met_et",&metet);
  b1b.SetBranchAddress("met_X",&metX);
  b1b.SetBranchAddress("met_Y",&metY);
  b1b.SetBranchAddress("rho", &rho);
  b1b.SetBranchAddress("photon_eta", &photoneta);
  b1b.SetBranchAddress("photon_phi", &photonphi);
  b1b.SetBranchAddress("photon_showershape", &photonshowershape);
  b1b.SetBranchAddress("photon_pixelseed", &photonpixelseed);
  b1b.SetBranchAddress("photon_chargeiso", &photonchargeiso);
  b1b.SetBranchAddress("photon_neutraliso", &photonneutraliso);
  b1b.SetBranchAddress("photon_photoniso", &photonphotoniso);
  b1b.SetBranchAddress("jet_pt", &jetpt);
  //b1a.SetBranchAddress("jet_pt_uncorrected", &jetptuncorrected);
  b1b.SetBranchAddress("jet_eta", &jeteta);
  b1b.SetBranchAddress("jet_phi", &jetphi);
  //b1a.SetBranchAddress("jet_eta_uncorrected", &jetetauncorrected);
  //b1a.SetBranchAddress("jet_phi_uncorrected", &jetphiuncorrected);
   // do similar for met

  //Int_t runnoprompt, eventnoprompt;
  //// variable declaration for tree2 ///////
  UInt_t          runnoprompt;
  Int_t           luminoprompt;
  ULong64_t       eventnoprompt;
  Int_t           verticesnprompt;
  Int_t           photonsnprompt;
  Float_t         metetprompt;
  Float_t         metXprompt;
  Float_t         metYprompt;
  Float_t         rhoprompt;
  vector<float>   *photoneprompt = 0;
  vector<float>   *photonetaprompt = 0;
  vector<float>   *photonphiprompt = 0;
  vector<float>   *photonshowershapeprompt = 0;
  vector<int>     *photonpixelseedprompt = 0;
  vector<float>   *photonchargeisoprompt = 0;
  vector<float>   *photonneutralisoprompt = 0;
  vector<float>   *photonphotonisoprompt = 0;
  vector<float>   *jetptprompt = 0;
  vector<float>   *jetptuncorrectedprompt = 0;
  vector<float>   *jetetaprompt = 0;
  vector<float>   *jetphiprompt = 0;
  vector<float>   *jetetauncorrectedprompt = 0;
  vector<float>   *jetphiuncorrectedprompt = 0;
  b2a.SetBranchAddress("runNo_prompt",&runnoprompt);
  b2a.SetBranchAddress("eventNo_prompt",&eventnoprompt);
  b2b.SetBranchAddress("photon_e_prompt", &photoneprompt);
  b2b.SetBranchAddress("lumiNo_prompt",&luminoprompt);
  b2b.SetBranchAddress("eventNo_prompt",&eventnoprompt);
  b2b.SetBranchAddress("vertices_n_prompt",&verticesnprompt);
  b2b.SetBranchAddress("photons_n_prompt",&photonsnprompt);
  b2b.SetBranchAddress("met_et_prompt",&metetprompt);
  b2b.SetBranchAddress("met_X_prompt",&metXprompt);
  b2b.SetBranchAddress("met_Y_prompt",&metYprompt);
  b2b.SetBranchAddress("rho_prompt", &rhoprompt);
  b2b.SetBranchAddress("photon_e_prompt", &photoneprompt);
  b2b.SetBranchAddress("photon_eta_prompt", &photonetaprompt);
  b2b.SetBranchAddress("photon_phi_prompt", &photonphiprompt);
  b2b.SetBranchAddress("photon_showershape_prompt", &photonshowershapeprompt);
  b2b.SetBranchAddress("photon_pixelseed_prompt", &photonpixelseedprompt);
  b2b.SetBranchAddress("photon_chargeiso_prompt", &photonchargeisoprompt);
  b2b.SetBranchAddress("photon_neutraliso_prompt", &photonneutralisoprompt);
  b2b.SetBranchAddress("photon_photoniso_prompt", &photonphotonisoprompt);
  b2b.SetBranchAddress("jet_pt_prompt", &jetptprompt);
  b2b.SetBranchAddress("jet_pt_uncorrected_prompt", &jetptuncorrectedprompt);
  b2b.SetBranchAddress("jet_eta_prompt", &jetetaprompt);
  b2b.SetBranchAddress("jet_phi_prompt", &jetphiprompt);
  b2b.SetBranchAddress("jet_eta_uncorrected_prompt", &jetetauncorrectedprompt);
  b2b.SetBranchAddress("jet_phi_uncorrected_prompt", &jetphiuncorrectedprompt);

  //// histograms ////
  TFile *FF = new TFile("PromptRecoElectronRawPfMet4.root","RECREATE");
  

  TH1F* h_met_rereco_diEMPtrhodiEMPtreweighted(new TH1F("h_met_rereco_diEMPtrhodiEMPtreweighted","rereco #slash{E}_{T}, diEMPt reweighted rho and diEMPt reweighted to candidate, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));

  TH1F *h_met_eeSample_DiEMPtRhoReweighted_ErrorFromDiEMPt[1000];
  TH1F *h_met_eeSample_DiEMPtRhoReweighted_ErrorFromRho[1000];

  
  TH1F* h_met_rereco_diEMPtreweighted(new TH1F("h_met_rereco_diEMPtreweighted","rereco #slash{E}_{T}, diEMPtreweighted to candidate, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
  TH1F* h_met_rereco_rhoreweighted(new TH1F("h_met_rereco_rhoreweighted","rereco #slash{E}_{T}, rhoreweighted to candidate, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
  TH1F* h_met_rereco_rhodiEMPtreweighted(new TH1F("h_met_rereco_rhodiEMPtreweighted","rereco #slash{E}_{T}, rhodiEMptreweighted to candidate, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.)); 
  TH1F* h_met_prompt(new TH1F("h_met_prompt","prompt all #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
  TH1F* h_met_rereco(new TH1F("h_met_rereco","rereco #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
  TH1F* h_met_difference(new TH1F("h_met_difference","rereco - prompt #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 100, -50., 50.));
  TH1F* h_met_percentdifference1(new TH1F("h_met_percentdifference1","(rereco - prompt)/prompt #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_met_percentdifference2(new TH1F("h_met_percentdifference2","(rereco - prompt)/rereco #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_metX_prompt(new TH1F("h_metX_prompt","prompt #slash{E}_{T} X component, Asymmetric Pt ;#slash{E}_{TX} (GeV);Events / GeV", 400, -200., 200.));
  TH1F* h_metX_rereco(new TH1F("h_metX_rereco","rereco #slash{E}_{T} X component, Asymmetric Pt ;#slash{E}_{TX} (GeV);Events / GeV", 400, -200., 200.));
  TH1F* h_metX_difference(new TH1F("h_metX_difference","rereco - prompt #slash{E}_{T} X component, Asymmetric Pt ;#slash{E}_{TX} (GeV);Events / GeV", 100, -50., 50.));
  TH1F* h_metX_percentdifference1(new TH1F("h_metX_percentdifference1","(rereco - prompt)/prompt #slash{E}_{T} X component, Asymmetric Pt ;#slash{E}_{TX} (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_metX_percentdifference2(new TH1F("h_metX_percentdifference2","(rereco - prompt)/rereco #slash{E}_{T} X component, Asymmetric Pt ;#slash{E}_{TX} (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_metY_prompt(new TH1F("h_metY_prompt","prompt #slash{E}_{T} Y component, Asymmetric Pt ;#slash{E}_{TY} (GeV);Events / GeV", 400, -200., 200.));
  TH1F* h_metY_rereco(new TH1F("h_metY_rereco","rereco #slash{E}_{T} Y component, Asymmetric Pt ;#slash{E}_{TY} (GeV);Events / GeV", 400, -200., 200.));
  TH1F* h_metY_difference(new TH1F("h_metY_difference","rereco - prompt #slash{E}_{T} Y component, Asymmetric Pt ;#slash{E}_{TY} (GeV);Events / GeV", 100, -50., 50.));
  TH1F* h_metY_percentdifference1(new TH1F("h_metY_percentdifference1","(rereco - prompt)/prompt #slash{E}_{T} Y component, Asymmetric Pt ;#slash{E}_{TY} (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_metY_percentdifference2(new TH1F("h_metY_percentdifference2","(rereco - prompt)/rereco #slash{E}_{T} Y component, Asymmetric Pt ;#slash{E}_{TY} (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_lead_photon_pt_prompt(new TH1F("h_lead_photon_pt_prompt","prompt lead photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 160, 0., 160.));
  TH1F* h_lead_photon_pt_rereco(new TH1F("h_lead_photon_pt_rereco","rereco lead photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 160, 0., 160.));
  TH1F* h_lead_photon_pt_difference(new TH1F("h_lead_photon_pt_difference","rereco - prompt lead photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 160, -20., 20.));
  TH1F* h_lead_photon_pt_percentdifference1(new TH1F("h_lead_photon_pt_percentdifference1","(rereco - prompt)/prompt lead photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_lead_photon_pt_percentdifference2(new TH1F("h_lead_photon_pt_percentdifference2","(rereco - prompt)/rereco lead photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_trail_photon_pt_prompt(new TH1F("h_trail_photon_pt_prompt","prompt trail photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 160, 0., 160.));
  TH1F* h_trail_photon_pt_rereco(new TH1F("h_trail_photon_pt_rereco","rereco trail photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 160, 0., 160.));
  TH1F* h_trail_photon_pt_difference(new TH1F("h_trail_photon_pt_difference","rereco - prompt trail photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 160, -20., 20.));
  TH1F* h_trail_photon_pt_percentdifference1(new TH1F("h_trail_photon_pt_percentdifference1","(rereco - prompt)/prompt trail photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_trail_photon_pt_percentdifference2(new TH1F("h_trail_photon_pt_percentdifference2","(rereco - prompt)/rereco trail photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_lead_photon_phi_prompt(new TH1F("h_lead_photon_phi_prompt","prompt lead photon #Phi, Asymmetric Pt ;#Phi ;Events", 80, -4., 4.));
  TH1F* h_lead_photon_phi_rereco(new TH1F("h_lead_photon_phi_rereco","rereco lead photon #Phi, Asymmetric Pt ;#Phi ;Events", 80, -4., 4.));
  TH1F* h_lead_photon_phi_difference(new TH1F("h_lead_photon_phi_difference","rereco - prompt lead photon #Phi, Asymmetric Pt ;#Phi ;Events", 240, -0.05, 0.05));
  TH1F* h_trail_photon_phi_prompt(new TH1F("h_trail_photon_phi_prompt","prompt trail photon #Phi, Asymmetric Pt ;#Phi ;Events ", 80, -4., 4.));
  TH1F* h_trail_photon_phi_rereco(new TH1F("h_trail_photon_phi_rereco","rereco trail photon #Phi, Asymmetric Pt ;#Phi ;Events ", 80, -4., 4.));
  TH1F* h_trail_photon_phi_difference(new TH1F("h_trail_photon_phi_difference","rereco - prompt trail photon #Phi, Asymmetric Pt ;#Phi ;Events ", 240, -0.05, 0.05));
  TH1F* h_lead_photon_eta_prompt(new TH1F("h_lead_photon_eta_prompt","prompt lead photon #eta, Asymmetric Pt ;#eta ;Events", 40, -2., 2.));
  TH1F* h_lead_photon_eta_rereco(new TH1F("h_lead_photon_eta_rereco","rereco lead photon #eta, Asymmetric Pt ;#eta ;Events", 40, -2., 2.));
  TH1F* h_lead_photon_eta_difference(new TH1F("h_lead_photon_eta_difference","rereco - prompt lead photon #eta, Asymmetric Pt ;#eta ;Events", 120, -0.05, 0.05));
  TH1F* h_trail_photon_eta_prompt(new TH1F("h_trail_photon_eta_prompt","prompt trail photon #eta, Asymmetric Pt ;#eta ;Events ", 40, -2., 2.));
  TH1F* h_trail_photon_eta_rereco(new TH1F("h_trail_photon_eta_rereco","rereco trail photon #eta, Asymmetric Pt ;#eta ;Events ", 40, -2., 2.));
  TH1F* h_trail_photon_eta_difference(new TH1F("h_trail_photon_eta_difference","rereco - prompt trail photon #eta, Asymmetric Pt ;#eta ;Events ", 120, -0.05, 0.05));
  TH1F* h_diEMPt_prompt(new TH1F("h_diEMPt_prompt","prompt diEMPt, Asymmetric Pt ;DiEM #P_T (GeV);Events / GeV", 250, 0., 250.));
  TH1F* h_diEMPt_rereco(new TH1F("h_diEMPt_rereco","rereco diEMPt, Asymmetric Pt ;DiEM #P_T (GeV);Events / GeV", 250, 0., 250.));
  TH1F* h_diEMPt_difference(new TH1F("h_diEMPt_difference","rereco - prompt diEMPt, Asymmetric Pt ;DiEM #P_T (GeV);Events / GeV", 240, -50., 50.));
  TH1F* h_diEMPt_percentdifference1(new TH1F("h_diEMPt_percentdifference1","(rereco - prompt)/prompt diEMPt, Asymmetric Pt ;DiEM #P_T (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_diEMPt_percentdifference2(new TH1F("h_diEMPt_percentdifference2","(rereco - prompt)/rereco diEMPt, Asymmetric Pt ;DiEM #P_T (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_reduced_met_prompt(new TH1F("h_reduced_met_prompt","prompt MET, Asymmetric Pt ;reduced MET (GeV);Events / GeV", 200, 0., 200.));
  TH1F* h_reduced_met_rereco(new TH1F("h_reduced_met_rereco","rereco MET, Asymmetric Pt ;reduced MET (GeV);Events / GeV", 200, 0., 200.));
  TH1F* h_reduced_met_difference(new TH1F("h_reduced_met_difference","rereco - prompt MET, Asymmetric Pt ;reduced MET (GeV);Events / GeV", 400, -100., 100.));
  TH1F* h_reduced_met_percentdifference1(new TH1F("h_reduced_met_percentdifference1","(rereco - prompt)/prompt MET, Asymmetric Pt ;reduced MET (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_reduced_met_percentdifference2(new TH1F("h_reduced_met_percentdifference2","(rereco - prompt)/rereco MET, Asymmetric Pt ;reduced MET (GeV);Events / GeV", 1000, -5., 5.));
  TH1F* h_jetpt_prompt(new TH1F("h_jetpt_prompt","prompt jet #P_T, Asymmetric Pt ; #P_T (GeV);Events / GeV", 250, 0., 250.));
  TH1F* h_jetpt_rereco(new TH1F("h_jetpt_rereco","rereco jet #P_T, Asymmetric Pt ; #P_T (GeV);Events / GeV", 250, 0., 250.));
  TH1F* h_jetpt_difference(new TH1F("h_jetpt_difference","rereco - prompt jet #P_T, Asymmetric Pt ; #P_T (GeV);Events / GeV", 400, -100., 100.));
  TH1F* h_jetphi_prompt(new TH1F("h_jetphi_prompt","prompt jet #Phi, Asymmetric Pt ; #Phi;Events ", 80, -4., 4.));
  TH1F* h_jetphi_rereco(new TH1F("h_jetphi_rereco","rereco jet #Phi, Asymmetric Pt ; #Phi;Events ", 80, -4., 4.));
  TH1F* h_jetphi_difference(new TH1F("h_jetphi_difference","rereco - prompt jet #Phi, Asymmetric Pt ; #Phi;Events ", 80, -1., 1.));
  TH1F* h_jeteta_prompt(new TH1F("h_jeteta_prompt","prompt jet #eta, Asymmetric Pt ; #eta;Events ", 40, -2., 2.));
  TH1F* h_jeteta_rereco(new TH1F("h_jeteta_rereco","rereco jet #eta, Asymmetric Pt ; #eta;Events ", 40, -2., 2.));
  TH1F* h_jeteta_difference(new TH1F("h_jeteta_difference","rereco - prompt jet #eta, Asymmetric Pt ; #eta;Events ", 40, -1., 1.));
  char *histNameee             = new char[50];
  char *histtitleee            = new char[50];
  char *histNameee2            = new char[50];
  char *histtitleee2           = new char[50];
  
   

   for (int d=0;d<1000; ++d) {
     sprintf(histNameee, "DiEMPtReweightedeeReReco%d",d+1);
     sprintf(histtitleee,"Reweighted ee MET with DIEMPT rereco %d",d+1);
     h_met_eeSample_DiEMPtRhoReweighted_ErrorFromDiEMPt[d]=new TH1F(histNameee,histtitleee,200,0.0,200.0);

     sprintf(histNameee2, "RhoReweightedeeReReco%d",d+1);
     sprintf(histtitleee2,"Reweighted ee MET with RHO rereco %d",d+1);
     h_met_eeSample_DiEMPtRhoReweighted_ErrorFromRho[d]=new TH1F(histNameee2,histtitleee2,200,0.0,200.0);

     
     
     
   }


  /// getting the reweighting from diEMPt and rho to candidate /////
  TFile *Fnew = new TFile("ElectronRatio.root","READ");
  TH1F *ElectronRho= (TH1F*)Fnew->Get("ElectronRho");
  TH1F *ElectronDiEMPt= (TH1F*)Fnew->Get("ElectronDiEMPt");
  
  TH1F *ElectrondiEMPtReweightedRho=(TH1F*)Fnew->Get("ElectrondiEMPtReweightedRho");



  



  TFile *F2 = new TFile("RandomRatioEE.root","READ");
  TH1F *EEDiEMPtRatio[1000]; //= (TH1F*)Fnew->Get("candidateeediEMPt");
  TH1F *EERhoRatio[1000];
  
  
   for(int g=0; g<1000; ++g){
     //printf("g:%d\n",g);
     char *nameGraph        = new char[50];
     char *nameGraph2       = new char[50];
    


     sprintf(nameGraph,"ElectronDiEMPtRandom%d",g+1);
     EEDiEMPtRatio[g] = (TH1F*)F2->Get(nameGraph);
     
     sprintf(nameGraph2,"ElectronRhoRandom%d",g+1);
     EERhoRatio[g] = (TH1F*)F2->Get(nameGraph2);
     
     
   }
 
  /////  getting the rereco event list ////////
  FILE *myfile;
  myfile=fopen("file_ele_ReReco4_Sorted.txt", "r");
  char line[50];
  int z1=0;
  const int linesize = 24256 ;
  ULong64_t ReRecoEventNo[linesize]={0};// rereco candidate: 308439, fake: 18010, electron: 2517353, first one 557950
  ULong64_t ReRecoEntry[linesize]={0};
  UInt_t ReRecoRunNo[linesize]={0};
  while(fgets(line, sizeof line, myfile) != NULL){
    fscanf(myfile,"%d\t%lld\t%lld", &ReRecoRunNo[z1], &ReRecoEventNo[z1], &ReRecoEntry[z1]);
    //printf("%d\t%lld\t%lld\n", ReRecoRunNo[z], ReRecoEventNo[z], ReRecoEntry[z]);
    z1++;
  }
  
  int   matched = 0;
  //Int_t nentries1 = (Int_t)b2a.GetEntries();
  //for (int i = 0; i<=nentries1; ++i) {
  //b2a.GetEntry(i);
  //b2b.GetEntry(i);
  
  //std::cout << "rereco events processed: " << i << std::endl;
  //b1b.GetEntry(i);
  //Int_t nentries2 = (Int_t)b2a.GetEntries();
  
  for(int f=0;f<linesize;++f){
    FILE *myfile1;
    myfile1=fopen("file_ele_Prompt4_Sorted.txt", "r");
    char line1[50];
    ULong64_t PromptEventNo=0;
    ULong64_t PromptEntry=0;
    UInt_t PromptRunNo=0;//387330
    int same=0;
    
    while(fgets(line1, sizeof line1, myfile1) != NULL){
      fscanf(myfile1,"%d\t%lld\t%lld", &PromptRunNo, &PromptEventNo, &PromptEntry);
    //printf("%d\t%lld\t%lld\n", ReRecoRunNo[z], ReRecoEventNo[z], ReRecoEntry[z]);
      
    
      if(ReRecoEventNo[f]==PromptEventNo && ReRecoRunNo[f]==PromptRunNo){same++;break;}
    }
  
    if(same==1){
      matched++;
      //if(matched>200000)break;
      if(matched%1000==0)std::cout << "matched events processed: " << matched << std::endl;
      b1a.GetEntry(ReRecoEntry[f]);
      b1b.GetEntry(ReRecoEntry[f]);
      b2a.GetEntry(PromptEntry);
      b2b.GetEntry(PromptEntry);
      
      float DIEMPT=findDiEMPt(photone->at(0),photoneta->at(0),photonphi->at(0),photone->at(1),photoneta->at(1),photonphi->at(1));
      float PxTotal = 0;
      float PyTotal = 0;
      for(size_t k = 0; k < photone->size(); ++k){
        float theta = 2*atan(exp(-photoneta->at(k)));
        float PX    = photone->at(k)*sin(theta)*cos(photonphi->at(k)); 
        float PY    = photone->at(k)*sin(theta)*sin(photonphi->at(k));
        PxTotal    += PX;
        PyTotal    += PY;
      }
      float redX    = metX+PxTotal;
      float redY    = metY+PyTotal;
      float redMET  = sqrt(redX*redX+redY*redY);
      h_lead_photon_pt_rereco->Fill(photone->at(0)*2*atan(exp(-photoneta->at(0))));
      h_trail_photon_pt_rereco->Fill(photone->at(1)*2*atan(exp(-photoneta->at(1))));
      h_lead_photon_phi_rereco->Fill(photonphi->at(0));
      h_trail_photon_phi_rereco->Fill(photonphi->at(1));
      h_lead_photon_eta_rereco->Fill(photoneta->at(0));
      h_trail_photon_eta_rereco->Fill(photoneta->at(1));
      h_diEMPt_rereco->Fill(DIEMPT);
      /* reweighting */
      int Bin1=ElectronDiEMPt->FindBin(DIEMPT);
      float Weight1=ElectronDiEMPt->GetBinContent(Bin1);
      int Bin2=ElectronRho->FindBin(rho);
      float Weight2=ElectronRho->GetBinContent(Bin2);
      int Bin22=ElectrondiEMPtReweightedRho->FindBin(rho);
      float Weight22=ElectrondiEMPtReweightedRho->GetBinContent(Bin22);
      h_met_rereco_diEMPtreweighted->Fill(metet,Weight1);
      h_met_rereco_rhoreweighted->Fill(metet,Weight2);
      h_met_rereco_rhodiEMPtreweighted->Fill(metet, Weight1*Weight2);

      h_met_rereco_diEMPtrhodiEMPtreweighted->Fill(metet, Weight22*Weight1);

      for( int s=0 ; s<1000; ++s){
        float Weightprop = EEDiEMPtRatio[s]->GetBinContent(Bin1);
        h_met_eeSample_DiEMPtRhoReweighted_ErrorFromDiEMPt[s]->Fill(metet, Weightprop);
      }
      for( int t=0 ; t<1000; ++t){
        float WeightFromRho = EERhoRatio[t]->GetBinContent(Bin22);
        h_met_eeSample_DiEMPtRhoReweighted_ErrorFromRho[t]->Fill(metet, WeightFromRho);
      }


      h_met_rereco->Fill(metet);
      h_metX_rereco->Fill(metX);
      h_metY_rereco->Fill(metY);
      h_reduced_met_rereco->Fill(redMET);
      for(size_t s=0; s<jetpt->size(); ++s){
        float dr=0;
        int count = 0;
        for(size_t j=0; j<photone->size(); ++j){
          dr=dRCalc(jeteta->at(s),jetphi->at(s),photoneta->at(j),photonphi->at(j));
          if(dr<0.6)count++;
        }
        if(count == 0){
          h_jetpt_rereco->Fill(jetpt->at(s));
          h_jeteta_rereco->Fill(jeteta->at(s));
          h_jetphi_rereco->Fill(jetphi->at(s));
        }
      } 
      ///// working with promptreco /////
      float DIEMPTprompt=findDiEMPt(photoneprompt->at(0),photonetaprompt->at(0),photonphiprompt->at(0),photoneprompt->at(1),photonetaprompt->at(1),photonphiprompt->at(1));
      float PxTotalprompt = 0;
      float PyTotalprompt = 0;
      for(size_t k = 0; k < photoneprompt->size(); ++k){
        float thetaprompt = 2*atan(exp(-photonetaprompt->at(k)));
        float PXprompt    = photoneprompt->at(k)*sin(thetaprompt)*cos(photonphiprompt->at(k)); 
        float PYprompt    = photoneprompt->at(k)*sin(thetaprompt)*sin(photonphiprompt->at(k));
        PxTotalprompt    += PXprompt;
        PyTotalprompt    += PYprompt;
      }
      float redXprompt    = metXprompt+PxTotalprompt;
      float redYprompt    = metYprompt+PyTotalprompt;
      float redMETprompt  = sqrt(redXprompt*redXprompt+redYprompt*redYprompt);
      h_lead_photon_pt_prompt->Fill(photoneprompt->at(0)*2*atan(exp(-photonetaprompt->at(0))));
      h_trail_photon_pt_prompt->Fill(photoneprompt->at(1)*2*atan(exp(-photonetaprompt->at(1))));
      h_lead_photon_phi_prompt->Fill(photonphiprompt->at(0));
      h_trail_photon_phi_prompt->Fill(photonphiprompt->at(1));
      h_lead_photon_eta_prompt->Fill(photonetaprompt->at(0));
      h_trail_photon_eta_prompt->Fill(photonetaprompt->at(1));
      h_diEMPt_prompt->Fill(DIEMPTprompt);
      h_met_prompt->Fill(metetprompt);
      h_metX_prompt->Fill(metXprompt);
      h_metY_prompt->Fill(metYprompt);
      h_reduced_met_prompt->Fill(redMETprompt);
      for(size_t p=0; p<jetptprompt->size(); ++p){
        float drprompt=0;
        int countprompt = 0;
        for(size_t j=0; j<photoneprompt->size(); ++j){
          drprompt=dRCalc(jetetaprompt->at(p),jetphiprompt->at(p),photonetaprompt->at(j),photonphiprompt->at(j));
          if(drprompt<0.6)countprompt++;
        }
        if(countprompt == 0){
          h_jetpt_prompt->Fill(jetptprompt->at(p));
          h_jeteta_prompt->Fill(jetetaprompt->at(p));
          h_jetphi_prompt->Fill(jetphiprompt->at(p));
        }
      }
        //// working with the differences ////
      h_lead_photon_pt_difference->Fill(photone->at(0)*2*atan(exp(-photoneta->at(0)))-photoneprompt->at(0)*2*atan(exp(-photonetaprompt->at(0))));
      h_lead_photon_pt_percentdifference1->Fill((photone->at(0)*2*atan(exp(-photoneta->at(0)))-photoneprompt->at(0)*2*atan(exp(-photonetaprompt->at(0))))/(photoneprompt->at(0)*2*atan(exp(-photonetaprompt->at(0)))));
      h_lead_photon_pt_percentdifference2->Fill((photone->at(0)*2*atan(exp(-photoneta->at(0)))-photoneprompt->at(0)*2*atan(exp(-photonetaprompt->at(0))))/(photone->at(0)*2*atan(exp(-photoneta->at(0)))));
      h_trail_photon_pt_percentdifference1->Fill((photone->at(1)*2*atan(exp(-photoneta->at(1)))-photoneprompt->at(1)*2*atan(exp(-photonetaprompt->at(1))))/(photoneprompt->at(1)*2*atan(exp(-photonetaprompt->at(1)))));
      h_trail_photon_pt_percentdifference2->Fill((photone->at(1)*2*atan(exp(-photoneta->at(1)))-photoneprompt->at(1)*2*atan(exp(-photonetaprompt->at(1))))/(photone->at(1)*2*atan(exp(-photoneta->at(1)))));
      h_trail_photon_pt_difference->Fill(photone->at(1)*2*atan(exp(-photoneta->at(1)))-photoneprompt->at(1)*2*atan(exp(-photonetaprompt->at(1))));
      h_lead_photon_phi_difference->Fill(photonphi->at(0)-photonphiprompt->at(0));
      h_trail_photon_phi_difference->Fill(photonphi->at(1)-photonphiprompt->at(1));
      h_lead_photon_eta_difference->Fill(photoneta->at(0)-photonetaprompt->at(0));
      h_trail_photon_eta_difference->Fill(photoneta->at(1)-photonetaprompt->at(1));
      h_diEMPt_difference->Fill(DIEMPT-DIEMPTprompt);
      h_diEMPt_percentdifference1->Fill((DIEMPT-DIEMPTprompt)/DIEMPTprompt);
      h_diEMPt_percentdifference2->Fill((DIEMPT-DIEMPTprompt)/DIEMPT);
      h_met_difference->Fill(metet-metetprompt);
      h_met_percentdifference1->Fill((metet-metetprompt)/metetprompt);
      h_met_percentdifference2->Fill((metet-metetprompt)/metet);
      h_metX_difference->Fill(metX-metXprompt);
      h_metX_percentdifference1->Fill((metX-metXprompt)/metXprompt);
      h_metX_percentdifference2->Fill((metX-metXprompt)/metX);
      h_metY_difference->Fill(metY-metYprompt);
      h_metY_percentdifference1->Fill((metY-metYprompt)/metYprompt);
      h_metY_percentdifference2->Fill((metY-metYprompt)/metY);
      h_reduced_met_difference->Fill(redMET-redMETprompt);
      h_reduced_met_percentdifference1->Fill((redMET-redMETprompt)/redMETprompt);
      h_reduced_met_percentdifference2->Fill((redMET-redMETprompt)/redMET);
      for(size_t s=0; s<jetpt->size(); ++s){
        float dr=0;
        int count = 0;
        for(size_t j=0; j<photone->size(); ++j){
          dr=dRCalc(jeteta->at(s),jetphi->at(s),photoneta->at(j),photonphi->at(j));
          if(dr<0.6)count++;
        }
        bool found = 0;
        for(size_t v=0; v<jetptprompt->size(); ++v){
          float drx=0;
          int countx = 0;
          for(size_t w=0; w<photoneprompt->size(); ++w){
            drx=dRCalc(jetetaprompt->at(v),jetphiprompt->at(v),photonetaprompt->at(w),photonphiprompt->at(w));
            if(drx<0.6)countx++;
          }
          if(count == 0 && countx == 0 && !found){
            h_jetpt_difference->Fill(jetpt->at(s)-jetptprompt->at(v));
            h_jeteta_difference->Fill(jeteta->at(s)-jetetaprompt->at(v));
            h_jetphi_difference->Fill(jetphi->at(s)-jetphiprompt->at(v));
            found =1;
          }
        }
      } 
    }
    fclose(myfile1);
  }
  fclose(myfile);
  FF->cd();
  FF->Write();

}
