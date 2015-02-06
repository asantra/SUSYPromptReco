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

void TwoTreeAll9(){

  ////// get trees  ///////
  TChain b1a("tree");
  TChain b1b("tree");
  TChain b2a("NNtree");
  TChain b2b("NNtree");
  //b1a.Add("PhotonAll_Tree.root");// contains corrected pfMet from rereco
  //b1b.Add("PhotonAll_Tree.root");
  b1a.Add("PhotonAll_PfMet_ReReco.root");// contains raw pfMet from rereco
  b1b.Add("PhotonAll_PfMet_ReReco.root");
  b2a.Add("PhotonAll_Prompt.root");
  b2b.Add("PhotonAll_Prompt.root");

  ///// variable declaration for tree 1 ///////
  UInt_t          runno;
  Int_t           lumino;
  ULong64_t       eventno;
  Int_t           verticesn;
  //Int_t           photonsn;
  Float_t         metet;
  Float_t         metX;
  Float_t         metY;
  Float_t         rho;
  /*vector<float>   *photone = 0;
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
  //vector<float>   *jetphiuncorrected = 0;*/
  b1a.SetBranchAddress("runNo",&runno);
  b1a.SetBranchAddress("eventNo",&eventno);
  //b1b.SetBranchAddress("photon_e", &photone);
  b1b.SetBranchAddress("lumiNo",&lumino);
  b1b.SetBranchAddress("eventNo",&eventno);
  b1b.SetBranchAddress("vertices_n",&verticesn);
  //b1b.SetBranchAddress("photons_n",&photonsn);
  b1b.SetBranchAddress("met_et",&metet);
  b1b.SetBranchAddress("met_X",&metX);
  b1b.SetBranchAddress("met_Y",&metY);
  b1b.SetBranchAddress("rho", &rho);
  /*b1b.SetBranchAddress("photon_eta", &photoneta);
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
  //b1a.SetBranchAddress("jet_phi_uncorrected", &jetphiuncorrected);*/
   // do similar for met

  //Int_t runnoprompt, eventnoprompt;
  //// variable declaration for tree2 ///////
  UInt_t          runnoprompt;
  Int_t           luminoprompt;
  ULong64_t       eventnoprompt;
  Int_t           verticesnprompt;
  //Int_t           photonsnprompt;
  Float_t         metetprompt;
  Float_t         metXprompt;
  Float_t         metYprompt;
  Float_t         rhoprompt;
  /*vector<float>   *photoneprompt = 0;
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
  vector<float>   *jetphiuncorrectedprompt = 0;*/
  b2a.SetBranchAddress("runNo_prompt",&runnoprompt);
  b2a.SetBranchAddress("eventNo_prompt",&eventnoprompt);
  //b2b.SetBranchAddress("photon_e_prompt", &photoneprompt);
  b2b.SetBranchAddress("lumiNo_prompt",&luminoprompt);
  b2b.SetBranchAddress("eventNo_prompt",&eventnoprompt);
  b2b.SetBranchAddress("vertices_n_prompt",&verticesnprompt);
  //b2b.SetBranchAddress("photons_n_prompt",&photonsnprompt);
  b2b.SetBranchAddress("met_et_prompt",&metetprompt);
  b2b.SetBranchAddress("met_X_prompt",&metXprompt);
  b2b.SetBranchAddress("met_Y_prompt",&metYprompt);
  b2b.SetBranchAddress("rho_prompt", &rhoprompt);
  /*b2b.SetBranchAddress("photon_e_prompt", &photoneprompt);
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
  b2b.SetBranchAddress("jet_phi_uncorrected_prompt", &jetphiuncorrectedprompt);*/

  //// histograms ////
  TFile *FF = new TFile("PromptRecoAllRawPfMet9.root","RECREATE");
   
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
  /*TH1F* h_lead_photon_pt_prompt(new TH1F("h_lead_photon_pt_prompt","prompt lead photon pt, Asymmetric Pt ;#P_T (GeV);Events / GeV", 160, 0., 160.));
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
  TH1F* h_jeteta_difference(new TH1F("h_jeteta_difference","rereco - prompt jet #eta, Asymmetric Pt ; #eta;Events ", 40, -1., 1.));*/

  FILE *myfile;
  myfile=fopen("file_All_ReReco9_Sorted.txt", "r");
  char line[50];
  const int linesize =  384842;
  int z1=0;
  ULong64_t ReRecoEventNo[linesize]={0};// rereco candidate: 308439, fake: 18010, electron: 2517353, first one 557950
  ULong64_t ReRecoEntry[linesize]={0};
  UInt_t ReRecoRunNo[linesize]={0};
  while(fgets(line, sizeof line, myfile) != NULL){
    fscanf(myfile,"%d\t%lld\t%lld", &ReRecoRunNo[z1], &ReRecoEventNo[z1], &ReRecoEntry[z1]);
    //printf("%d\t%lld\t%lld\n", ReRecoRunNo[z], ReRecoEventNo[z], ReRecoEntry[z]);
    z1++;
  }
  //std:: cout << "let's see" << std::endl;
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
    myfile1=fopen("file_All_Prompt9_Sorted.txt", "r");
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
      if(matched%100==0)std::cout << "matched events processed: " << matched << std::endl;
      b1a.GetEntry(ReRecoEntry[f]);
      b1b.GetEntry(ReRecoEntry[f]);
      b2a.GetEntry(PromptEntry);
      b2b.GetEntry(PromptEntry);
      //b2a.GetEntry(i);
      //printf("eventNoprompt: %llu , runNoprompt: %u, PromptEventNo: %llu, PromptRunNo: %u, eventNo: %llu, runNo: %u, ReRecoEvent: %llu, ReRecoRun:%u \n",eventnoprompt,runnoprompt,PromptEventNo, PromptRunNo,eventno, runno, ReRecoEventNo[f], ReRecoRunNo[f]);
    //std::cout << "prompt events processed: " << j << std::endl;
    
    //b2b.GetEntry(j);
      
        
      //b1b.GetEntry(i);
        //b2b.GetEntry(j);
        /* /////// Photon, electron and fake Selection /////// */
//         if(matched == 98){
//           std::cout << "runno1: " << runno << " and runnoprompt1: " << runnoprompt << std::endl;
//           std::cout << "eventno1: " << eventno << " and eventnoprompt1: " << eventnoprompt << std::endl;
//           std::cout << "met: " << metet << " and metX: " << metX << std::endl;
//           std::cout << "metprompt: " << metetprompt << " and metXprompt: " << metXprompt << std::endl;
//           printf("%d\t%lld\t%lld\n", ReRecoRunNo[f], ReRecoEventNo[f], ReRecoEntry[f]);
//         }
      
        /*if(matched==98){
          std::cout << "rho: " << rho << " and rhoprompt: " << rhoprompt << std::endl;
          //std::cout << "photone: " << photone->at(0) << " and photoneprompt: " << photoneprompt->at(0) << std::endl;
          //std::cout << "photonphi: " << photonphi->at(0) << " and photonphiprompt: " << photonphiprompt->at(0) << std::endl;
        }*/
        ///// histograms for rereco ///////
        /*if(matched==98 or matched==113 or matched==239 or matched==255 or matched==273 or matched==287 or matched==309)continue;
        if(matched==338 or matched==434 or matched==449 or matched==454 or matched==524 or matched==597 or matched==614)continue;
        if(matched==620 or matched==735 or matched==776 or matched==781 or matched==784 or matched==821 or matched==846)continue;
        if(matched==854 or matched==892 or matched==968 or matched==1002 or matched==1039 or matched==1056 or matched==1064)continue;
        if(matched==1093 or matched==1117 or matched==1269 or matched==1301 or matched==1310 or matched==1317 or matched==1364)continue;
        if(matched==1385 or matched==1407 or matched==1566 or matched==1648 or matched==1649 or matched==1684 or matched==1696)continue;
        if(matched==1708 or matched==1769 or matched==1783 or matched==1795 or matched==1811 or matched==1837 or matched==1842)continue;
        if(matched==1865 or matched==1874 or matched==1877 or matched==1899 or matched==1977 or matched==2021 or matched==2065 or matched==2122)continue;*/
        //if(matched==2148 or matched==2163 or matched==2185 or matched==2276 or matched==2282 or matched==2283 or matched==2295)continue;
      /*float DIEMPT=findDiEMPt(photone->at(0),photoneta->at(0),photonphi->at(0),photone->at(1),photoneta->at(1),photonphi->at(1));
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
      h_diEMPt_rereco->Fill(DIEMPT);*/
      h_met_rereco->Fill(metet);
      h_metX_rereco->Fill(metX);
      h_metY_rereco->Fill(metY);
      /*h_reduced_met_rereco->Fill(redMET);
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
      h_diEMPt_prompt->Fill(DIEMPTprompt);*/
      h_met_prompt->Fill(metetprompt);
      h_metX_prompt->Fill(metXprompt);
      h_metY_prompt->Fill(metYprompt);
      /*h_reduced_met_prompt->Fill(redMETprompt);
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
      h_trail_photon_pt_difference->Fill(photone->at(1)*2*atan(exp(-photoneta->at(1)))-photoneprompt->at(1)*2*atan(exp(-photonetaprompt->at(1))));
      h_lead_photon_phi_difference->Fill(photonphi->at(0)-photonphiprompt->at(0));
      h_trail_photon_phi_difference->Fill(photonphi->at(1)-photonphiprompt->at(1));
      h_lead_photon_eta_difference->Fill(photoneta->at(0)-photonetaprompt->at(0));
      h_trail_photon_eta_difference->Fill(photoneta->at(1)-photonetaprompt->at(1));
      h_diEMPt_difference->Fill(DIEMPT-DIEMPTprompt);
      h_diEMPt_percentdifference1->Fill((DIEMPT-DIEMPTprompt)/DIEMPTprompt);
      h_diEMPt_percentdifference2->Fill((DIEMPT-DIEMPTprompt)/DIEMPT);*/
      h_met_difference->Fill(metet-metetprompt);
      h_met_percentdifference1->Fill((metet-metetprompt)/metetprompt);
      h_met_percentdifference2->Fill((metet-metetprompt)/metet);
      h_metX_difference->Fill(metX-metXprompt);
      h_metX_percentdifference1->Fill((metX-metXprompt)/metXprompt);
      h_metX_percentdifference2->Fill((metX-metXprompt)/metX);
      h_metY_difference->Fill(metY-metYprompt);
      h_metY_percentdifference1->Fill((metY-metYprompt)/metYprompt);
      h_metY_percentdifference2->Fill((metY-metYprompt)/metY);
      /*h_reduced_met_difference->Fill(redMET-redMETprompt);
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
      }*/ 
    }
    fclose(myfile1);
  }
  fclose(myfile);
  FF->cd();
  FF->Write();

}
