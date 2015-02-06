// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.cc
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.cc,v 1.15 2012/08/31 11:33:53 bfrancis Exp $
//

#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>

#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>

#include "SusyEventAnalyzer.h"
#include "SusyEventPrinter.h"

//#include "JetCorrectorParameters.h"
//#include "FactorizedJetCorrector.h"

#include "../jec/JetMETObjects/interface/JetCorrectorParameters.h"
#include "../jec/JetMETObjects/interface/FactorizedJetCorrector.h"


template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}


void SusyEventAnalyzer::InitializePerEvent() {

}


bool SusyEventAnalyzer::isSameObject(TLorentzVector& p1, TLorentzVector& p2) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < 0.5) return true;
  return false;
}


void
photonEffectiveAreas(double _eta, double* _effA)
{
  double& effACH(_effA[0]);
  double& effANH(_effA[1]);
  double& effAPh(_effA[2]);

  // Source: CutBasedPhotonID2012 twiki
  if(_eta < 1.){
    effACH = 0.012;
    effANH = 0.03;
    effAPh = 0.148;
  }
  else if(_eta < 1.479){
    effACH = 0.010;
    effANH = 0.057;
    effAPh = 0.13;
  }
}

float SusyEventAnalyzer::d0correction(TVector3& beamSpot, susy::Track& track) const {

  float d0 = track.d0() - beamSpot.X()*std::sin(track.phi()) + beamSpot.Y()*std::cos(track.phi());
  return d0;
}


bool SusyEventAnalyzer::PassTrigger(TString path) {
  bool pass = false;
  for(susy::TriggerMap::iterator it = event->hltMap.begin(); it != event->hltMap.end(); it++) {
    if(it->first.Contains(path) && (int(it->second.second)) ) {
      pass = true;
      break;
    }
  }
  return pass;
}


bool SusyEventAnalyzer::PassTriggers() {
  bool pass = false;
  for(std::vector<TString>::iterator it = hltNames.begin(); it != hltNames.end(); it++) {
    if(PassTrigger(*it)) {
      pass = true;
      break;
    }
  }
  return pass;
}



void SusyEventAnalyzer::Loop() {

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  std::cout << "total events in files  : " << nentries << std::endl;

  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

  std::cout << "events to be processed : " << processNEvents << std::endl; 


  if(printLevel > 0) std::cout << "Initialize event counters." << std::endl;
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;

  int nFiltered = 0;
  TTree* filterTree = 0;

  if(enableFilter) {
    TFile* filterFile = new TFile(filtered_file_name,"RECREATE");
    filterTree = (TTree*) fChain->GetTree()->CloneTree(0);
    filterTree->SetAutoSave();
  }


  // open hist file and define histograms

  TFile* fout = new TFile("hist_"+ds+".root","RECREATE");

  fout->cd();

  /*TH1F* h_vtxZ = new TH1F("vtxZ","Z position of the primary vertex;Z (cm);Events",100,-50.0,50.0);
  TH1F* h_bsZ = new TH1F("bsZ","Z position of the beam spot;Z (cm);Events",100,-50.0,50.0);
  TH1F* h_met = new TH1F("met","missing transverse energy;#slash{E}_{T} (GeV);Events",200,0.0,1000.0);
  TH1F* h_sumEt = new TH1F("sumEt","Scalar sum of all calorimeter energy;#sigmaE_{T} (GeV);Events",200,0.0,2000.0);*/

  /* Ntuple vectors */
    
  photon_e_prompt            = new std::vector<float>();
  photon_eta_prompt          = new std::vector<float>();
  photon_phi_prompt          = new std::vector<float>();
  photon_showershape_prompt  = new std::vector<float>();
  photon_pixelseed_prompt    = new std::vector<int>(); 
  photon_chargeiso_prompt    = new std::vector<float>();
  photon_neutraliso_prompt   = new std::vector<float>();
  photon_photoniso_prompt    = new std::vector<float>();
  //photon_worstiso     = new std::vector<float>();
  jet_pt_prompt              = new std::vector<float>();
  jet_pt_uncorrected_prompt  = new std::vector<float>();
  jet_eta_prompt             = new std::vector<float>();
  jet_phi_prompt             = new std::vector<float>();
  jet_eta_uncorrected_prompt = new std::vector<float>();
  jet_phi_uncorrected_prompt = new std::vector<float>();


  TTree* NNtree = new TTree("NNtree", "GMSB SUSY PromptReco tree");

  NNtree->Branch("runNo_prompt",   & runNo_prompt,   "runNo_prompt/I");
  NNtree->Branch("lumiNo_prompt",  & lumiNo_prompt,  "lumiNo_prompt/I");
  NNtree->Branch("eventNo_prompt", & eventNo_prompt, "eventNo_prompt/l");

  NNtree->Branch("vertices_n_prompt",  & vertices_n_prompt,  "vertices_n_prompt/I");
  NNtree->Branch("photons_n_prompt",   & photons_n_prompt,   "photons_n_prompt/I");

  NNtree->Branch("met_et_prompt", & met_et_prompt, "met_et_prompt/F");
  NNtree->Branch("met_X_prompt",  & met_X_prompt,  "met_X_prompt/F");
  NNtree->Branch("met_Y_prompt",  & met_Y_prompt,  "met_Y_prompt/F");
  NNtree->Branch("rho_prompt",    & rho_prompt,    "rho_prompt/F");

    
  NNtree->Branch("photon_e_prompt",           "vector<float>", & photon_e_prompt);
  NNtree->Branch("photon_eta_prompt",         "vector<float>", & photon_eta_prompt);
  NNtree->Branch("photon_phi_prompt",         "vector<float>", & photon_phi_prompt);
  NNtree->Branch("photon_showershape_prompt", "vector<float>", & photon_showershape_prompt);
  NNtree->Branch("photon_pixelseed_prompt",   "vector<int>",   & photon_pixelseed_prompt);
  NNtree->Branch("photon_chargeiso_prompt",   "vector<float>", & photon_chargeiso_prompt);
  NNtree->Branch("photon_neutraliso_prompt",  "vector<float>", & photon_neutraliso_prompt);
  NNtree->Branch("photon_photoniso_prompt",   "vector<float>", & photon_photoniso_prompt);
  //Ntree->Branch("photon_worstiso",    "vector<float>", & photon_worstiso);
  NNtree->Branch("jet_pt_prompt",             "vector<float>", & jet_pt_prompt);
  NNtree->Branch("jet_pt_uncorrected_prompt", "vector<float>", & jet_pt_uncorrected_prompt);
  NNtree->Branch("jet_eta_prompt",            "vector<float>", & jet_eta_prompt);
  NNtree->Branch("jet_phi_prompt",            "vector<float>", & jet_phi_prompt);
  NNtree->Branch("jet_eta_uncorrected_prompt","vector<float>", & jet_eta_uncorrected_prompt);
  NNtree->Branch("jet_phi_uncorrected_prompt","vector<float>", & jet_phi_uncorrected_prompt);
    


  // to check duplicated events
  std::map<int, std::set<int> > allEvents;

  // start event looping

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {

    if(printLevel > 0) std::cout << "Get the tree contents." << std::endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      std::cout << int(jentry) << " events processed with run="
		<< event->runNumber << ", event=" << event->eventNumber << std::endl;
    }

    ////////// Clear Flat Ntuple ///////////
      
    photon_e_prompt->clear();     photon_showershape_prompt->clear();       photon_pixelseed_prompt->clear();   photon_eta_prompt->clear();
    photon_phi_prompt->clear();   photon_chargeiso_prompt->clear();         photon_neutraliso_prompt->clear();
    photon_photoniso_prompt->clear();  /*photon_worstiso->clear();*/ jet_pt_prompt->clear();             jet_eta_prompt->clear();
    jet_phi_prompt->clear();      jet_pt_uncorrected_prompt->clear();       jet_phi_uncorrected_prompt->clear();jet_eta_uncorrected_prompt->clear();


    if(printLevel > 0) std::cout << "Initialize any global variables to be reset per event." << std::endl;

    InitializePerEvent();


    if(printLevel > 0) std::cout << "Apply good run list." << std::endl;

    
    // uncomment this to use the Json file to flag good data (or bad depending on your outlook)    
    // if(!isInJson(event->runNumber,event->luminosityBlockNumber)) continue;

    // uncomment this to print all ntuple variables
    //Print(*event);

    if(printLevel > 0) std::cout << "Check duplicated events for data only." << std::endl;

    bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
    if(event->isRealData && duplicateEvent) continue;
 

    // remove events filtered by optional met filters
    if(!event->passMetFilters()) continue;

    if(printLevel > 0) std::cout << "Setup object vectors." << std::endl;

    // classify photon objects

    // loose objects have all standard cuts except for isolation
    std::vector<susy::Photon*>   loose_photons;
    std::vector<susy::Photon*>   allphotons;

    // tight objects hava isolation cuts applied on top of loose objects
    std::vector<susy::Photon*>   tight_photons;

    // same as tight except for nPixelSeeds > 0
    std::vector<susy::Photon*>   ele_photons;

    // same as tight except for reversing either trackIso or sigmaIetaIeta
    std::vector<susy::Photon*>   fake_photons;

    std::vector<susy::CaloJet*>  caloJets;
    std::vector<susy::PFJet*>    pfJets;

    if(printLevel > 0) std::cout << "Find primary vertex in the event." << std::endl;

    TVector3* primVtx = 0;
    if(event->vertices.size() > 0) primVtx = &(event->vertices[0].position);

    //if(primVtx) h_vtxZ->Fill(primVtx->Z());
    //h_bsZ->Fill(event->beamSpot.Z());


    if(printLevel > 0) std::cout << "Find loose and tight photons in the event." << std::endl;

    std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");

    if(phoMap != event->photons.end()) {

      for(std::vector<susy::Photon>::iterator it = phoMap->second.begin();
	  it != phoMap->second.end(); it++) {

	// fiducial cuts. Look for only barrel now
	if(!it->isEB()) continue;

	// Et cuts, 25 GeV for trailing photons. Will apply tighter for the leading one.
	if(it->momentum.Et() < 25.0) continue;

        // optional Spike cleaning
        if(it->r9 > 1.0) continue;

        // H/E (in trigger, 0.15 for EB, 0.10 for EE)
        bool heCut = (it->hadronicOverEm < 0.05);
        
        // sigma_ietaieta (in trigger 0.014 for EB, 0.034 for EE)
        bool sIetaCut = (it->sigmaIetaIeta < 0.013);

        double effA[3];
        photonEffectiveAreas(fabs(it->momentum.Eta()), effA);

        // Ecal Isolation
        bool ecalIsoCut = (it->ecalRecHitSumEtConeDR04 < 4.2 + 0.006 * it->momentum.Et());

        // Hcal Isolation
        bool hcalIsoCut = (it->hcalTowerSumEtConeDR04() < 2.2 + 0.0025 * it->momentum.Et());

        // Track Isolation
        bool trackIsoCut = (it->trkSumPtHollowConeDR04 < 2.0 + 0.001 * it->momentum.Et());

        bool pixelCut = (it->nPixelSeeds == 0);

        // loose & tight ID variables
        bool looseCut = heCut && ecalIsoCut && hcalIsoCut;
        bool tightCut = looseCut && pixelCut && sIetaCut && trackIsoCut;
        bool eleClass  = looseCut && !pixelCut && sIetaCut && trackIsoCut;
        bool fakeClass = looseCut && pixelCut && !(sIetaCut && trackIsoCut);

        allphotons.push_back(&*it);

        if(looseCut) {
          loose_photons.push_back(&*it);
        }
        if(tightCut) {
          tight_photons.push_back(&*it);
        }
        if(eleClass) {
          ele_photons.push_back(&*it);
        }
        if(fakeClass) {
          fake_photons.push_back(&*it);
        }

      }// for photon
    }// else

    // sort photons by Et
    std::sort(loose_photons.begin(),loose_photons.end(),EtGreater<susy::Photon>);
    std::sort(allphotons.begin(),allphotons.end(),EtGreater<susy::Photon>);
    std::sort(tight_photons.begin(),tight_photons.end(),EtGreater<susy::Photon>);
    std::sort(ele_photons.begin(),ele_photons.end(),EtGreater<susy::Photon>);
    std::sort(fake_photons.begin(),fake_photons.end(),EtGreater<susy::Photon>);

    if(allphotons.size() != 0 && allphotons[0]->momentum.Pt() < 40.) allphotons.clear();

    /////// MET //////////
    std::map<TString, susy::MET>::iterator met_it = event->metMap.find("pfMet");
    if(met_it == event->metMap.end()) {
      std::cout << "MET map is not available!!!" << std::endl;
      continue;
    }
    susy::MET* met = &(met_it->second);

  
    ///////// Working with electrons, fakes and candidate samples ///////////
      
      
    runNo_prompt        = event->runNumber;
    lumiNo_prompt       = event->luminosityBlockNumber;
    eventNo_prompt      = event->eventNumber;
    vertices_n_prompt   = event->vertices.size();
    rho_prompt          = event->rho25;
    met_et_prompt       = met->met();
    met_X_prompt        = met->metX();
    met_Y_prompt        = met->metY();
    
    if(allphotons.size()>=2){
      photons_n_prompt = allphotons.size();
      for (size_t i = 0; i < allphotons.size(); i++) {
	photon_e_prompt->push_back(allphotons[i]->momentum.E());
        photon_eta_prompt->push_back(allphotons[i]->momentum.Eta());
	photon_phi_prompt->push_back(allphotons[i]->momentum.Phi());
        //cout << "photon phi: " << allPhotons[i]->momentum.Phi() << endl;
        photon_showershape_prompt->push_back(allphotons[i]->sigmaIetaIeta);
        photon_chargeiso_prompt->push_back(allphotons[i]->chargedHadronIso);
        photon_neutraliso_prompt->push_back(allphotons[i]->neutralHadronIso);
        photon_photoniso_prompt->push_back(allphotons[i]->photonIso);
        //photon_worstiso->push_back(allphotons[i]->worstOtherVtxChargedHadronIso);
        //cout << "photon worst iso: " << allPhotons[i]->worstOtherVtxChargedHadronIso << "for entry:" << iEntry << endl;
        photon_pixelseed_prompt->push_back(allphotons[i]->nPixelSeeds);
      }
    }
      

    if(printLevel > 0) std::cout << "Find caloJets in the event." << std::endl;
      
    std::map<TString,susy::CaloJetCollection>::iterator caloJets_it = event->caloJets.find("ak5");

    if(caloJets_it != event->caloJets.end()){

      susy::CaloJetCollection& jetColl = caloJets_it->second;

      for(std::vector<susy::CaloJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {

	std::map<TString,Float_t>::iterator s_it = it->jecScaleFactors.find("L2L3");
	if (s_it == it->jecScaleFactors.end()) {
	  std::cout << "JEC is not available for this jet!!!" << std::endl;
	  continue;
	}
	float scale = s_it->second;

        if(printLevel > 2) std::cout << "CaloJet stored (" << scale << ")" << std::endl;

	TLorentzVector corrP4 = scale * it->momentum;

	if(std::abs(corrP4.Eta()) > 3.0) continue;

	bool same = false;

	for(std::vector<susy::Photon*>::iterator m_it = tight_photons.begin();
	    m_it != tight_photons.end(); m_it++){
	  if(isSameObject(corrP4,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	//	if(pt < 20) continue;

	caloJets.push_back(&*it);

      }// for jet
    }// else

    std::sort(caloJets.begin(),caloJets.end(),EtGreater<susy::CaloJet>);


    if(printLevel > 0) std::cout << "Find pfJets in the event." << std::endl;

    std::vector<susy::PFJet const*> allJets;
    std::vector<susy::PFJet const*> allJets_uncorrected;
      
    std::map<TString,susy::PFJetCollection>::iterator pfJets_it = event->pfJets.find("ak5");
    if(pfJets_it == event->pfJets.end()){
      if(event->pfJets.size() > 0) std::cout << "JetCollection is not available!!!" << std::endl;
    }
    else {

      susy::PFJetCollection& jetColl = pfJets_it->second;

      for(std::vector<susy::PFJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {

	std::map<TString,Float_t>::iterator s_it = it->jecScaleFactors.find("L2L3");
	if (s_it == it->jecScaleFactors.end()) {
	  std::cout << "JEC is not available for this jet!!!" << std::endl;
	  continue;
	}
	float scale = s_it->second;

        if(printLevel > 2) std::cout << "PFJet stored (" << scale << ")" << std::endl;

        if(it->momentum.Pt() > 30 && it->momentum.Eta()<=3.0)allJets_uncorrected.push_back(&*it);

	TLorentzVector corrP4 = scale * it->momentum;

	if(std::abs(corrP4.Eta()) > 3.0) continue;

	bool same = false;

	for(std::vector<susy::Photon*>::iterator m_it = tight_photons.begin();
	    m_it != tight_photons.end(); m_it++){
	  if(isSameObject(corrP4,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
        pfJets.push_back(&*it);
	if(same) continue;

	//	if(pt < 20) continue;

	
        if(it->momentum.Pt() > 30)allJets.push_back(&*it);

      }// for jet
    }// else

    std::sort(pfJets.begin(),pfJets.end(),EtGreater<susy::PFJet>);
    std::sort(allJets.begin(),allJets.end(),EtGreater<susy::PFJet>);
    std::sort(allJets_uncorrected.begin(),allJets_uncorrected.end(),EtGreater<susy::PFJet>);

    if(allJets_uncorrected.size() > 0){
        for(unsigned int j = 0; j < allJets_uncorrected.size(); ++j){
          jet_pt_uncorrected_prompt->push_back(allJets_uncorrected[j]->momentum.Pt());
          jet_eta_uncorrected_prompt->push_back(allJets_uncorrected[j]->momentum.Eta());
          jet_phi_uncorrected_prompt->push_back(allJets_uncorrected[j]->momentum.Phi());
      }
    }
    if(allJets.size() > 0){
      for(unsigned int j = 0; j < allJets.size(); ++j){
        jet_pt_prompt  ->push_back(allJets[j]->momentum.Pt());
        //cout << "jet_pt"<<allJets[j]->momentum.Pt()<<"for entry:"<<iEntry << endl;
        jet_eta_prompt ->push_back(allJets[j]->momentum.Eta());
        //cout << "jet_eta"<<allJets[j]->momentum.Eta()<<"for entry:"<<iEntry << endl;
        jet_phi_prompt ->push_back(allJets[j]->momentum.Phi());
      }
    }


    if(printLevel > 0) std::cout << "Apply trigger selection in the event." << std::endl;

    bool passHLT = (useTrigger ? PassTriggers() : true);

    if(printLevel > 0) std::cout << "Select which met will be used in the event." << std::endl;

    

    if(printLevel > 0) {
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "              event summary" << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "loose_photons     : " << loose_photons.size() << std::endl;
      std::cout << "tight_photons     : " << tight_photons.size() << std::endl;
      std::cout << "ele_photons       : " << ele_photons.size() << std::endl;
      std::cout << "fake_photons      : " << fake_photons.size() << std::endl;
      std::cout << "caloJets          : " << caloJets.size() << std::endl;
      std::cout << "pfJets            : " << pfJets.size() << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "met               : " << met->met() << std::endl;
    } 


    if(printLevel > 0) std::cout << "Apply event level cuts from now on..." << std::endl;


    // filter conditions

    if(enableFilter) {
      bool filterThis = (loose_photons.size() > 0);
      if(filterThis) {
	nFiltered++;
	filterTree->Fill();
      }
    }// if(enableFilter)


    // event counter

    nCnt[0]++; // total number of events

    if(!passHLT) continue;

    nCnt[1]++;
 
    NNtree->Fill();

    if(loose_photons.size() == 0) continue;

    nCnt[2]++;

    //h_met->Fill(met->met());
    //h_sumEt->Fill(met->sumEt);


    // two photons
    if(tight_photons.size() >= 2) {
      nCnt[3]++;
    }

    // one photon + one electron
    if(tight_photons.size() >= 1 && ele_photons.size() >= 1) {
      nCnt[4]++;
    }

    // two electrons
    if(ele_photons.size() >= 2) {
      nCnt[5]++;
    }

    // one photon + one fake
    if(tight_photons.size() >= 1 && fake_photons.size() >= 1) {
      nCnt[6]++;
    }

    // two fakes
    if(fake_photons.size() >= 2) {
      nCnt[7]++;
    }


    if(met->met() < 50.0){ 

    nCnt[8]++;
    }
  
    

  } // for jentry

  
  // end of event loop and print summary

  std::cout << " ----------------- Job Summary ----------------- " << std::endl;
  std::cout << " Total events            : " << nCnt[0] << std::endl;
  std::cout << " HLT passed              : " << nCnt[1] << " (" << nCnt[1]/float(nCnt[0]) << ") wrt total events" << std::endl;
  std::cout << " loose_photons > 0       : " << nCnt[2] << " (" << nCnt[2]/float(nCnt[1]) << ") wrt HLT" << std::endl;
  std::cout << " gg events               : " << nCnt[3] << " (" << nCnt[3]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ge events               : " << nCnt[4] << " (" << nCnt[4]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ee events               : " << nCnt[5] << " (" << nCnt[5]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " gf events               : " << nCnt[6] << " (" << nCnt[6]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ff events               : " << nCnt[7] << " (" << nCnt[7]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " met > 50 GeV            : " << nCnt[8] << " (" << nCnt[8]/float(nCnt[1]) << ")" << std::endl;

  if(enableFilter) {
    std::cout << " --------------- Filtered events --------------- " << std::endl;
    std::cout << " filtered events         : " << nFiltered << " (" << nFiltered/float(nCnt[0]) << ")" << std::endl;
  }
  std::cout << " ----------------------------------------------- " << std::endl;

  // close the output file

  fout->cd();
  fout->Write();
  fout->Close();

  if(enableFilter) {
    filterTree->GetCurrentFile()->cd();
    filterTree->GetCurrentFile()->Write();
    filterTree->GetCurrentFile()->Close();
  }

}

