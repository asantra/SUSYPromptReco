{
  TTree *tree1, *tree2; //pointers to your 3 Trees
  
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PhotonAll_Prompt.root");
  if (!f || !f->IsOpen()) {
    f = new TFile("PhotonAll_Prompt.root");
  }
  f->GetObject("NNtree",tree1);
  //Init(tree1);
  TFile *f2 = (TFile*)gROOT->GetListOfFiles()->FindObject("PhotonAll_Tree.root");
  if (!f2 || !f2->IsOpen()) {
     f2 = new TFile("PhotonAll_Tree.root");
  }
  f2->GetObject("tree",tree2);
  //Init(tree2);
  TFile *F2 = new TFile("ReRecoPrompt.root","RECREATE");
  TList *list = new TList;
  list->Add(tree1);
  list->Add(tree2);
  TTree *newtree = TTree::MergeTrees(list);
  newtree->SetName("newtree");
  F2->cd();
  F2->Write();

}
