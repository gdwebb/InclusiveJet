class StJetCandidate;
class StJetVertex;
class StChain *chain;
#include <map>
#include <TObject.h>

void RunInclusiveJetCode(int nentries = 1e3,
			 const char* jetfile  = "/gpfs01/star/i_bnl/gdwebb/Run9/500GeVJets/Data/10103041/st_physics_*.jets.root",
			 const char* skimfile = "/gpfs01/star/i_bnl/gdwebb/Run9/500GeVJets/Data/10103041/st_physics_*.skim.root",
			 const char* uefile = "/gpfs01/star/i_bnl/gdwebb/Run9/500GeVJets/Data/10103041/st_physics_*.ue.root",
			 const char* outfile = "test500GeVjets10103041data.root",
			 int isEmbd = 0,
			 int isFzd = 0)
{
  cout << "nentries = " << nentries << endl;
  cout << "jetfile  = " << jetfile  << endl;
  cout << "skimfile = " << skimfile << endl;
  cout << "uefile = " << uefile << endl;
  cout << "outfile = " << outfile << endl;

  // Load libraries
  gROOT->Macro("loadMuDst.C");
  gSystem->Load("StJetEvent"); 
  gSystem->Load("StJetSkimEvent");
  gSystem->Load("StInclusiveJetMaker");

  chain = new StChain("StChain");
  // chain->SetDebug(1);
  // Open jet & skim files
  TChain* jetChain = new TChain("jet");
  TChain* skimChain = new TChain("jetSkimTree");
  TChain* ueChain = new TChain("ue");

  jetChain->Add(jetfile);
  skimChain->Add(skimfile);
  ueChain->Add(uefile);
  StInclusiveJetMaker* inclusivejet = new StInclusiveJetMaker("StInclusiveJetMaker","AntiKtR060NHits12",jetChain,skimChain,ueChain,isEmbd,isFzd,outfile);

  chain->Init();
  chain->PrintInfo();
  chain->ls(3);
  // Event loop
  Int_t ntotal = 0 ;
  for (int iEntry = 0; iEntry < nentries; ++iEntry) {
    if (jetChain->GetEvent(iEntry) <= 0 || skimChain->GetEvent(iEntry) <= 0 || ueChain->GetEvent(iEntry) <= 0) break;
    cout << "*****************************************************" << endl;
    cout << "Working on eventNumber:\t" << iEntry << "\tof:\t" << ntotal << endl;
    cout << "*****************************************************" << endl;
    chain->Clear();
    // Should not be null
 
    int iret = chain->Make();
    ntotal++;
  } // End event loop
  cout << "GOT THROUGH MAKE" << endl;
  chain->Finish();
  
  cout << "*************************************************" << endl;
  cout << "total number of events  " << ntotal << endl;
  cout << "*************************************************" << endl;
}

