//*-- Author : Victor Perevoztchikov
// 
// $Id: StInclusiveJetMaker.cxx,v 1.18 2007/10/27 17:42:59 fine Exp $
// $Log: StInclusiveJetMaker.cxx,v $
// Revision 1.18  2007/10/27 17:42:59  fine
// replace the obsolete class name
//
// Revision 1.17  2007/04/28 17:57:13  perev
// Redundant StChain.h removed
//
// Revision 1.16  2006/12/19 21:59:16  fine
// replace the class StMessMgr forward declaration with the real declaration and adjust StInclusiveJetMaker to show how to use logger
//
// Revision 1.15  2002/04/28 01:28:36  jeromel
// Reshaped comments for doxygen. Hopefully, users will propagate this good
// habit.
//
// Revision 1.14  2000/06/23 16:50:07  fisyak
// remove params
//
// Revision 1.13  1999/12/19 16:07:01  perev
// Add README
//
// Revision 1.12  1999/07/15 13:57:44  perev
// cleanup
//
// Revision 1.11  1999/07/10 22:59:16  fine
// Some comments have been introduced to show html docs
//


#include "StInclusiveJetMaker.h"
#include "TDataSetIter.h"
#include "StDAQMaker/StDAQReader.h"
#include "TChain.h"
//std
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>

//root
#include "TTree.h"
#include "TFriendElement.h"
#include "TFile.h"
#include "TSeqCollection.h"
#include "TH2F.h"
#include <map>
#include "TLorentzVector.h"
//STAR
#include "StMessMgr.h"

// StSpinPool
#include "StSpinPool/StJets/StJet.h"
#include "StSpinPool/StJets/StJets.h"
#include "StSpinPool/StJetSkimEvent/StJetSkimEvent.h"
#include "StSpinPool/StJetSkimEvent/StPythiaEvent.h"
#include "StSpinPool/StJetEvent/StJetEvent.h"
#include "StSpinPool/StJetEvent/StUeEvent.h"
#include "StSpinPool/StJetEvent/StJetVertex.h"
#include "StSpinPool/StJetEvent/StJetCandidate.h"
#include "StSpinPool/StJetEvent/StJetTrack.h"
#include "StSpinPool/StJetEvent/StJetTower.h"

ClassImp(StInclusiveJetMaker)

//_____________________________________________________________________________
/// TLA constructor
/*!
  const char *name -  the name of this constructor
  The first comment lines after the opening bracket
  ({) of a member function are considered as a member function description 
  See <A HREF="http://root.cern.ch/root/Documentation.html"> ROOT HTML documentation </A>

 */
  StInclusiveJetMaker::StInclusiveJetMaker(const char *name, const char *jetBranchname, TChain* jetChain, TChain *skimChain, int isEmbed, int isPythia, const char* outputname):StMaker(name), mjetChain(jetChain), mskimChain(skimChain), mjetname(jetBranchname), mIsEmbed(isEmbed), mIsPythia(isPythia){
    myDijetfile = TString(outputname);
}


//_____________________________________________________________________________
/// This is TLA destructor
/*!
  The first comment lines after the opening bracket
  ({) of a member function are considered as a member function description 
  
  The first comment lines after the opening bracket
  ({) of a member function are considered as a member function description 
  see: <A HREF="http://root.cern.ch/root/Documentation.html"> ROOT HTML documentation </A> 

 */
StInclusiveJetMaker::~StInclusiveJetMaker(){
  //
}


//_____________________________________________________________________________
/// Init - is a first method the top level StChain calls to initialize all its makers 
Int_t StInclusiveJetMaker::Init(){
  // Create tables
  // Create Histograms
  outfile = new TFile(myDijetfile,"recreate");
 
  jetEvent = 0;
  pjetEvent = 0;
  skimEvent =0;
  // Data and Embedding have Jet Trees and Skim Trees. Jets created from the pythia.root files only have jet trees and no skim Trees
  // The goal is to make this code versitle to take also take Pythia jets as input. So if we are using pythia jet trees mIsPythia == 1 
  if(mIsPythia == 0){
    mjetChain->SetBranchAddress("AntiKtR060NHits12", &jetEvent);
    if(mIsEmbed == 1){  mjetChain->SetBranchAddress("AntiKtR060Particle",&pjetEvent);}
    mskimChain->SetBranchAddress("skimEventBranch", &skimEvent);
  }
  else{
    mjetChain->SetBranchAddress("AntiKtR060Particle", &jetEvent); 
  }
  inclujetTree = new TTree(Form("Inclujets_%s",mjetname),Form("Inclujets_%s",mjetname));
  inclujetTree->Branch("jets",&jets,"pT/F:eta:phi:Rt:Et:sumtrackpT:sumtowerEt:detEta:dR:y:numtracks/I:numtowers:eventId:geoFlagJP1:geoFlagJP2:geoFlagAdj"); 
  inclujetTree->Branch("jets_particle",&jets_particle,"pT/F:eta:phi:Rt:Et:sumtrackpT:sumtowerEt:detEta:dR:y:numtracks/I:numtowers:eventId:geoFlagJP1:geoFlagJP2:geoFlagAdj"); 
  inclujetTree->Branch("event",&event,"flagAJP/I:flagJP1:flagJP2:eventId:verz/F:partonpT:x1:x2:prescaleJP1:prescaleJP2:prescaleAJP");
  inclujetTree->Branch("event_particle",&event_particle,"flagAJP/I:flagJP1:flagJP2:eventId:verz/F:partonpT:x1:x2:prescaleJP1:prescaleJP2:prescaleAJP"); 

  inclujetTracks = new TTree(Form("Inclujettracks_%s",mjetname),Form("Inclujettracks_%s",mjetname));
  inclujetTracks->Branch("jettracks",&jettracks,"eventId/I:pt/F:eta:phi:jetpT");

  inclujetTowers = new TTree(Form("Inclujettowers_%s",mjetname),Form("Inclujettowers_%s",mjetname));
  inclujetTowers->Branch("jettowers",&jettowers,"eventId/I:energy/F:eta:phi:jetpT");
  hnevents = new TH1F("nevents","nevents",2,0,2);
  return StMaker::Init();
}
//_____________________________________________________________________________
/// Make - this method is called in loop for each event
Int_t StInclusiveJetMaker::Make(){  
  TObjArray* jetArray = new TObjArray();
  hnevents->Fill(1.0);  
  //*******Initialization of the data or detector or pythia level dijet variables**********
  jets.pT = -1; 
  jets.eta = -100 ;
  jets.y = -100;
  jets.phi = 100;
  jets.Rt = 0;
  jets.Et = -1;
  jets.sumtrackpT = -1;
  jets.sumtowerEt = -1;
  jets.detEta = -100;
  jets.dR = 1000;
  jets.numtracks = 0;
  jets.numtowers = 0;
  jets.geoFlagJP1 = 0;
  jets.geoFlagJP2 = 0;
  jets.geoFlagAdj = 0;
 
  
  event.flagJP1 = event.flagJP2 = event.flagAJP  = 0;
  event.eventId = -1;
  event.verz = 1e3;
  event.partonpT = -1; 
  event.x1 = -1; 
  event.x2 = -1;
  event.prescaleJP1 = 1;
  event.prescaleJP2 = 1; 
  event.prescaleAJP = 1;
  // **************************************************************************************
  
  //------ Initialization of the particle level dijet variables -----
  R_min3 = R_min4 = 1000;  
  jets_particle.pT = -1; 
  jets_particle.eta = -100 ;
  jets_particle.y = -100 ;
  jets_particle.phi = 100;
  jets_particle.Rt = 0;
  jets_particle.Et = -1;
  jets_particle.sumtrackpT = -1;
  jets_particle.sumtowerEt = -1;
  jets_particle.detEta = -100;
  jets_particle.numtracks = 0;
  jets_particle.numtowers = 0;
  
 
  
  event_particle.flagJP1 = event_particle.flagJP2 = event_particle.flagAJP  = 0;
  event_particle.verz = 1e3;
  event_particle.partonpT = -1; 
  event_particle.x1 = -1; 
  event_particle.x2 = -1;
  event_particle.prescaleJP1 = 1;
  event_particle.prescaleJP2 = 1; 
  event_particle.prescaleAJP = 1;
  // -----------------------------------------
  

  if(mIsPythia == 0) {
    cout << "IN MAKE" << endl;
    if(!jetEvent) return kStOK;
    if(!skimEvent) return kStOK;
    // Enforce event synchronization
    assert(jetEvent->runId() == skimEvent->runId() && jetEvent->eventId() == skimEvent->eventId());   
    event.eventId = jetEvent->eventId();
    
    //     JP1 trigger
    StJetSkimTrig* trigJP1 = skimEvent->trigger(230410);
    if (!trigJP1) trigJP1 = skimEvent->trigger(230410);
    //     JP2 trigger
    StJetSkimTrig* trigJP2 = skimEvent->trigger(230411);
    if (!trigJP2) trigJP2 = skimEvent->trigger(230411);
    //     BHT3 trigger
    StJetSkimTrig* trigBHT3 = skimEvent->trigger(230531);
    if (!trigBHT3) trigBHT3 = skimEvent->trigger(230531);
    //     AJP trigger
    StJetSkimTrig* trigAJP = skimEvent->trigger(230420);
    if (!trigAJP) trigAJP = skimEvent->trigger(230420);
    
    event.prescaleJP1 = skimEvent->trigHeader(230410)->prescale;
    event.prescaleJP2 = skimEvent->trigHeader(230411)->prescale;
    event.prescaleAJP = skimEvent->trigHeader(230420)->prescale;
    
    barrelJetPatches0 = skimEvent->barrelJetPatchesAboveTh(0);
    barrelJetPatches = skimEvent->barrelJetPatchesAboveTh(1);
    barrelJetPatches2 = skimEvent->barrelJetPatchesAboveTh(2);
     
    if(mIsEmbed == 1){    
      const StPythiaEvent *pythiaEvent = skimEvent->mcEvent();
      event.x1 = pythiaEvent->x1();
      event.x2 = pythiaEvent->x2();
      event.partonpT = pythiaEvent->pt();
      cout << "x1: " << event.x1 << " x2: " << event.x2 << endl;
    }   
    else{
      event.x1 = 0;
      event.x2 = 0; 
      event.partonpT = 0;
    }
    
    vertex = jetEvent->vertex(0);
    // cout << "number of verties:" << jetEvent->numberOfVertices() << endl;
    // for(int i = 0; i < jetEvent->numberOfVertices(); i++){

    //   cout << "vertex " << i << " rank: " << vertex->ranking() << endl;
    // }
    if (!vertex) return kStOK;
    if (mIsEmbed == 1){ pvertex = pjetEvent->vertex();}
    
    if(vertex->ranking() < 0.0 ) return kStOK;
    //    if(vertex->numberOfJets() < 2) return kStOK; // Require at least 2 jets in each event
    event.verz = vertex->position().Z();    
    if (mIsEmbed == 1){ event_particle.verz = pvertex->position().Z();}
    
    if(trigJP1){
      if (mIsEmbed == 1 || trigJP1->didFire()){
	if(trigJP1->shouldFire() == 1){
	  event.flagJP1 = 1;
	  cout << "Check JP1" << endl;
	}
      }
    }
    if(trigJP2){
      if (mIsEmbed == 1 || trigJP2->didFire()){
	if(trigJP2->shouldFire() == 1){
	  event.flagJP2 = 1;
	  cout << "Check JP2" << endl; 
	}
      }
    }
    if(trigAJP){
      if (mIsEmbed ==1  || trigAJP->didFire()){
	if(trigAJP->shouldFire() == 1){
	  event.flagAJP = 1;
	  //cout << "Check AJP" << endl;
	}
      }
    }
        
    for (int iJet = 0; iJet < vertex->numberOfJets(); ++iJet) { //Loop over jets
      jet = vertex->jet(iJet);  
      if(event.flagAJP == 1 || event.flagJP2 == 1 || event.flagJP1 == 1){ jetArray->Add(jet);}    
    } // End loop over detector jets
        
    numEntries =  jetArray->GetEntries();

    jets.eventId = jetEvent->eventId();
    if(numEntries != 0 )
      {

	for(int iJet = 0; iJet < numEntries; iJet++){ //loop over all the jets associated with the highest ranked vertex
	  jethold = (StJetCandidate*) jetArray->At(iJet);
	  
	  
	  jet_y = 0.5*TMath::Log( (jethold->E() + jethold->pz()) / (jethold->E() - jethold->pz()) );
	  
	  
	  //---------Geometric Trigger/ Matching to a JP----------
	    if (matchedToJetPatch(jethold,barrelJetPatches)){jets.geoFlagJP1 = 1;}
	    if (matchedToJetPatch(jethold,barrelJetPatches2)){jets.geoFlagJP2 = 1;}
	    //  if (matchedToAdjacentJetPatch(jethold,barrelJetPatches0)){jets.geoFlagAdj = 1;}

	    // ----------------------------------------------------
	 
	    jets.pT = jethold->pt();
	    jets.eta = jethold->eta();
	    jets.phi = jethold->phi();
	    jets.y = 0.5*TMath::Log( (jethold->E() + jethold->pz()) / (jethold->E() - jethold->pz()) );
	    jets.Rt = jethold->rt();
	    jets.Et = jethold->E()/cosh(jethold->eta());
	    jets.sumtrackpT = jethold->sumTrackPt();
	    jets.sumtowerEt = jethold->sumTowerPt();
	    jets.detEta = jethold->detEta();
	    jets.numtracks = jethold->numberOfTracks();
	    jets.numtowers = jethold->numberOfTowers();
	    // Loop over tracks
	    for (int iTrack = 0; iTrack < jethold->numberOfTracks(); ++iTrack) {
	      jtrack = jethold->track(iTrack);
	      jettracks.eventId = jetEvent->eventId();
	      jettracks.jetpT = jethold->pt();
	      jettracks.pt = jtrack->pt();
	      jettracks.eta = jtrack->eta();
	      jettracks.phi = jtrack->phi();
	      inclujetTracks->Fill();
	      
	    }
	    for (int iTower = 0; iTower < jethold->numberOfTowers(); ++iTower) {
	      jtower = jethold->tower(iTower);
	      jettowers.eventId = jetEvent->eventId();
	      jettowers.jetpT = jethold->pt();
	      jettowers.energy = jtower->energy();
	      jettowers.eta = jtower->eta();
	      jettowers.phi = jtower->phi();
	      inclujetTowers->Fill();
	      
	    }

	 
	    if (mIsEmbed == 1){
	      cout << "The number of the jets" << pvertex->numberOfJets() <<  endl;
	      R_min3 = 1000;
	      for(int jJet = 0; jJet < pvertex->numberOfJets(); jJet++){
		cout << "jJet " << jJet << endl;
		pjet = pvertex->jet(jJet); 
		
		diffR3 = diffParticleDetector(jethold,pjet);
		cout << "dR 3: "  << R_min3 << endl; 
		if(R_min3 >  diffR3) {
		  R_min3 =  diffR3;
		  // cout << "dR 3 in loop: "  << R_min3 << endl; 
		  jets_particle.dR = R_min3;
		  jets_particle.pT = pjet->pt(); 
		  jets_particle.eta = pjet->eta() ;
		  jets_particle.phi = pjet->phi();
		  jets_particle.y = 0.5*TMath::Log( (pjet->E() + pjet->pz()) / (pjet->E() - pjet->pz()) );
		  jets_particle.Rt = pjet->rt();
		  jets_particle.Et = pjet->E()/cosh(pjet->eta());
		  jets_particle.sumtrackpT = pjet->sumTrackPt();
		  jets_particle.sumtowerEt = pjet->sumTowerPt();
		  jets_particle.detEta = pjet->detEta();
		  jets_particle.numtracks = pjet->numberOfTracks();
		  jets_particle.numtowers = pjet->numberOfTowers();
		}
	      }
	    } // End of IsEmbed if Statement	    
	    inclujetTree->Fill(); 
	} // End of Loop over Jets 
      }
  }
  else if(mIsPythia == 1){
    cout << "IN MAKE" << endl;
    assert(jetEvent);
    
    // Enforce event synchronization
    
    event.eventId = jetEvent->eventId();
    
    
    event.prescaleJP1 = 1;
    event.prescaleJP2 = 1;
    event.prescaleAJP = 1;
    
    vertex = jetEvent->vertex();
    if (!vertex) return kStOK;   
      
    // if(vertex->ranking() < 0.0 ) return kStOK;   
    //    if(vertex->numberOfJets() < 2) return kStOK; // Require at least 2 jets in each event
    event.verz = vertex->position().Z();    
    //    if (mIsEmbed == 1){ event_particle.verz = pvertex->position().Z();}
        
    event.flagAJP = 1; event.flagJP2 = 1; event.flagJP1 = 1;
    

    for (int iJet = 0; iJet < vertex->numberOfJets(); ++iJet) { //Loop over jets
      jet = vertex->jet(iJet);  
      if(event.flagAJP == 1 || event.flagJP2 == 1 || event.flagJP1 == 1){ jetArray->Add(jet);}    
    } // End loop over detector jets
    numEntries =  jetArray->GetEntries();
    
    jets.eventId = jetEvent->eventId(); 
    if(numEntries != 0 )
      {
	for(int iJet = 0; iJet < numEntries; iJet++){
	  
	  //  	    jethold = (StJetCandidate*) jetArray->At(0);
	  jethold = (StJetCandidate*) jetArray->At(iJet);
	  
	  jet_y = 0.5*TMath::Log( (jethold->E() + jethold->pz()) / (jethold->E() - jethold->pz()) );

	    
	 
  	    jets.pT = jethold->pt();
  	    jets.eta = jethold->eta();
  	    jets.phi = jethold->phi();
  	    jets.y = 0.5*TMath::Log( (jethold->E() + jethold->pz()) / (jethold->E() - jethold->pz()) );
  	    jets.Rt = jethold->rt();
  	    jets.Et = jethold->E()/cosh(jethold->eta());
  	    jets.sumtrackpT = jethold->sumTrackPt();
  	    jets.sumtowerEt = jethold->sumTowerPt();
  	    jets.detEta = jethold->detEta();
  	    jets.numtracks = jethold->numberOfTracks();
  	    jets.numtowers = jethold->numberOfTowers();
	    
	    
  	    inclujetTree->Fill(); 
	}
      }
  }
  
  else{cout << "Non valid input for mIsPythia " << endl;}
  delete jetArray;
  return kStOK;
}

  Int_t StInclusiveJetMaker::Finish(){
  outfile->Write();
  outfile->Close();
  delete outfile;
  return kStOk;  
}

void StInclusiveJetMaker::BubbleSort(TObjArray* array, Int_t size)
{

 
  int last = size-2;
  int isChanged = 1;
  while(last >=0 && isChanged)
    {
      isChanged =0;
      for (int k =0; k <= last; k++)
	{
	  
 	  StJetCandidate* jet0  = (StJetCandidate*) array->At(k);
 	  StJetCandidate* jet1  = (StJetCandidate*) array->At(k+1);
	  StJetCandidate* temp;
 	  if(jet0->pt() < jet1->pt())
 	    {    
	      temp = jet0;
	      jet0=jet1;
	      jet1=temp;
	      //  Swap(&jet0,&jet1);
 	      isChanged = 1;
	      array->AddAt(jet0, k);
	      array->AddAt(jet1, k+1);
	    }
	}
      last--;
    }
}//end bubbleSort


Double_t StInclusiveJetMaker:: diffParticleDetector(StJetCandidate *detjet, StJetCandidate *pjet){
 
  Double_t diffdR, diffphi, diffeta;
  
  diffphi = TVector2::Phi_mpi_pi(pjet->phi() - detjet->phi());
  diffeta = pjet->eta() - detjet->eta();
  diffdR  = sqrt(diffphi*diffphi + diffeta*diffeta);
  return diffdR;   
 
}

Bool_t StInclusiveJetMaker:: matchedToJetPatch(const StJetCandidate* jethold,
		       const map<int,int>& barrelJetPatches){

  for (map<int,int>::const_iterator it = barrelJetPatches.begin(); it != barrelJetPatches.end(); ++it) {
    int id = it->first;
    int adc = it->second;
    float eta, phi;
    if(StJetCandidate::getBarrelJetPatchEtaPhi(id,eta,phi)){
    float deta = jethold->detEta() - eta;
    float dphi = TVector2::Phi_mpi_pi(jethold->phi() - phi);
    if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
    }
  }
 
  return false;
}

// Bool_t StInclusiveJetMaker:: matchedToAdjacentJetPatch(const StJetCandidate* jethold,
// 		       const map<int,int>& barrelJetPatches){
//   //  cout << "Size: " << barrelJetPatches.size() <<endl;
//   const int num = 18;
//   int jpID[num] = {99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99};
//   int i= 0;
//   Double_t aveEta, avePhi = 0;
//   for (map<int,int>::const_iterator it = barrelJetPatches.begin(); it != barrelJetPatches.end(); ++it) {
//     int id = it->first;
//     int adc = it->second;
//     jpID[i] = id; 
//     i++;
//   }
//   if( CheckAdjacent(jpID,barrelJetPatches.size(),aveEta,avePhi, jethold)== true) return true;
//    else{
//     return false;
//   }
// }
// Bool_t StInclusiveJetMaker::CheckAdjacent(Int_t jpID[], Int_t num, Double_t aveEta, Double_t avePhi, const StJetCandidate* jethold2)
// {
//   //cout << "jethold2 eta :  " << jethold2->detEta() << " phi: " << jethold2->phi() <<   endl;
//   float eta1, phi1 = 0;
//   float eta2, phi2 = 0;
//   for(int i =0; i<num; i++) // loop through array of JP ids that are above thr
//     {
//       //cout << "jpID: " << jpID[i] << endl;
//       for(int j = 0; j < num; j++){ // loop to compare elements within the array
// 	if(jpID[i] > 0 && jpID[i] < 5){ // JPs  1-4
// 	  if(jpID[j]== jpID[i]+6){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0); // convert back to -pi to pi
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }	    
// 	  }
// 	  if(jpID[j]==jpID[i]+12){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	  if(jpID[j]==jpID[i]+1){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	    }
// 	  }
// 	  if(jpID[j]==jpID[i]-1){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	}
// 	if(jpID[i] > 6 && jpID[i] < 11){ // JPs  7-10
// 	  if(jpID[j]== jpID[i]+6){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	  if(jpID[j]==jpID[i]-6){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	  if(jpID[j]==jpID[i]+1){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	  if(jpID[j]==jpID[i]-1){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	}
// 	if(jpID[i] > 12  && jpID[i] < 17){ // JPs  13-16
// 	  if(jpID[j]== jpID[i]-6){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	  if(jpID[j]==jpID[i]-12){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	  if(jpID[j]==jpID[i]+1){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	  if(jpID[j]==jpID[i]-1){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	}
//       }
//       if(jpID[i] == 0){
// 	for(int j = 0; j < num; j++){ // loop to compare elements within the array
// 	  if(jpID[j] == 6){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	    aveEta  = (eta1+ eta2) / 2.0 ;
// 	    phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	    phi2 = TVector2::Phi_0_2pi(phi2); 
// 	    avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	    //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	    float deta = jethold2->detEta() - aveEta;
// 	    float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	    //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	    if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
// 	  if(jpID[j] == 12){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  } 
// 	  if(jpID[j] == 5){ 
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }	
// 	}// end of id ==0 
// 	if(jpID[i] == 5){
// 	for(int j = 0; j < num; j++){ // loop to compare elements within the array
// 	  if(jpID[j] == 5){
// 	 if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){  
// 	    aveEta  = (eta1+ eta2) / 2.0 ;
// 	    avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	    phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	    phi2 = TVector2::Phi_0_2pi(phi2); 
// 	    //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	    float deta = jethold2->detEta() - aveEta;
// 	    float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	    //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	    if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	  }
// 	  if(jpID[j] == 17){
	  
// 	    aveEta  = (eta1+ eta2) / 2.0 ;
// 	    phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	    phi2 = TVector2::Phi_0_2pi(phi2); 
// 	    avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	    //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	    float deta = jethold2->detEta() - aveEta;
// 	    float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	    //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	    if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	  } 	
// 	}
//       }// end of id == 5
//       if(jpID[i] == 6){
// 	for(int j = 0; j < num; j++){ // loop to compare elements within the array
// 	  if(jpID[j] == 11){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	    //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }
	 
// 	if(jpID[j] == 12){ 
// 	  if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	  aveEta  = (eta1+ eta2) / 2.0 ;
// 	  phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	  phi2 = TVector2::Phi_0_2pi(phi2); 
// 	  avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	  //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	  float deta = jethold2->detEta() - aveEta;
// 	  float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	  //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	  if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	  }
// 	}	
// 	}
//       }// end of id == 6 
      
//       if(jpID[i] == 11){
// 	for(int j = 0; j < num; j++){ // loop to compare elements within the array
// 	  if(jpID[j] == 17){
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      //cout << "phi1: " << phi1 << "phi2: " << phi2 << endl;
// 	      //cout << "eta1: " << eta1 << "eta2: " << eta2 << endl;
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi(( phi1 +phi2) / 2);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  } 	  
// 	}
//       }// end of id ==11 
//       if(jpID[i] == 12){
// 	for(int j = 0; j < num; j++){ // loop to compare elements within the array
// 	  if(jpID[j] == 17){ 
// 	    if(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[i],eta1,phi1) &&(StJetCandidate::getBarrelJetPatchEtaPhi(jpID[j],eta2,phi2))){ 
// 	      aveEta  = (eta1+ eta2) / 2.0 ;
// 	      phi1 = TVector2::Phi_0_2pi(phi1); // convert from 0 to 2pi
// 	      phi2 = TVector2::Phi_0_2pi(phi2); 
// 	      avePhi = TVector2::Phi_mpi_pi((phi1 + phi2) / 2.0);
// 	      //cout << "<eta>: " <<  aveEta << "<phi>: " << avePhi << endl;
// 	      float deta = jethold2->detEta() - aveEta;
// 	      float dphi = TVector2::Phi_mpi_pi(jethold2->phi() - avePhi);
// 	      //cout << "deta: " << deta <<"dphi: " << dphi << endl;
// 	      if (fabs(deta) < 0.6 && fabs(dphi) < 0.6) return true;
// 	    }
// 	  }	
// 	}
//       }// end of id ==12
//       if(jpID[i] == 17){
// 	//cout << "Should have been determine " << endl; 
//       } 
      
//     }
//   return false;
//   //cout << "testing" << endl;
// }



