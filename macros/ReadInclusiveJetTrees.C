///////////////////////////////////////////////////////////////
//
// QualityAnalysis.C
//   Macro for processing the dijet trees created from StDijetMaker and 
//   creating the necessary plots to compare Data and MC.
//
//   Author: Grant Webb, University of Kentucky
//           Jan 27, 2010
//
//////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "TCanvas.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TProfile.h"
#include "TArrayD.h"
#include "TCut.h"
#include "TPolyMarker.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include <Stiostream.h>
#endif

// Main Routine:
//void ReadDiJetTrees(const char* infile, const char* outfile);

//Helper Functions:

// Global Parameters: 
const int nfip=128;
char filename[nfip];


static const Int_t npTbins = 13; // jet pT 18% resolutions
const Double_t npTbins0[npTbins+1] = {10.0,11.8,13.9,16.4,19.4,22.9,27.0,31.9,37.6,44.4,52.3,61.8,72.9,86.0};

//static const Int_t npTbins = 25;
//const Double_t npTbins0[npTbins+1] = {10.0,10.9,11.8,12.8,14.0,15.2,16.5,17.9,19.5,21.2,23.0,25.0,27.2,29.2,32.1,34.9,38.0,41.3,44.9,48.8,53.0,57.6,62.6,68.1,74.0,80.4};

void ReadInclusiveJetTrees(
			   //const char* infile = "/gpfs01/star/i_uky/gdwebb/InclusiveJetAnalysis/500GeV/R060/InclusivejetTrees/inclujetAna_10103041.inclujets.root",
			   const char* infile = "/star/u/gdwebb/InclusiveJet/test500GeVjets10103041data.root",
			   const char* outname = "testInclusive500GeVjets10103041.root"){
  
  struct protoDijet {
    
    float pT,corr_pT,eta,phi,Rt,Et,sumtrackpT,sumtowerEt,detEta,dR,y,area,ueDensity;
    int numtracks ,numtowers,eventId,geoFlagJP1,geoFlagJP2,geoFlagAdj;    
  } jets;
  struct events{
    Int_t flagAJP, flagJP1, flagJP2, eventId;
    Float_t verz, partonpT, x1,x2,prescaleJP1,prescaleJP2,prescaleAJP;
    
  } event;
  struct jetTracks{
    Int_t eventId;
    Float_t pt, eta, phi, jetpT;
  } jettracks;
  struct jetTowers{
    Int_t eventId;
    Float_t energy, eta, phi, jetpT;
  } jettowers;
  
  TFile *outfile = new TFile(outname,"RECREATE");
  TChain jetchain("Inclujets_AntiKtR060NHits12");
  jetchain.Add(infile);  
  
  TChain trackchain("Inclujettracks_AntiKtR060NHits12");
  trackchain.Add(infile); 
  
  TChain towerchain("Inclujettowers_AntiKtR060NHits12");
  towerchain.Add(infile); 

  // Set Branch Addresses 
  jetchain.SetBranchAddress("jets",&jets); 
  jetchain.SetBranchAddress("event",&event);
  
  trackchain.SetBranchAddress("jettracks",&jettracks);
  trackchain.SetBranchAddress("jettowers",&jettowers);
 
  // Create Data histograms 
  

  //---------Before Cuts-----------------
  TH1F* hjetpt_bfCut = new TH1F("hjetpt_bfCut","hjetpt_bfCut",npTbins,npTbins0);
  TH1F* hjet_corrpt_bfCut = new TH1F("hjet_corrpt_bfCut","hjet_corrpt_bfCut",npTbins,npTbins0);
  TH1F* hjeteta_bfCut = new TH1F("hjeteta_bfCut","hjeteta_bfCut",66,-1.5,1.5);
  TH1F* hjety_bfCut = new TH1F("hjety_bfCut","hjety_bfCut",66,-1.5,1.5);
  TH1F* hjetphi_bfCut = new TH1F("hjetphi_bfCut","hjetphi_bfCut",120,-3.14159,3.14159);
  TH1F* hjetRt_bfCut = new TH1F("hjetRt_bfCut","hjetRt_bfCut",24,0,1.2);
  TH1F* hjetnumTracks_bfCut = new TH1F("hjetnumTracks_bfCut","hjetnumTracks_bfCut",80,0,80);
  TH1F* hjetnumTowers_bfCut = new TH1F("hjetnumTowers_bfCut","hjetnumTowers_bfCut",100,0,100);
  TH1F* hjetsumTrackpT_bfCut = new TH1F("hjetsumTrackpT_bfCut","hjetsumTrackpT_bfCut",100,0,100);  
  TH1F* hjetsumTowerEt_bfCut = new TH1F("hjetsumTowerEt_bfCut","hjetsumTowerEt_bfCut",100,0,100);  

  TH1F* hverZ_bfCut = new TH1F("hverZ_bfCut","hverZ_bfCut",240,-120,120);
  TH2F* hdetEtaVsjety = new TH2F("hdetEtaVsjety","hdetEtaVsjety",66,-1.0,1.0,66,-1.0,1.0); 
  //-------------------------------------

  //***********After Cuts****************
  TH1F* hdphi_afCut = new TH1F("hdphi_afCut","hdphi_afCut",100,-6,6);  
  TH1F* hverZ_afCut = new TH1F("hverZ_afCut","hverZ_afCut",240,-120,120);
  TH1F* hjetpt_afCut = new TH1F("hjetpt_afCut","hjetpt_afCut",npTbins,npTbins0);
  TH1F* hjet_corrpt_afCut = new TH1F("hjet_corrpt_afCut","hjet_corrpt_afCut",npTbins,npTbins0);
  TH1F* hjeteta_afCut = new TH1F("hjeteta_afCut","hjeteta_afCut",66,-1.5,1.5);
  TH1F* hjety_afCut = new TH1F("hjety_afCut","hjety_afCut",66,-1.5,1.5);
  TH1F* hjetphi_afCut = new TH1F("hjetphi_afCut","hjetphi_afCut",120,-3.14159,3.14159);
  TH1F* hjetRt_afCut = new TH1F("hjetRt_afCut","hjetRt_afCut",24,0,1.2);
  TH1F* hjetnumTracks_afCut = new TH1F("hjetnumTracks_afCut","hjetnumTracks_afCut",80,0,80);
  TH1F* hjetnumTowers_afCut = new TH1F("hjetnumTowers_afCut","hjetnumTowers_afCut",100,0,100);
  TH1F* hjetsumTrackpT_afCut = new TH1F("hjetsumTrackpT_afCut","hjetsumTrackpT_afCut",100,0,100);  
  TH1F* hjetsumTowerEt_afCut = new TH1F("hjetsumTowerEt_afCut","hjetsumTowerEt_afCut",100,0,100);  
  //************************************

  //---------Trigger Histograms------------


  TH1F *hprescaleJP1 = new TH1F("hprescaleJP1","hprescaleJP1",200,-100,100);
  TH1F *hprescaleJP2 = new TH1F("hprescaleJP2","hprescaleJP2",200,-100,100);
  TH1F *hprescaleAJP = new TH1F("hprescaleAJP","hprescaleAJP",200,-100,100);

  TH1F* hVerZJP1 = new TH1F("hVerZJP1","hVerZJP1",120,-60,60);
  TH1F* hVerZJP2 = new TH1F("hVerZJP2","hVerZJP2",120,-60,60);
  TH1F* hVerZAJP = new TH1F("hVerZAJP","hVerZAJP",120,-60,60);
  
  TH1F* hVerZJP1_noCut = new TH1F("hVerZJP1_noCut","hVerZJP1_noCut",240,-120,120);
  TH1F* hVerZJP2_noCut = new TH1F("hVerZJP2_noCut","hVerZJP2_noCut",240,-120,120);
  TH1F* hVerZAJP_noCut = new TH1F("hVerZAJP_noCut","hVerZAJP_noCut",240,-120,120);

  TH2F* heta_yvspTJP1 = new TH2F("heta_yvspTJP1","heta_yvspTJP1",npTbins,npTbins0,66,-1.0,1.0);  
  TH2F* heta_yvspTJP2 = new TH2F("heta_yvspTJP2","heta_yvspTJP2",npTbins,npTbins0,66,-1.0,1.0);  
  TH2F* heta_yvspTAJP = new TH2F("heta_yvspTAJP","heta_yvspTAJP",npTbins,npTbins0,66,-1.0,1.0);  

  TH1F* hjetptJP1 = new TH1F("hjetptJP1","hjetptJP1",npTbins,npTbins0);
  TH1F* hjetetaJP1 = new TH1F("hjetetaJP1","hjetetaJP1",66,-1.5,1.5);
  TH1F* hjetyJP1 = new TH1F("hjetyJP1","hjetyJP1",66,-1.5,1.5);
  TH1F* hjetphiJP1 = new TH1F("hjetphiJP1","hjetphiJP1",120,-3.14159,3.14159);
  TH1F* hjetRtJP1 = new TH1F("hjetRtJP1","hjetRtJP1",24,0,1.2);
  TH1F* hjetnumTracksJP1 = new TH1F("hjetnumTracksJP1","hjetnumTracksJP1",80,0,80);
  TH1F* hjetnumTowersJP1 = new TH1F("hjetnumTowersJP1","hjetnumTowersJP1",100,0,100);
  TH1F* hjetsumTrackpTJP1 = new TH1F("hjetsumTrackpTJP1","hjetsumTrackpTJP1",100,0,100);  
  TH1F* hjetsumTowerEtJP1 = new TH1F("hjetsumTowerEtJP1","hjetsumTowerEtJP1",100,0,100);
   
  TH1F* hjetptJP2 = new TH1F("hjetptJP2","hjetptJP2",npTbins,npTbins0);
  TH1F* hjet_corrptJP2 = new TH1F("hjet_corrptJP2","hjet_corrptJP2",npTbins,npTbins0);
  TH1F* hjetetaJP2 = new TH1F("hjetetaJP2","hjetetaJP2",66,-1.5,1.5);
  TH1F* hjetyJP2 = new TH1F("hjetyJP2","hjetyJP2",66,-1.5,1.5);
  TH1F* hjetphiJP2 = new TH1F("hjetphiJP2","hjetphiJP2",120,-3.14159,3.14159);
  TH1F* hjetRtJP2 = new TH1F("hjetRtJP2","hjetRtJP2",24,0,1.2);
  TH1F* hjetnumTracksJP2 = new TH1F("hjetnumTracksJP2","hjetnumTracksJP2",80,0,80);
  TH1F* hjetnumTowersJP2 = new TH1F("hjetnumTowersJP2","hjetnumTowersJP2",100,0,100);
  TH1F* hjetsumTrackpTJP2 = new TH1F("hjetsumTrackpTJP2","hjetsumTrackpTJP2",100,0,100);
  TH1F* hjetsumTowerEtJP2 = new TH1F("hjetsumTowerEtJP2","hjetsumTowerEtJP2",100,0,100);
 
  TH1F* hjetptAJP = new TH1F("hjetptAJP","hjetptAJP",npTbins,npTbins0);
  TH1F* hjetetaAJP = new TH1F("hjetetaAJP","hjetetaAJP",66,-1.5,1.5);
  TH1F* hjetyAJP = new TH1F("hjetyAJP","hjetyAJP",66,-1.5,1.5);
  TH1F* hjetphiAJP = new TH1F("hjetphiAJP","hjetphiAJP",120,-3.14159,3.14159);
  TH1F* hjetRtAJP = new TH1F("hjetRtAJP","hjetRtAJP",24,0,1.2);
  TH1F* hjetnumTracksAJP = new TH1F("hjetnumTracksAJP","hjetnumTracksAJP",80,0,80);
  TH1F* hjetnumTowersAJP = new TH1F("hjetnumTowersAJP","hjetnumTowersAJP",100,0,100);
  TH1F* hjetsumTrackpTAJP = new TH1F("hjetsumTrackpTAJP","hjetsumTrackpTAJP",100,0,100);
  TH1F* hjetsumTowerEtAJP = new TH1F("hjetsumTowerEtAJP","hjetsumTowerEtAJP",100,0,100);  
  TH1F* hjetTrackpt = new TH1F("hjetTrackpt","hjetTrackpt",150,0,150); 
 
    //    if(TMath::Abs(jets3.eta) > 0.8 || TMath::Abs(jets4.eta) > 0.8 ) continue; // detector eta cut
  // --------------------------------------
  
  Int_t nevents = jetchain.GetEntries();
  cout << "nevents: " << nevents << endl;
  for(Int_t i=0; i< nevents;i++){
    
    if(i%1000 == 0 )cout << "nevent: " << i << endl; 
    
    jetchain.GetEvent(i);
    
    
    
   
    jetResol = -100;
    //    cout << "jetRatio: " << jetRatio << endl;
    //---------Filling the Before Cuts Histograms----------


   
   
    hjetpt_bfCut->Fill(jets.pT);
    hjet_corrpt_bfCut->Fill(jets.corrpT);
    hjeteta_bfCut->Fill(jets.eta);
    hjety_bfCut->Fill(jets.y);
    hjetphi_bfCut->Fill(jets.phi);
    hjetRt_bfCut->Fill(jets.Rt);
    hjetnumTracks_bfCut->Fill(jets.numtracks);
    hjetnumTowers_bfCut->Fill(jets.numtowers);
    hjetsumTrackpT_bfCut->Fill(jets.sumtrackpT);
    hjetsumTowerEt_bfCut->Fill(jets.sumtowerEt);   
    hverZ_bfCut->Fill(event.verz); 
   
    if(TMath::Abs(jets.detEta) > 0.7) continue; // detector eta cut
    
    if(TMath::Abs(event.verz) < 50){
      hdetEtaVsjety->Fill(jets.detEta,jets.y);
    }

    if(TMath::Abs(jets.y) > 0.8 ) continue; // rapidity cut

    if(jets.Rt > 0.95) continue; // Rt Cut
    
    //*****************_VERTEZ Z BEFORE CUT *********************
    if(event.flagJP1 ==1){
      if(jets.geoFlagJP1 == 1 ){ // matched to a least one JP ;
	hVerZJP1_noCut->Fill(event.verz);
      }
    }
    if(event.flagJP2 ==1){
      if(jets.geoFlagJP2 == 1){ // matched to a least one JP ;
	hVerZJP2_noCut->Fill(event.verz);
      }
    }
    if(event.flagAJP ==1){
      if(jets.geoFlagAdj == 1){ // matched to a least one JP ;
	hVerZAJP_noCut->Fill(event.verz);
      }
    }
    // **********************************************************
    
    if(TMath::Abs(event.verz) > 50) continue; // Z vertex Cut  
    //    if(jets.pT < 11.0) continue; // jet pT cut on the data. 
    



    // Trigger Requirements    
    if(event.flagJP1== 1){
      //cout << "JP1 fired" << endl;
      if(jets.geoFlagJP1 == 1){ // matched to a least one JP ;
	

	hprescaleJP1->Fill(event.prescaleJP1);

	hjetptJP1->Fill(jets.pT);
	hjetetaJP1->Fill(jets.eta);
	hjetyJP1->Fill(jets.y);
	hjetphiJP1->Fill(jets.phi);
	hjetRtJP1->Fill(jets.Rt);
	hjetnumTracksJP1->Fill(jets.numtracks);
	hjetnumTowersJP1->Fill(jets.numtowers);
	hjetsumTrackpTJP1->Fill(jets.sumtrackpT);
	hjetsumTowerEtJP1->Fill(jets.sumtowerEt);

	hVerZJP1->Fill(event.verz);
	heta_yvspTJP1->Fill(jets.pT,jets.y - jets.eta);

      }
    }
    
    if(event.flagJP2 == 1 ){
      //      cout << "JP2 fired" <<  endl;
      if(jets.geoFlagJP2 == 1){ // matched to a least one JP ;


	hprescaleJP2->Fill(event.prescaleJP2);

	hjetptJP2->Fill(jets.pT);
	hjet_corrptJP2->Fill(jets.corrpT);
	hjetetaJP2->Fill(jets.eta);
	hjetyJP2->Fill(jets.y);
	hjetphiJP2->Fill(jets.phi);
	hjetRtJP2->Fill(jets.Rt);
	hjetnumTracksJP2->Fill(jets.numtracks);
	hjetnumTowersJP2->Fill(jets.numtowers);
	hjetsumTrackpTJP2->Fill(jets.sumtrackpT);
	hjetsumTowerEtJP2->Fill(jets.sumtowerEt);
	hVerZJP2->Fill(event.verz);
	heta_yvspTJP2->Fill(jets.pT,jets.y - jets.eta);	

   


      } // End of GeoJP2 trigger
    } // End of soft/hard trigger
    if(event.flagAJP == 1){
      //      cout << "AJP fired" << endl;
   
      if(jets.geoFlagAdj == 1){ // matched to a least one JP ;
	

	hprescaleAJP->Fill(event.prescaleAJP);

	hjetptAJP->Fill(jets.pT);
	hjetetaAJP->Fill(jets.eta);
	hjetyAJP->Fill(jets.y);
	hjetphiAJP->Fill(jets.phi);
	hjetRtAJP->Fill(jets.Rt);
	hjetnumTracksAJP->Fill(jets.numtracks);
	hjetnumTowersAJP->Fill(jets.numtowers);
	hjetsumTrackpTAJP->Fill(jets.sumtrackpT);
	hjetsumTowerEtAJP->Fill(jets.sumtowerEt);

	hVerZAJP->Fill(event.verz);
	heta_yvspTJP2->Fill(jets.pT,jets.y - jets.eta);
      }
    }
   
    // Loop over track Chain
    for (int j =0; j <  trackchain.GetEntries(); j++){
      trackchain.GetEvent(j);
      
      if(jets.eventId == jettracks.eventId && jets.pT == jettracks.jetpT){
	cout << "track jet pT: " << jettracks.jetpT <<  "jets pT " <<  jets.pT<< endl;
	test ++;
	hjetTrackpt->Fill(jettracks.pt);
      }  
    } // End of Loop over tracks
    
    // Loop over tower Chain
    for (int j =0; j <  towerchain.GetEntries(); j++){
      towerchain.GetEvent(j);
      
      if(jets.eventId == jettowers.eventId && jets.pT == jettowers.jetpT){
	cout << "tower jet pT: " << jettowers.jetpT <<  "jets energy" <<  jets.energy<< endl;
	test ++;
	hjetTrackpt->Fill(jettracks.pt);
      }  
    } // End of Loop over tower
    //---------Filling the After Cuts Histograms----------

    
    hverZ_afCut->Fill(event.verz);
    // cout << "jetpTawaysideAJP: " << jetpTawaySideJP1 << endl;
   
    hjetpt_afCut->Fill(jets.pT);
    hjet_corrpt_afCut->Fill(jets.corrpT);
    hjeteta_afCut->Fill(jets.eta);
    hjety_afCut->Fill(jets.y);
    hjetphi_afCut->Fill(jets.phi);
    hjetRt_afCut->Fill(jets.Rt);
    hjetnumTracks_afCut->Fill(jets.numtracks);
    hjetnumTowers_afCut->Fill(jets.numtowers);
    hjetsumTrackpT_afCut->Fill(jets.sumtrackpT);
    hjetsumTowerEt_afCut->Fill(jets.sumtowerEt);    
    //-----------------------------------------------------
  }
  outfile->Write();
  outfile->Close();
  delete outfile;
}

void comboALL(TFile *t[13],char pta[100],TH1F *totH,double Norm[12]){

  float MinbXsec=28.12;
  const int lowP=1;
  TH1F *pt[13];
  for (int i=lowP;i<13;i++){
    pt[i]=(TH1F*) t[i]->Get(pta);
    // pt[i]->Scale(1);
    pt[i]->Sumw2();
    //    cout << "number bins: " << pt[i]->GetXaxis()->GetNbins() << endl;
    pt[i]->Scale(1/Norm[i]);
    totH->Add(pt[i]);
  }
}


