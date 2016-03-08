// $Id: StInclusiveJetMaker.h,v 1.15 2003/09/10 19:47:43 perev Exp $

#ifndef STAR_StInclusiveJetMaker_hh
#define STAR_StInclusiveJetMaker_hh

/*!
 *                                                                     
 * \class  StInclusiveJetMaker
 * \author Grant Webb
 * \date   2010/10/26
 * \brief  virtual base class for Maker
 *

 *
 */                                                                      

#ifndef StMaker_H
#include "StMaker.h"
#endif

class StJetCandidate;
class StJetEvent;
class StUeEvent;
class StJetSkimEvent;
class StJetTrack;
class StJetTower;
class StPythiaEvent;
class StJetVertex;
class TChain;
class TH1F;
class TH2F;
// You may forward declare other classes if you have data-members
// used in pointer-only context by using declaration like
// class St_SomeExternClass;
//
// You do need in such simple case to add the include file
// (and compilation is much faster).

class StInclusiveJetMaker : public StMaker {
 private:
  // Private method declaration if any
 
 protected:
  // Protected method if any

 public: 
  StInclusiveJetMaker(const char *name,const char *jetBranchname, TChain *jetChain, TChain* skimChain, TChain* ueChain, int isEmbed, int isPythia, const char *outputfile);
  virtual       ~StInclusiveJetMaker();
  virtual Int_t Init();
  virtual Int_t  Make();
  virtual Int_t  Finish();
  void BubbleSort(TObjArray *array, Int_t size);
  
  TChain* mjetChain; 
  TChain* mskimChain;
  TChain* mUeChain;
  const char *mjetname;
  int mIsEmbed,mIsPythia;
  StJetEvent* jetEvent, *pjetEvent;
  StJetSkimEvent* skimEvent;
  StUeEvent *ueEvent_transP, *ueEvent_transM;
  StUeEvent *pueEvent_transP, *pueEvent_transM;
  StJetTrack * jtrack;
  StJetTower *jtower;
  StJetVertex* vertex, *pvertex;
  TTree *inclujetTree, *inclujetTracks, *inclujetTowers;
  TFile * outfile;
  TString myDijetfile;
  TH1F *hnevents;
  StJetCandidate *jet, *jethold , *pjet;
  Double_t jetmass3sq, jetmass4sq, R_min3, R_min4, diffR3,diffR4, jet_y;
  Int_t flagPjet3, flagPjet4;
  map<int,int> barrelJetPatches;
  map<int,int> barrelJetPatches2;
  map<int,int> barrelJetPatches0;
  /* struct Dijet{ */

  /*   StJetCandidate *jet3, *jet4; */

  /* } dijet, pdijet; */

  struct protoDijet {
    
    float pT, corr_pT,eta, phi, Rt, Et, sumtrackpT, sumtowerEt, detEta,dR, y, area, ueDensity;
    int numtracks, numtowers, eventId, geoFlagJP1, geoFlagJP2, geoFlagAdj;
    
  } jets,jets_particle ;
  struct events{
    Int_t flagAJP, flagJP1, flagJP2, eventId;
    Float_t verz, partonpT, x1,x2,prescaleJP1,prescaleJP2,prescaleAJP;

  } event, event_particle;

  struct jetTracks{
    Int_t eventId;
    Float_t pt, eta, phi, jetpT;
  } jettracks;
  struct jetTowers{
    Int_t eventId;
    Float_t energy, eta, phi, jetpT;
  } jettowers;
  Double_t diffParticleDetector(StJetCandidate *detjet, StJetCandidate *pjet);
  Bool_t matchedToJetPatch(const StJetCandidate *jet, const map<int,int>& barrelJetPatches);
  Bool_t matchedToAdjacentJetPatch(const StJetCandidate *jet, const map<int,int>& barrelJetPatches);
  Bool_t CheckAdjacent(Int_t jpID[], Int_t num, Double_t aveEta, Double_t avePhi, const StJetCandidate *jettest);
  float dphi, dphiNoSort,deta, regionUEarea;
  Int_t numEntries;
  /// Displayed on session exit, leave it as-is please ...
  /* virtual const char *GetCVS() const { */
  /*   static const char cvs[]="Tag $Name:  $ $Id: StInclusiveJetMaker.h,v 1.15 2003/09/10 19:47:43 perev Exp $ built "__DATE__" "__TIME__ ;  */
  /*   return cvs; */
  /* } */

  ClassDef(StInclusiveJetMaker,0)   //StAF chain virtual base class for Makers
};

#endif


// $Log: StInclusiveJetMaker.h,v $
// Revision 1.15  2003/09/10 19:47:43  perev
// ansi corrs
//
// Revision 1.14  2002/11/26 23:49:40  jeromel
// Small modif after Art's note ... doxygen issue + cleanup
//
// Revision 1.13  2002/04/28 01:28:36  jeromel
// Reshaped comments for doxygen. Hopefully, users will propagate this good
// habit.
//
// Revision 1.12  1999/09/24 22:03:09  perev
// Add InitRun & FinishRun to template maker
//
