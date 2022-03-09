#include "FairMCEventHeader.h"
#include "FairRunAna.h"
#include "TChain.h"
#include "TFile.h"
#include "TAxis.h"
#include "TF1.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdVertex.h"
#include "MpdMCTrack.h"
#include "MpdTofMatchingData.h"
#include "MpdTofHit.h"
#include "MpdKalmanFilter.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdParticle.h"
#include "MpdHelix.h"
#include "MpdTpcKalmanTrack.h"

using namespace std;

void tree() {

  // input tree setup
  TChain* inTree = new TChain("mpdsim");
  inTree->AddFile("rec.root");

  FairMCEventHeader* MCHeader = new FairMCEventHeader();
  MpdEvent* MPDEvent = new MpdEvent();
  TClonesArray* MCTracks = new TClonesArray("MpdMCTrack");
  TClonesArray* mpdKalmanTracks = new TClonesArray("MpdTpcKalmanTrack");
  TClonesArray* vertexes = new TClonesArray("MpdVertex");
  TClonesArray* tofMatchings = new TClonesArray("MpdTofMatchingData");
  TClonesArray* tofHits = new TClonesArray("MpdTofHit");

  inTree->SetBranchAddress("MCEventHeader.", &MCHeader);
  inTree->SetBranchAddress("MCTrack", &MCTracks);
  inTree->SetBranchAddress("MPDEvent.", &MPDEvent);
  inTree->SetBranchAddress("TpcKalmanTrack", &mpdKalmanTracks);
  inTree->SetBranchAddress("Vertex", &vertexes);
  inTree->SetBranchAddress("TOFMatching",&tofMatchings);
  inTree->SetBranchAddress("TOFHit",&tofHits);

  // Kalman filter setup
  TFile* inFileSim = new TFile("sim.root");
  inFileSim->Get("FairGeoParSet");
  FairRunAna ana;
  MpdKalmanFilter::Instance("KF")->Init();
  MpdTpcKalmanFilter* recoTpc = new MpdTpcKalmanFilter("TPC Kalman filter");
  recoTpc->SetSectorGeo(MpdTpcSectorGeo::Instance());
  recoTpc->FillGeoScheme();

  TF1 fMeanTPC("fMeanTPC","1./([0]*x+[2])+[1]",0.,5);
  TF1 fSigmTPC("fSigmTPC","pol5",0.,5);
  TF1 fMeanTOF("fMeanTOF","pol0",0.,5);
  TF1 fSigmTOF("fSigmTOF","pol2",0.,5);

  fMeanTPC.SetParameters(-9.351e-03,4.328e+03,-3.895e-04);
  fSigmTPC.SetParameters(338.877,-502.381,1572.066,-1857.843,917.432,-161.531);
  fMeanTOF.SetParameters(0);
  fSigmTOF.SetParameters(7.60503e-04,-4.54220e-04,2.55733e-02);
  
  // output tree setup
  TFile* convFile = new TFile("conv.root","recreate");
  TTree* convTree = new TTree("convTree","convTree");
  const Int_t nMcMax   = 10000;
  const Int_t nRcMax   = 10000;
  const Int_t nPairMax = 20000;
  Int_t   nPrimary=0;
  Int_t   nMc=0;
  Int_t   nRc=0;
  Int_t   nPair=0;
  Float_t bMc=0;
  Float_t xVtxMc=0;
  Float_t yVtxMc=0;
  Float_t zVtxMc=0;
  Float_t xVtxRc=0;
  Float_t yVtxRc=0;
  Float_t zVtxRc=0;
  Float_t vxMc[nMcMax]={0};
  Float_t vyMc[nMcMax]={0};
  Float_t vzMc[nMcMax]={0};
  Float_t vPxMc[nMcMax]={0};
  Float_t vPyMc[nMcMax]={0};
  Float_t vPzMc[nMcMax]={0};
  Int_t   vMotherId[nMcMax]={0};
  Int_t   vMcTrackId[nMcMax]={0};
  Int_t   vPdgCode[nMcMax]={0};
  Float_t vDcaX[nRcMax]={0};
  Float_t vDcaY[nRcMax]={0};
  Float_t vDcaZ[nRcMax]={0};
  Float_t vPxRc[nRcMax]={0};
  Float_t vPyRc[nRcMax]={0};
  Float_t vPzRc[nRcMax]={0};
  Float_t vPxRcRefit[nRcMax]={0};
  Float_t vPyRcRefit[nRcMax]={0};
  Float_t vPzRcRefit[nRcMax]={0};
  Int_t   vTofFlag[nRcMax]={0};
  Float_t vMass2[nRcMax]={0};
  Float_t vDedx[nRcMax]={0};
  Int_t   vMcTrackRef[nRcMax]={0};
  Short_t vNhits[nRcMax]={0};
  Bool_t  vEdgeCut[nRcMax]={0};
  Char_t  vCharge[nRcMax]={0};
  Float_t vTofWeight[nRcMax]={0};
  Bool_t  vTofTrueMatch[nRcMax]={0};
  Float_t vPairX[nPairMax]={0};
  Float_t vPairY[nPairMax]={0};
  Float_t vPairZ[nPairMax]={0};
  Float_t vPairPx[nPairMax]={0};
  Float_t vPairPy[nPairMax]={0};
  Float_t vPairPz[nPairMax]={0};
  Float_t vPairMass[nPairMax]={0};
  Float_t vPairChi2[nPairMax]={0};
  Float_t vPairOpening[nPairMax]={0};
  Int_t   vPairIndex1[nPairMax]={0};
  Int_t   vPairIndex2[nPairMax]={0};
  Float_t phiEpMc=0;

  convTree->Branch("nMc",&nMc);
  convTree->Branch("nRc",&nRc);
  convTree->Branch("nPair",&nPair);
  convTree->Branch("bMc",&bMc);
  convTree->Branch("phiEpMc",&phiEpMc);
  convTree->Branch("xVtxMc",&xVtxMc);
  convTree->Branch("yVtxMc",&yVtxMc);
  convTree->Branch("zVtxMc",&zVtxMc);
  convTree->Branch("xVtxRc",&xVtxRc);
  convTree->Branch("yVtxRc",&yVtxRc);
  convTree->Branch("zVtxRc",&zVtxRc);
  convTree->Branch("vPxMc",&vPxMc,"vPxMc[nMc]/F");
  convTree->Branch("vPyMc",&vPyMc,"vPyMc[nMc]/F");
  convTree->Branch("vPzMc",&vPzMc,"vPzMc[nMc]/F");
  convTree->Branch("vxMc",&vxMc,"vxMc[nMc]/F");
  convTree->Branch("vyMc",&vyMc,"vyMc[nMc]/F");
  convTree->Branch("vzMc",&vzMc,"vzMc[nMc]/F");
  convTree->Branch("vPdgCode",&vPdgCode,"vPdgCode[nMc]/I");
  convTree->Branch("vMotherId",&vMotherId,"vMotherId[nMc]/I");
  convTree->Branch("vMcTrackId",&vMcTrackId,"vMcTrackId[nMc]/I");
  convTree->Branch("vDcaX",&vDcaX,"vDcaX[nRc]/F");
  convTree->Branch("vDcaY",&vDcaY,"vDcaY[nRc]/F");
  convTree->Branch("vDcaZ",&vDcaZ,"vDcaZ[nRc]/F");
  convTree->Branch("vPxRc",&vPxRc,"vPxRc[nRc]/F");
  convTree->Branch("vPyRc",&vPyRc,"vPyRc[nRc]/F");
  convTree->Branch("vPzRc",&vPzRc,"vPzRc[nRc]/F");
  convTree->Branch("vPxRcRefit",&vPxRcRefit,"vPxRcRefit[nRc]/F");
  convTree->Branch("vPyRcRefit",&vPyRcRefit,"vPyRcRefit[nRc]/F");
  convTree->Branch("vPzRcRefit",&vPzRcRefit,"vPzRcRefit[nRc]/F");
  convTree->Branch("vTofFlag",&vTofFlag,"vTofFlag[nRc]/I");
  convTree->Branch("vMass2",&vMass2,"vMass2[nRc]/F");
  convTree->Branch("vDedx",&vDedx,"vDedx[nRc]/F");
  convTree->Branch("vMcTrackRef",&vMcTrackRef,"vMcTrackRef[nRc]/I");
  convTree->Branch("vNhits",&vNhits,"vNhits[nRc]/S");
  convTree->Branch("vEdgeCut",&vEdgeCut,"vEdgeCut[nRc]/O");
  convTree->Branch("vCharge",&vCharge,"vCharge[nRc]/B");
  convTree->Branch("vTofWeight",&vTofWeight,"vTofWeight[nRc]/F");
  convTree->Branch("vTofTrueMatch",&vTofTrueMatch,"vTofTrueMatch[nRc]/O");
  convTree->Branch("vPairX",&vPairX,"vPairX[nPair]/F");
  convTree->Branch("vPairY",&vPairY,"vPairY[nPair]/F");
  convTree->Branch("vPairZ",&vPairZ,"vPairZ[nPair]/F");
  convTree->Branch("vPairPx",&vPairPx,"vPairPx[nPair]/F");
  convTree->Branch("vPairPy",&vPairPy,"vPairPy[nPair]/F");
  convTree->Branch("vPairPz",&vPairPz,"vPairPz[nPair]/F");
  convTree->Branch("vPairMass",&vPairMass,"vPairMass[nPair]/F");
  convTree->Branch("vPairChi2",&vPairChi2,"vPairChi2[nPair]/F");
  convTree->Branch("vPairOpening",&vPairOpening,"vPairOpening[nPair]/F");
  convTree->Branch("vPairIndex1",&vPairIndex1,"vPairIndex1[nPair]/I");
  convTree->Branch("vPairIndex2",&vPairIndex2,"vPairIndex2[nPair]/I");

  // useful vectors
  TVector3 vRc,pRc,pRcRefit;
  TClonesArray* vPart = new TClonesArray("MpdParticle");
  vector<MpdParticle*> vpair;
  TLorentzVector p4pair;

  // cuts
  Int_t minTPChits = 10;
  Float_t maxEta =  6;
  Float_t minEta = -6;
  Float_t ptMax = 2.5;
  Float_t mEl = 0.000511;


  int tofMatchingIndices[nRcMax];

  // Loop over events
  Int_t nEvents = inTree->GetEntries();

  for (Int_t ev=0;ev<nEvents;ev++){
    inTree->GetEntry(ev);
    phiEpMc = MCHeader->GetRotZ();

    bMc    = MCHeader->GetB();
    xVtxMc = MCHeader->GetX();
    yVtxMc = MCHeader->GetY();
    zVtxMc = MCHeader->GetZ();

    MpdVertex*vertex = (MpdVertex*) vertexes->First();
    TVector3 primaryVertex;
    vertex->Position(primaryVertex);
    xVtxRc = primaryVertex.X();
    yVtxRc = primaryVertex.Y();
    zVtxRc = primaryVertex.Z();

    nRc=0;
    vPart->Clear();

    TClonesArray* MpdGlobalTracks = MPDEvent->GetGlobalTracks();
    Int_t nRcTracks = MpdGlobalTracks->GetEntriesFast();

    // build indices from reconstructed global tracks to tof matchings
    for (int i=0;i<nRcTracks;i++) tofMatchingIndices[i]=-1;
    for (int i=0;i<tofMatchings->GetEntriesFast();i++){
      MpdTofMatchingData* matching = (MpdTofMatchingData*) tofMatchings->UncheckedAt(i);
      tofMatchingIndices[matching->GetKFTrackIndex()]=i;
    }

    TBits isRcStored;

    for (Int_t i=0;i<nRcTracks;i++){
      MpdTrack*       rcTrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack* kfTrack = (MpdTpcKalmanTrack*) mpdKalmanTracks->UncheckedAt(i);
      Int_t nHits    = rcTrack->GetNofHits();
      Float_t eta    = rcTrack->GetEta();
      if (nHits<minTPChits)  continue;
      if (!(eta>minEta && eta<maxEta)) continue;
      Int_t mcTrackId = rcTrack->GetID();
      Float_t pt     = TMath::Abs(rcTrack->GetPt());
      Float_t theta  = rcTrack->GetTheta();
      Bool_t edgeCut = rcTrack->GetEdgeCut();
      Char_t charge  = rcTrack->GetCharge() > 0 ? 1 : -1;
      Int_t tofFlag  = rcTrack->GetTofFlag();
      Float_t mass2  = rcTrack->GetTofMass2();
      Float_t dedx   = rcTrack->GetdEdXTPC();
      Float_t p      = theta>1e-10 ? pt/sin(theta) : 0;
      Float_t dcax = rcTrack->GetDCAX();
      Float_t dcay = rcTrack->GetDCAY();
      Float_t dcaz = rcTrack->GetDCAZ();
      if (edgeCut) continue;
      if (!(pt>0.05)) continue;

      // check TOF matching info
      float tofWeight = -1;
      bool tofTrueMatch = 0;

      if (tofMatchingIndices[i]>=0){
        MpdTofMatchingData* matching = (MpdTofMatchingData*) tofMatchings->UncheckedAt(tofMatchingIndices[i]);
        tofWeight = matching->GetWeight();
        MpdTofHit* tofHit = (MpdTofHit*) tofHits->UncheckedAt(matching->GetTofHitIndex());
        vector<int> links;
        tofHit->getLinks(MpdTofUtils::mcTrackIndex,links);
        tofTrueMatch = 0;
        for (unsigned int iMcTrackInHit=0;iMcTrackInHit<links.size();iMcTrackInHit++){
          if (links[iMcTrackInHit]==mcTrackId) tofTrueMatch = 1;
        }
      }
      float meanTPC = fMeanTPC.Eval(p);
      float sigmTPC = fSigmTPC.Eval(p);
      float meanTOF = fMeanTOF.Eval(p);
      float sigmTOF = fSigmTOF.Eval(p);
      if (dedx-meanTPC> 4*sigmTPC) continue;
      if (dedx-meanTPC<-4*sigmTPC) continue;
      if ((tofFlag==2 || tofFlag==6) && mass2-meanTOF>4*sigmTOF) continue;
      vTofFlag[nRc]=tofFlag;
      vMass2[nRc]=mass2;
      vDedx[nRc]=dedx;
      MpdTpcKalmanTrack kfTrackRefit = *kfTrack;
      kfTrackRefit.SetDirection(MpdKalmanTrack::kInward);
      bool isRefitOk = recoTpc->Refit(&kfTrackRefit,mEl,charge,kFALSE,1);
      if (!isRefitOk) continue;
      pRcRefit = kfTrackRefit.Momentum3();
      pRc.SetPtThetaPhi(pt,theta,rcTrack->GetPhi());
      vPxRc[nRc]         = pRc[0];
      vPyRc[nRc]         = pRc[1];
      vPzRc[nRc]         = pRc[2];
      vPxRcRefit[nRc]    = pRcRefit[0];
      vPyRcRefit[nRc]    = pRcRefit[1];
      vPzRcRefit[nRc]    = pRcRefit[2];
      vDcaX[nRc]         = dcax;
      vDcaY[nRc]         = dcay;
      vDcaZ[nRc]         = dcaz;
      vMcTrackRef[nRc]   = mcTrackId;
      vNhits[nRc]        = nHits;
      vEdgeCut[nRc]      = edgeCut;
      vCharge[nRc]       = charge;
      vTofWeight[nRc]    = tofWeight;
      vTofTrueMatch[nRc] = tofTrueMatch;
      isRcStored.SetBitNumber(mcTrackId);
      MpdParticle* part = new ((*vPart)[nRc]) MpdParticle(kfTrackRefit,i,mEl);
      nRc++;
      if (nRc>nRcMax) { printf("warning: increase nRcMax\n"); return; }
    }

    // store mothers and grandmothers of stored reconstructed particles
    for (Int_t i=MCTracks->GetEntriesFast()-1;i>=0;i--){
      if (isRcStored.TestBitNumber(i)==0) continue;
      MpdMCTrack* mcTrack = (MpdMCTrack*) MCTracks->UncheckedAt(i);
      Int_t iMother = mcTrack->GetMotherId();
      if (iMother==-1) continue;
      isRcStored.SetBitNumber(iMother);
    }

    nMc=0;
    for (Int_t i=0;i<MCTracks->GetEntriesFast();i++){
      MpdMCTrack* mcTrack = (MpdMCTrack*) MCTracks->UncheckedAt(i);
      vxMc[nMc]=mcTrack->GetStartX();
      vyMc[nMc]=mcTrack->GetStartY();
      vzMc[nMc]=mcTrack->GetStartZ();
      // if (!isRcStored.TestBitNumber(i) && vxMc[nMc]*vxMc[nMc]+vyMc[nMc]*vyMc[nMc]>40.*40.) continue;
      Int_t pdgCode = mcTrack->GetPdgCode();
      Int_t motherId = mcTrack->GetMotherId();
      if (!isRcStored.TestBitNumber(i) && !(pdgCode==11 || pdgCode==-11 || pdgCode==22 || pdgCode==111 || pdgCode==221 || motherId<0)) continue;
      vPxMc[nMc]=mcTrack->GetPx();
      vPyMc[nMc]=mcTrack->GetPy();
      vPzMc[nMc]=mcTrack->GetPz();
      vPdgCode[nMc]=mcTrack->GetPdgCode();
      vMotherId[nMc]=mcTrack->GetMotherId();
      vMcTrackId[nMc]=i;
      nMc++;
      if (nMc>nMcMax) { printf("warning: increase nMcMax\n"); return; }
    }

    nPair=0;
    for (Int_t in=0;in<vPart->GetEntriesFast();in++){
      MpdParticle* pn = (MpdParticle*) vPart->UncheckedAt(in);
      if (pn->GetCharge()>0) continue;
      for (Int_t ip=0;ip<vPart->GetEntriesFast();ip++){
        MpdParticle* pp = (MpdParticle*) vPart->UncheckedAt(ip);
        if (pp->GetCharge()<0) continue;
        vpair.clear();
        vpair.push_back(new MpdParticle(*pn));
        vpair.push_back(new MpdParticle(*pp));

        MpdParticle pair;

        Double_t chi2 = pair.BuildMother(vpair);
        Double_t m = pair.GetMass();

        if (m<0 || m>0.2) {
          delete vpair[0];
          delete vpair[1];
          continue;
        }

        pRc = pair.Momentum3();
        vRc.SetXYZ(pair.Getx()(0,0),pair.Getx()(1,0),pair.Getx()(2,0));
        TVector3 p1atDCA = vpair[0]->Momentum3();
        TVector3 p2atDCA = vpair[1]->Momentum3();
        vPairX[nPair]      = vRc[0];
        vPairY[nPair]      = vRc[1];
        vPairZ[nPair]      = vRc[2];
        vPairPx[nPair]     = pRc[0];
        vPairPy[nPair]     = pRc[1];
        vPairPz[nPair]     = pRc[2];
        vPairMass[nPair]   = pair.GetMass();
        vPairChi2[nPair]   = chi2;
        vPairIndex1[nPair] = in;
        vPairIndex2[nPair] = ip;
        vPairOpening[nPair] = p1atDCA.Angle(p2atDCA);
        nPair++;
        delete vpair[0];
        delete vpair[1];
        if (nPair>nPairMax) { printf("warning: increase nPairMax\n"); return; }
      }
    }
    printf("Event: %i",ev);
    printf(" nMc=%i ",nMc);
    printf(" nRc=%i ",nRc);
    printf(" nPair=%i ",nPair);
    printf("\n");
    convTree->Fill();
  }
  convTree->Write();
  convFile->Close();
}

