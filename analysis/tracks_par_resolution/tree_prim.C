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
#include "TH2D.h"
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
#include "MpdKfV0Fitter.h"

using namespace std;

//--------------------------------------------------

void tree_prim()
{

  // input tree setup
  TChain *inTree = new TChain("mpdsim");
  inTree->AddFile("rec.root");

  FairMCEventHeader *MCHeader = new FairMCEventHeader();
  MpdEvent *MPDEvent = new MpdEvent();
  TClonesArray *MCTracks = new TClonesArray("MpdMCTrack");
  TClonesArray *mpdKalmanTracks = new TClonesArray("MpdTpcKalmanTrack");
  TClonesArray *vertexes = new TClonesArray("MpdVertex");
  TClonesArray *tofMatchings = new TClonesArray("MpdTofMatchingData");
  TClonesArray *tofHits = new TClonesArray("MpdTofHit");

  inTree->SetBranchAddress("MCEventHeader.", &MCHeader);
  inTree->SetBranchAddress("MCTrack", &MCTracks);
  inTree->SetBranchAddress("MPDEvent.", &MPDEvent);
  inTree->SetBranchAddress("TpcKalmanTrack", &mpdKalmanTracks);
  inTree->SetBranchAddress("Vertex", &vertexes);
  inTree->SetBranchAddress("TOFMatching", &tofMatchings);
  inTree->SetBranchAddress("TOFHit", &tofHits);

  // Kalman filter setup
  TFile *inFileSim = new TFile("sim.root");
  inFileSim->Get("FairGeoParSet");
  FairRunAna ana;
  MpdKalmanFilter::Instance("KF")->Init();
  MpdTpcKalmanFilter *recoTpc = new MpdTpcKalmanFilter("TPC Kalman filter");
  recoTpc->SetSectorGeo(MpdTpcSectorGeo::Instance());
  recoTpc->FillGeoScheme();

  TF1 fMeanTPC("fMeanTPC", "1./([0]*x+[2])+[1]", 0., 5);
  TF1 fSigmTPC("fSigmTPC", "pol5", 0., 5);
  TF1 fMeanTOF("fMeanTOF", "pol0", 0., 5);
  TF1 fSigmTOF("fSigmTOF", "pol2", 0., 5);

  fMeanTPC.SetParameters(-1.12616e-02, 4.22476e+03, -1.92951e-04);
  fSigmTPC.SetParameters(295.256, -667.891, 2341.81, -3323.45, 2060.76, -470.015);
  fMeanTOF.SetParameters(0);
  fSigmTOF.SetParameters(-6.32897e-06, 0.00289299, 0.0296647);

  // output tree setup
  TFile *convFile = new TFile("conv.prim.root", "recreate");
  TTree *convTree = new TTree("convTree", "convTree");
  const Int_t nMcMax = 10000;
  const Int_t nRcMax = 10000;
  const Int_t nPairMax = 100;
  Int_t nPrimary = 0;
  Int_t nMc = 0;
  Int_t nRc = 0;
  Int_t nPair = 0;
  Float_t bMc = 0;
  Float_t xVtxMc = 0;
  Float_t yVtxMc = 0;
  Float_t zVtxMc = 0;
  Float_t xVtxRc = 0;
  Float_t yVtxRc = 0;
  Float_t zVtxRc = 0;
  Float_t vxMc[nMcMax] = {0};
  Float_t vyMc[nMcMax] = {0};
  Float_t vzMc[nMcMax] = {0};
  Float_t vPxMc[nMcMax] = {0};
  Float_t vPyMc[nMcMax] = {0};
  Float_t vPzMc[nMcMax] = {0};
  Int_t vMotherId[nMcMax] = {0};
  Int_t vMcTrackId[nMcMax] = {0};
  Int_t vPdgCode[nMcMax] = {0};
  Float_t vDcaX[nRcMax] = {0};
  Float_t vDcaY[nRcMax] = {0};
  Float_t vDcaZ[nRcMax] = {0};
  Float_t vPxRc[nRcMax] = {0};
  Float_t vPyRc[nRcMax] = {0};
  Float_t vPzRc[nRcMax] = {0};
  Float_t vPxRcRefit[nRcMax] = {0};
  Float_t vPyRcRefit[nRcMax] = {0};
  Float_t vPzRcRefit[nRcMax] = {0};
  Int_t vTofFlag[nRcMax] = {0};
  Float_t vMass2[nRcMax] = {0};
  Float_t vDedx[nRcMax] = {0};
  Int_t vMcTrackRef[nRcMax] = {0};
  Short_t vNhits[nRcMax] = {0};
  Bool_t vEdgeCut[nRcMax] = {0};
  Char_t vCharge[nRcMax] = {0};
  Float_t vTofWeight[nRcMax] = {0};
  Bool_t vTofTrueMatch[nRcMax] = {0};
  Float_t vPairX[nPairMax] = {0};
  Float_t vPairY[nPairMax] = {0};
  Float_t vPairZ[nPairMax] = {0};
  Float_t vPairPx[nPairMax] = {0};
  Float_t vPairPy[nPairMax] = {0};
  Float_t vPairPz[nPairMax] = {0};
  Float_t vPairMass[nPairMax] = {0};
  Float_t vPairChi2[nPairMax] = {0};
  Float_t vPairOpening[nPairMax] = {0};
  Int_t vPairIndex1[nPairMax] = {0};
  Int_t vPairIndex2[nPairMax] = {0};
  Float_t phiEpMc = 0;

  Float_t vParR[nRcMax] = {0};
  Float_t vParZ[nRcMax] = {0};
  Float_t vParPhi[nRcMax] = {0};
  Float_t vParLam[nRcMax] = {0};
  Float_t vParQPt[nRcMax] = {0};

  Float_t vChi2KF[nRcMax] = {0};
  Float_t vCovR[nRcMax] = {0};
  Float_t vCovZ[nRcMax] = {0};
  Float_t vCovPhi[nRcMax] = {0};
  Float_t vCovLam[nRcMax] = {0};
  Float_t vCovQPt[nRcMax] = {0};

  convTree->Branch("nMc", &nMc);
  convTree->Branch("nRc", &nRc);
  convTree->Branch("nPair", &nPair);
  convTree->Branch("bMc", &bMc);
  convTree->Branch("phiEpMc", &phiEpMc);
  convTree->Branch("xVtxMc", &xVtxMc);
  convTree->Branch("yVtxMc", &yVtxMc);
  convTree->Branch("zVtxMc", &zVtxMc);
  convTree->Branch("xVtxRc", &xVtxRc);
  convTree->Branch("yVtxRc", &yVtxRc);
  convTree->Branch("zVtxRc", &zVtxRc);
  convTree->Branch("vPxMc", &vPxMc, "vPxMc[nMc]/F");
  convTree->Branch("vPyMc", &vPyMc, "vPyMc[nMc]/F");
  convTree->Branch("vPzMc", &vPzMc, "vPzMc[nMc]/F");
  convTree->Branch("vxMc", &vxMc, "vxMc[nMc]/F");
  convTree->Branch("vyMc", &vyMc, "vyMc[nMc]/F");
  convTree->Branch("vzMc", &vzMc, "vzMc[nMc]/F");
  convTree->Branch("vPdgCode", &vPdgCode, "vPdgCode[nMc]/I");
  convTree->Branch("vMotherId", &vMotherId, "vMotherId[nMc]/I");
  convTree->Branch("vMcTrackId", &vMcTrackId, "vMcTrackId[nMc]/I");
  convTree->Branch("vDcaX", &vDcaX, "vDcaX[nRc]/F");
  convTree->Branch("vDcaY", &vDcaY, "vDcaY[nRc]/F");
  convTree->Branch("vDcaZ", &vDcaZ, "vDcaZ[nRc]/F");
  convTree->Branch("vPxRc", &vPxRc, "vPxRc[nRc]/F");
  convTree->Branch("vPyRc", &vPyRc, "vPyRc[nRc]/F");
  convTree->Branch("vPzRc", &vPzRc, "vPzRc[nRc]/F");
  convTree->Branch("vPxRcRefit", &vPxRcRefit, "vPxRcRefit[nRc]/F");
  convTree->Branch("vPyRcRefit", &vPyRcRefit, "vPyRcRefit[nRc]/F");
  convTree->Branch("vPzRcRefit", &vPzRcRefit, "vPzRcRefit[nRc]/F");
  convTree->Branch("vTofFlag", &vTofFlag, "vTofFlag[nRc]/I");
  convTree->Branch("vMass2", &vMass2, "vMass2[nRc]/F");
  convTree->Branch("vDedx", &vDedx, "vDedx[nRc]/F");
  convTree->Branch("vMcTrackRef", &vMcTrackRef, "vMcTrackRef[nRc]/I");
  convTree->Branch("vNhits", &vNhits, "vNhits[nRc]/S");
  convTree->Branch("vEdgeCut", &vEdgeCut, "vEdgeCut[nRc]/O");
  convTree->Branch("vCharge", &vCharge, "vCharge[nRc]/B");
  convTree->Branch("vTofWeight", &vTofWeight, "vTofWeight[nRc]/F");
  convTree->Branch("vTofTrueMatch", &vTofTrueMatch, "vTofTrueMatch[nRc]/O");
  convTree->Branch("vPairX", &vPairX, "vPairX[nPair]/F");
  convTree->Branch("vPairY", &vPairY, "vPairY[nPair]/F");
  convTree->Branch("vPairZ", &vPairZ, "vPairZ[nPair]/F");
  convTree->Branch("vPairPx", &vPairPx, "vPairPx[nPair]/F");
  convTree->Branch("vPairPy", &vPairPy, "vPairPy[nPair]/F");
  convTree->Branch("vPairPz", &vPairPz, "vPairPz[nPair]/F");
  convTree->Branch("vPairMass", &vPairMass, "vPairMass[nPair]/F");
  convTree->Branch("vPairChi2", &vPairChi2, "vPairChi2[nPair]/F");
  convTree->Branch("vPairOpening", &vPairOpening, "vPairOpening[nPair]/F");
  convTree->Branch("vPairIndex1", &vPairIndex1, "vPairIndex1[nPair]/I");
  convTree->Branch("vPairIndex2", &vPairIndex2, "vPairIndex2[nPair]/I");

  convTree->Branch("vChi2KF", &vChi2KF, "vChi2KF[nRc]/F");
  convTree->Branch("vParR",   &vParR,   "vParR[nRc]/F");
  convTree->Branch("vParZ",   &vParZ,   "vParZ[nRc]/F");
  convTree->Branch("vParPhi", &vParPhi, "vParPhi[nRc]/F");
  convTree->Branch("vParLam", &vParLam, "vParLam[nRc]/F");
  convTree->Branch("vParQPt", &vParQPt, "vParQPt[nRc]/F");

  convTree->Branch("vChi2KF", &vChi2KF, "vChi2KF[nRc]/F");
  convTree->Branch("vCovR",   &vCovR,   "vCovR[nRc]/F");
  convTree->Branch("vCovZ",   &vCovZ,   "vCovZ[nRc]/F");
  convTree->Branch("vCovPhi", &vCovPhi, "vCovPhi[nRc]/F");
  convTree->Branch("vCovLam", &vCovLam, "vCovLam[nRc]/F");
  convTree->Branch("vCovQPt", &vCovQPt, "vCovQPt[nRc]/F");

  //--------------------------------------------------

  // useful vectors
  TVector3 vRc, pRc, pRcRefit;
  TClonesArray *vPart = new TClonesArray("MpdParticle");
  vector<MpdParticle *> vpair;
  TLorentzVector p4pair;

  // cuts
  Int_t minTPChits = 10;
  Float_t maxEta = 1.5;
  Float_t minEta = -1.5;
  Float_t ptMax = 2.5;
  Float_t mEl = 0.000511;

  int tofMatchingIndices[nRcMax];

  // Loop over events
  Int_t nEvents = inTree->GetEntries();

  for (Int_t ev = 0; ev < nEvents; ev++)
  {
    inTree->GetEntry(ev);
    phiEpMc = MCHeader->GetRotZ();

    bMc = MCHeader->GetB();
    xVtxMc = MCHeader->GetX();
    yVtxMc = MCHeader->GetY();
    zVtxMc = MCHeader->GetZ();

    MpdVertex *vertex = (MpdVertex *)vertexes->First();
    TVector3 primaryVertex;
    vertex->Position(primaryVertex);
    xVtxRc = primaryVertex.X();
    yVtxRc = primaryVertex.Y();
    zVtxRc = primaryVertex.Z();

    nRc = 0;

    TClonesArray *MpdGlobalTracks = MPDEvent->GetGlobalTracks();
    Int_t nRcTracks = MpdGlobalTracks->GetEntriesFast();
    int nMcTracks = MCTracks->GetEntriesFast();

    // build indices from reconstructed global tracks to tof matchings
    for (int i = 0; i < nRcTracks; i++)
      tofMatchingIndices[i] = -1;
    for (int i = 0; i < tofMatchings->GetEntriesFast(); i++)
    {
      MpdTofMatchingData *matching = (MpdTofMatchingData *)tofMatchings->UncheckedAt(i);
      tofMatchingIndices[matching->GetKFTrackIndex()] = i;
    }

    TBits isRcStored;

    int n_wrong = 0;
    for (Int_t i = 0; i < nRcTracks; i++)
    {
      MpdTrack *rcTrack = (MpdTrack *)MpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *kfTrack = (MpdTpcKalmanTrack *)mpdKalmanTracks->UncheckedAt(i);
      if (kfTrack->GetNofWrong() > 0)
      {
        n_wrong++;
        continue;
      }
      Int_t nHits = rcTrack->GetNofHits();
      Float_t eta = rcTrack->GetEta();
      if (nHits < minTPChits)
        continue;
      if (!(eta > minEta && eta < maxEta))
        continue;
      Int_t mcTrackId = rcTrack->GetID();
      Float_t pt = TMath::Abs(rcTrack->GetPt());
      Float_t theta = rcTrack->GetTheta();
      Bool_t edgeCut = rcTrack->GetEdgeCut();
      Char_t charge = rcTrack->GetCharge() > 0 ? 1 : -1;
      Int_t tofFlag = rcTrack->GetTofFlag();
      Float_t mass2 = rcTrack->GetTofMass2();
      Float_t dedx = rcTrack->GetdEdXTPC();
      Float_t p = theta > 1e-10 ? pt / sin(theta) : 0;
      Float_t dcax = rcTrack->GetDCAX();
      Float_t dcay = rcTrack->GetDCAY();
      Float_t dcaz = rcTrack->GetDCAZ();
      if (edgeCut)
        continue;
      if (!(pt > 0.05))
        continue;

      // check TOF matching info
      float tofWeight = -1;
      bool tofTrueMatch = 0;

      if (tofMatchingIndices[i] >= 0)
      {
        MpdTofMatchingData *matching = (MpdTofMatchingData *)tofMatchings->UncheckedAt(tofMatchingIndices[i]);
        tofWeight = matching->GetWeight();
        MpdTofHit *tofHit = (MpdTofHit *)tofHits->UncheckedAt(matching->GetTofHitIndex());
        vector<int> links;
        tofHit->getLinks(MpdTofUtils::mcTrackIndex, links);
        tofTrueMatch = 0;
        for (unsigned int iMcTrackInHit = 0; iMcTrackInHit < links.size(); iMcTrackInHit++)
        {
          if (links[iMcTrackInHit] == mcTrackId)
            tofTrueMatch = 1;
        }
      }
      float meanTPC = fMeanTPC.Eval(p);
      float sigmTPC = fSigmTPC.Eval(p);
      float meanTOF = fMeanTOF.Eval(p);
      float sigmTOF = fSigmTOF.Eval(p);
      vTofFlag[nRc] = tofFlag;
      vMass2[nRc] = mass2;
      auto* mcTrack = (MpdMCTrack *)MCTracks->UncheckedAt(mcTrackId);
      double mc_mass = mcTrack->GetMass();
      vDedx[nRc] = dedx;
      MpdTpcKalmanTrack kfTrackRefit = *kfTrack;
      kfTrackRefit.SetDirection(MpdKalmanTrack::kInward);
      bool isRefitOk = recoTpc->Refit(&kfTrackRefit, mc_mass, charge, kFALSE, 1);
      if (!isRefitOk)
        continue;
      pRcRefit = kfTrackRefit.Momentum3();
      pRc.SetPtThetaPhi(pt, theta, rcTrack->GetPhi());
      vPxRc[nRc] = pRc[0];
      vPyRc[nRc] = pRc[1];
      vPzRc[nRc] = pRc[2];
      vPxRcRefit[nRc] = pRcRefit[0];
      vPyRcRefit[nRc] = pRcRefit[1];
      vPzRcRefit[nRc] = pRcRefit[2];
      vDcaX[nRc] = dcax;
      vDcaY[nRc] = dcay;
      vDcaZ[nRc] = dcaz;
      vMcTrackRef[nRc] = mcTrackId;
      vNhits[nRc] = nHits;
      vEdgeCut[nRc] = edgeCut;
      vCharge[nRc] = charge;
      vTofWeight[nRc] = tofWeight;
      vTofTrueMatch[nRc] = tofTrueMatch;
      isRcStored.SetBitNumber(mcTrackId);
      // get fit parameters
      vParR[nRc] = kfTrackRefit.GetParam(0);
      vParZ[nRc] = kfTrackRefit.GetParam(1);
      vParPhi[nRc] = kfTrackRefit.GetParam(2);
      vParLam[nRc] = kfTrackRefit.GetParam(3);
      vParQPt[nRc] = kfTrackRefit.GetParam(4);
      vChi2KF[nRc] = kfTrackRefit.GetChi2();
      vCovR[nRc] = (*kfTrackRefit.GetCovariance())(0, 0);
      vCovZ[nRc] = (*kfTrackRefit.GetCovariance())(1, 1);
      vCovPhi[nRc] = (*kfTrackRefit.GetCovariance())(2, 2);
      vCovLam[nRc] = (*kfTrackRefit.GetCovariance())(3, 3);
      vCovQPt[nRc] = (*kfTrackRefit.GetCovariance())(4, 4);
      nRc++;
      if (nRc > nRcMax)
      {
        printf("warning: increase nRcMax\n");
        return;
      }
    }

    // store mothers and grandmothers of stored reconstructed particles
    for (Int_t i = MCTracks->GetEntriesFast() - 1; i >= 0; i--)
    {
      if (isRcStored.TestBitNumber(i) == 0)
        continue;
      MpdMCTrack *mcTrack = (MpdMCTrack *)MCTracks->UncheckedAt(i);
      Int_t iMother = mcTrack->GetMotherId();
      isRcStored.SetBitNumber(iMother);
    }

    nMc = 0;
    for (Int_t i = 0; i < MCTracks->GetEntriesFast(); i++)
    {
      MpdMCTrack *mcTrack = (MpdMCTrack *)MCTracks->UncheckedAt(i);
      vxMc[nMc] = mcTrack->GetStartX();
      vyMc[nMc] = mcTrack->GetStartY();
      vzMc[nMc] = mcTrack->GetStartZ();
      if (!isRcStored.TestBitNumber(i)  && vxMc[nMc] * vxMc[nMc] + vyMc[nMc] * vyMc[nMc] > 40. * 40.)
        continue;
      Int_t pdgCode = mcTrack->GetPdgCode();
      Int_t motherId = mcTrack->GetMotherId();
      if (!isRcStored.TestBitNumber(i) && !(pdgCode == 11 || pdgCode == -11 || pdgCode == 22 ||
                                            pdgCode == 221 || pdgCode == -211 || pdgCode == 211 ||
                                            pdgCode == 2212 || pdgCode == -2212) )
      {
        continue;
      }
      vPxMc[nMc] = mcTrack->GetPx();
      vPyMc[nMc] = mcTrack->GetPy();
      vPzMc[nMc] = mcTrack->GetPz();
      vPdgCode[nMc] = mcTrack->GetPdgCode();
      vMotherId[nMc] = mcTrack->GetMotherId();
      vMcTrackId[nMc] = i;
      nMc++;
      if (nMc > nMcMax)
      {
        printf("warning: increase nMcMax\n");
        return;
      }
    }
    printf("Event: %i", ev);
    printf(" nMc=%i ", nMc);
    printf(" nRc=%i ", nRc);
    printf(" n_wrong=%i ", n_wrong);
    printf("\n");
    convTree->Fill();
  }
  convTree->Write();
  convFile->Close();
}
