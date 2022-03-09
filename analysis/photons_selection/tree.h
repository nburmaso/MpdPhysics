#ifndef tree_h
#define tree_h 1
#include "TF1.h"
#include "TTree.h"

Double_t mPi0 = 0.1349766;
Double_t mEta = 0.547862;
Double_t mEl  = 0.000511;

const Int_t nMcMax   = 10000;
const Int_t nRcMax   = 10000;
const Int_t nPairMax = 100000;
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
//Float_t vSigmaX[nRcMax]={0};
//Float_t vSigmaY[nRcMax]={0};
//Float_t vSigmaZ[nRcMax]={0};
Float_t vPxRc[nRcMax]={0};
Float_t vPyRc[nRcMax]={0};
Float_t vPzRc[nRcMax]={0};
Float_t vPxRcRefit[nRcMax]={0};
Float_t vPyRcRefit[nRcMax]={0};
Float_t vPzRcRefit[nRcMax]={0};
Int_t   vTofFlag[nRcMax]={0};
Float_t vMass2[nRcMax]={0};
Float_t vTofWeight[nRcMax]={0};
Bool_t  vTofTrueMatch[nRcMax]={0};
Float_t vDedx[nRcMax]={0};
//Float_t vProbEl[nRcMax]={0};
//Float_t vProbPi[nRcMax]={0};
//Float_t vProbKa[nRcMax]={0};
//Float_t vProbPr[nRcMax]={0};
Int_t   vMcTrackRef[nRcMax]={0};
Short_t vNhits[nRcMax]={0};
Char_t  vCharge[nRcMax]={0};
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
//Float_t vZdcEnergy[90];
Float_t phiEpMc=0;

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

void SetBranchAddresses(TTree* tree){
  tree->SetBranchAddress("phiEpMc",&phiEpMc);
//  tree->SetBranchAddress("vZdcEnergy",&vZdcEnergy);
  tree->SetBranchAddress("nMc",&nMc);
  tree->SetBranchAddress("nRc",&nRc);
  tree->SetBranchAddress("nPair",&nPair);
//  tree->SetBranchAddress("nPrimary",&nPrimary);
  tree->SetBranchAddress("bMc",&bMc);
  tree->SetBranchAddress("xVtxMc",&xVtxMc);
  tree->SetBranchAddress("yVtxMc",&yVtxMc);
  tree->SetBranchAddress("zVtxMc",&zVtxMc);
  tree->SetBranchAddress("vxMc",&vxMc);
  tree->SetBranchAddress("vyMc",&vyMc);
  tree->SetBranchAddress("vzMc",&vzMc);
  tree->SetBranchAddress("xVtxRc",&xVtxRc);
  tree->SetBranchAddress("yVtxRc",&yVtxRc);
  tree->SetBranchAddress("zVtxRc",&zVtxRc);
  tree->SetBranchAddress("vPxMc",&vPxMc);
  tree->SetBranchAddress("vPyMc",&vPyMc);
  tree->SetBranchAddress("vPzMc",&vPzMc);
  tree->SetBranchAddress("vPdgCode",&vPdgCode);
  tree->SetBranchAddress("vMotherId",&vMotherId);
  tree->SetBranchAddress("vMcTrackId",&vMcTrackId);
  tree->SetBranchAddress("vDcaX",&vDcaX);
  tree->SetBranchAddress("vDcaY",&vDcaY);
  tree->SetBranchAddress("vDcaZ",&vDcaZ);
//  tree->SetBranchAddress("vSigmaX",&vSigmaX);
//  tree->SetBranchAddress("vSigmaY",&vSigmaY);
//  tree->SetBranchAddress("vSigmaZ",&vSigmaZ);
  tree->SetBranchAddress("vPxRc",&vPxRc);
  tree->SetBranchAddress("vPyRc",&vPyRc);
  tree->SetBranchAddress("vPzRc",&vPzRc);
  tree->SetBranchAddress("vPxRcRefit",&vPxRcRefit);
  tree->SetBranchAddress("vPyRcRefit",&vPyRcRefit);
  tree->SetBranchAddress("vPzRcRefit",&vPzRcRefit);
  tree->SetBranchAddress("vTofFlag",&vTofFlag);
  tree->SetBranchAddress("vMass2",&vMass2);
  tree->SetBranchAddress("vTofWeight",&vTofWeight);
  tree->SetBranchAddress("vTofTrueMatch",&vTofTrueMatch);
  tree->SetBranchAddress("vDedx",&vDedx);
//  tree->SetBranchAddress("vProbEl",&vProbEl);
//  tree->SetBranchAddress("vProbPi",&vProbPi);
//  tree->SetBranchAddress("vProbKa",&vProbKa);
//  tree->SetBranchAddress("vProbPr",&vProbPr);
  tree->SetBranchAddress("vMcTrackRef",&vMcTrackRef);
  tree->SetBranchAddress("vNhits",&vNhits);
  tree->SetBranchAddress("vCharge",&vCharge);
  tree->SetBranchAddress("vPairX",&vPairX);
  tree->SetBranchAddress("vPairY",&vPairY);
  tree->SetBranchAddress("vPairZ",&vPairZ);
  tree->SetBranchAddress("vPairPx",&vPairPx);
  tree->SetBranchAddress("vPairPy",&vPairPy);
  tree->SetBranchAddress("vPairPz",&vPairPz);
  tree->SetBranchAddress("vPairMass",&vPairMass);
  tree->SetBranchAddress("vPairChi2",&vPairChi2);
  tree->SetBranchAddress("vPairOpening",&vPairOpening);
  tree->SetBranchAddress("vPairIndex1",&vPairIndex1);
  tree->SetBranchAddress("vPairIndex2",&vPairIndex2);

//  tree->SetBranchAddress("vParR",&vCovR);
//  tree->SetBranchAddress("vParZ",&vCovZ);
//  tree->SetBranchAddress("vParPhi",&vCovPhi);
//  tree->SetBranchAddress("vParLam",&vCovLam);
//  tree->SetBranchAddress("vParQPt",&vCovQPt);
//  tree->SetBranchAddress("vChi2KF",&vChi2KF);
//  tree->SetBranchAddress("vCovR",&vCovR);
//  tree->SetBranchAddress("vCovZ",&vCovZ);
//  tree->SetBranchAddress("vCovPhi",&vCovPhi);
//  tree->SetBranchAddress("vCovLam",&vCovLam);
//  tree->SetBranchAddress("vCovQPt",&vCovQPt);
}

TF1 fMeanTPC("fMeanTPC","1./([0]*x+[2])+[1]",0.,5);
TF1 fSigmTPC("fSigmTPC","pol5",0.,5);
TF1 fMeanTOF("fMeanTOF","pol0",0.,5);
TF1 fSigmTOF("fSigmTOF","pol2",0.,5);
TF1 fPidTPC("fPidTPC","x<[1] ? [0]*TMath::Gaus(x,[1],[2]) : [0]*TMath::Gaus(x,[1],[2]*[3])",0,10000);
TF1 fPidTOF("fPidTOF","[0]*TMath::Gaus(x,[1],[2])",-1,1);


float ElectronProbabilityTPC(float p, float dedx, bool edge = 0){
  float meanTPC = fMeanTPC.Eval(p);
  float sigmTPC = fSigmTPC.Eval(p);
  fPidTPC.SetParameters(1.,meanTPC,sigmTPC,1.5);
  if (dedx<0) return 0;
  if (dedx>10000) return 0;
  float total = fPidTPC.Integral(0,10000.);
  if (total<1e-10) return 0;
  if (edge==0) return fPidTPC.Integral(0,dedx)/total;
  else         return fPidTPC.Integral(dedx,10000.)/total;
  return 0;
}

float ElectronProbabilityTOF(float p, float mass2, int tofFlag){
  float meanTOF = fMeanTOF.Eval(p);
  float sigmTOF = fSigmTOF.Eval(p);
  fPidTOF.SetParameters(1.,meanTOF,sigmTOF);
  if (!(tofFlag==2 || tofFlag==6)) return 0.9999;
  if (mass2>1) return 0;
  float total = fPidTOF.Integral(-1.,1.,1e-13);
  float outside = fPidTOF.Integral(mass2,1.,1e-13);
  if (total<1e-20) return 0;
  return outside/total;
}

#endif

