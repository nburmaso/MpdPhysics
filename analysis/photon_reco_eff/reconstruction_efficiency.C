#include "Rtypes.h"
#include "TSystem.h"
#include <TVector3.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TStyle.h"
#include "tree.h"
#include "map"
#include "TFile.h"
#include "TF1.h"
using namespace std;

Int_t nHitsMin=20;
Float_t ptMin=0.05;
Float_t etaMax=1.0;
Float_t chi2Max=10;
Float_t psiMax=10;
Float_t rMin=10;
Float_t yPiMax=1.;

Float_t p0m = 0.014327;
Float_t p1m = 0.023507;
Float_t p0p = -0.49387;
Float_t p1p = -9.19766;
Float_t p2p =  0.03194;


Bool_t IsElectronSelected(TVector3 pe, Int_t nHits, Float_t dedx, Int_t tofFlag, Float_t mass2){
  if (nHits<nHitsMin) return 0;
  if (pe.Pt()<ptMin) return 0;
  if (fabs(pe.Eta())>etaMax) return 0;
  double p = pe.Mag();
  float meanTPC = fMeanTPC.Eval(p);
  float sigmTPC = fSigmTPC.Eval(p);
  float meanTOF = fMeanTOF.Eval(p);
  float sigmTOF = fSigmTOF.Eval(p);
  if (dedx-meanTPC> 4*sigmTPC) return 0;
  if (dedx-meanTPC<-2*sigmTPC) return 0;
  if ((tofFlag==2 || tofFlag==6) && mass2-meanTOF>4*sigmTOF) return 0;
  return 1;
}

Bool_t IsElectronPairSelected(Float_t m, Float_t chi2, Float_t r, Float_t openingAngle, Float_t pointingAngle, Float_t psi, Float_t pt){
  if (isnan(m)) return 0;
  if (pointingAngle>exp(p0p+p1p*pt)+p2p) return 0;
  if (psi>psiMax) return 0;
  if (m>p0m+p1m*pt) return 0;
  if (chi2>chi2Max) return 0;
  if (chi2<0) return 0;
  if (r<rMin) return 0;
  if (pt<0.1) return 0;
  return 1;
}

void reconstruction_efficiency() {
  TChain *tree = new TChain("convTree");
  tree->Add("conv.root");
  SetBranchAddresses(tree);
  int nEvents = tree->GetEntries();
  printf("nEvents = %i\n",nEvents);
  TH1D* h = new TH1D("h","",200,-100,100);
  TH1D* hPtAll  = new TH1D("hPtAll","",60,0,1.5);
  TH1D* hPtConv = new TH1D("hPtConv","",60,0,1.5);
  TH1D* hPtAcc  = new TH1D("hPtAcc","",60,0,1.5);
  TH1D* hPtRec  = new TH1D("hPtRec","",60,0,1.5);
  TH1D* hPtRecGa  = new TH1D("hPtRecGa","",60,0,1.5);
  TH1D* hEtaAll  = new TH1D("hEtaAll","",20,0.,1.);
  TH1D* hEtaConv = new TH1D("hEtaConv","",20,0.,1.);

  TH2D* hPtEtaAll    = new TH2D("hPtEtaAll","",60,0,1.5,100,-1.,1.);
  TH2D* hPtEtaConv   = new TH2D("hPtEtaConv","",60,0,1.5,100,-1.,1.);
  TH2D* hPtEtaAcc    = new TH2D("hPtEtaAcc","",60,0,1.5,100,-1.,1.);
  TH2D* hPtEtaRec    = new TH2D("hPtEtaRec","",60,0,1.5,100,-1.,1.);
  TH2D* hPtEtaRecGa  = new TH2D("hPtEtaRecGa","",60,0,1.5,100,-1.,1.);
  map<int,int> mcTrackId_to_iMc;

  fMeanTPC.SetParameters(-9.351e-03,4.328e+03,-3.895e-04);
  fSigmTPC.SetParameters(338.877,-502.381,1572.066,-1857.843,917.432,-161.531);
  fMeanTOF.SetParameters(0);
  fSigmTOF.SetParameters(7.60503e-04,-4.54220e-04,2.55733e-02);

  for (int iEv = 0; iEv < nEvents; iEv++) {
    if (iEv%100000==0) printf("Event: %i\n",iEv);
    tree->GetEntry(iEv);
    int nTracks = nMc;
    int nTracksPassed = 0;
    int nElectrons = 0;
    int nElectronsAccepted = 0;
    double pt = -1;
    double eta = -2;
    TLorentzVector v;
    for (int i=0;i<nTracks;i++){
      double x = vxMc[i];
      double y = vyMc[i];
      double z = vzMc[i];
      v.SetXYZT(vPxMc[i],vPyMc[i],vPzMc[i],0);
      double p = v.P();
      int moId = vMotherId[i];
      int pdgCode = vPdgCode[i];

      if (x*x+y*y>40*40) continue;
      if (fabs(z)>150) continue;
      if (moId>=1) continue; // analysing dedicated gamma production with gamma moId = 0
      if (p<0.01) continue;
      h->Fill(pdgCode);
      if (fabs(pdgCode)==11 && moId==0) nElectrons++;
      if (fabs(pdgCode)==22) pt = v.Pt();
      if (fabs(pdgCode)==22) eta = v.Eta();
      if (fabs(pdgCode)==11 && moId==0 && v.Pt()>0.05 && fabs(v.Eta())<1) nElectronsAccepted++;
      nTracksPassed++;
    }

    hPtAll->Fill(pt);
    if (nElectrons==2) hPtConv->Fill(pt);
    if (nElectronsAccepted==2) hPtAcc->Fill(pt);

    hPtEtaAll->Fill(pt,eta);
    if (nElectrons==2) hPtEtaConv->Fill(pt,eta);
    if (nElectronsAccepted==2) hPtEtaAcc->Fill(pt,eta);


    // accepted electrons
    mcTrackId_to_iMc.clear();
    for (Int_t i=0;i<nMc;i++) mcTrackId_to_iMc[vMcTrackId[i]]=i;

    int nElectronsReconstructed = 0;
    for (Int_t iRc=0;iRc<nRc;iRc++){
      Int_t iMc   = mcTrackId_to_iMc[vMcTrackRef[iRc]];
      Int_t pdg   = vPdgCode[iMc];
      Int_t moId  = vMotherId[iMc];
      if (moId!=0 && abs(pdg)!=11) continue;
      TVector3 p1e(vPxRcRefit[iRc],vPyRcRefit[iRc],vPzRcRefit[iRc]);
      if (vNhits[iRc]<nHitsMin) continue;
      if (p1e.Pt()<ptMin) continue;
      if (fabs(p1e.Eta())>etaMax) continue;
      nElectronsReconstructed++;
    }

    if (nElectronsReconstructed<2) continue;
    TVector3 vtx0(xVtxMc,yVtxMc,zVtxMc);
    Bool_t isPairAccepted = 0;
    for (Int_t iPair=0;iPair<nPair;iPair++){
      Int_t i1 = vPairIndex1[iPair];
      Int_t i2 = vPairIndex2[iPair];
      if (i1==i2) continue;
      TVector3 p1e(vPxRcRefit[i1],vPyRcRefit[i1],vPzRcRefit[i1]);
      TVector3 p2e(vPxRcRefit[i2],vPyRcRefit[i2],vPzRcRefit[i2]);
      if (!IsElectronSelected(p1e,vNhits[i1],vDedx[i1],vTofFlag[i1],vMass2[i1])) continue;
      if (!IsElectronSelected(p2e,vNhits[i2],vDedx[i2],vTofFlag[i2],vMass2[i2])) continue;
      TVector3 p12(vPairPx[iPair],vPairPy[iPair],vPairPz[iPair]);
      if (isnan(vPairMass[iPair])) continue;
      if (p12.Pt()<0.1) continue;
      hPtRec->Fill(pt);
      hPtEtaRec->Fill(pt,eta);
      TVector3 vtx(vPairX[iPair],vPairY[iPair],vPairZ[iPair]);
      Float_t openingAngle = vPairOpening[iPair];
      Float_t openingAngleInPolarDirection = p1e.Theta()-p2e.Theta();
      Float_t psi = fabs(TMath::ASin(openingAngleInPolarDirection/openingAngle));
      Float_t r = vtx.Perp();
      vtx-=vtx0;
      Float_t pointingAngle = vtx.Angle(p12);
      if (!IsElectronPairSelected(vPairMass[iPair],vPairChi2[iPair],r, openingAngle, pointingAngle, psi, p12.Pt())) continue;
      isPairAccepted=1;
    }
    if (isPairAccepted==0) continue;
    hPtRecGa->Fill(pt);
    hPtEtaRecGa->Fill(pt,eta);
  }

  TFile* f = new TFile("reconstruction_efficiency.root","recreate");
  h->Write();
  hPtConv->Write();
  hPtAll->Write();
  hPtAcc->Write();
  hEtaConv->Write();
  hEtaAll->Write();
  hPtRec->Write();
  hPtRecGa->Write();
  hPtEtaAll->Write();
  hPtEtaConv->Write();
  hPtEtaAcc->Write();
  hPtEtaRec->Write();
  hPtEtaRecGa->Write();
  f->Close();
}
