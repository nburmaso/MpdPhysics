#include "TClonesArray.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TFile.h"
#include "tree.h"
using namespace std;
using namespace TMath;
Float_t pi=Pi();

Int_t nHitsMin=20;
Float_t ptMin=0.05;
Float_t etaMax=1.0;
Float_t chi2Max=10;
Float_t psiMax= 10;;
Float_t rMin=5.;
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

Bool_t IsPhotonPairSelected(TLorentzVector p){
  Double_t y = p.Rapidity();
  if (fabs(y)>yPiMax) return 0;
  return 1;
}

map<int,int> mcTrackId_to_iMc;
void GetMcInfo(Int_t iRc, Int_t &mcTrackId, Int_t &iMc, Int_t &pdg, Int_t &moId, Int_t &iMoMc, Int_t &moPdg, Int_t &grId, Int_t &iGrMc, Int_t &grPdg){
  iMc   = -2;
  pdg   = -2;
  moId  = -2;
  iMoMc = -2;
  moPdg = -2;
  iGrMc = -2;
  grId  = -2;
  grPdg = -2;
  mcTrackId   = vMcTrackRef[iRc];
  if (mcTrackId_to_iMc.count(mcTrackId)==0) return;
  iMc  = mcTrackId_to_iMc[mcTrackId];
  pdg  = vPdgCode[iMc];
  moId = vMotherId[iMc];
  if (moId<0) return;
  if (mcTrackId_to_iMc.count(moId)==0) return;
  iMoMc = mcTrackId_to_iMc[moId];
  moPdg = vPdgCode[iMoMc];
  grId  = vMotherId[iMoMc];
  if (grId<0) return;
  if (mcTrackId_to_iMc.count(grId)==0) return;
  iGrMc = mcTrackId_to_iMc[grId];
  grPdg = vPdgCode[iGrMc];
}


void analysis(){
  fMeanTPC.SetParameters(-9.351e-03,4.328e+03,-3.895e-04);
  fSigmTPC.SetParameters(338.877,-502.381,1572.066,-1857.843,917.432,-161.531);
  fMeanTOF.SetParameters(0);
  fSigmTOF.SetParameters(7.60503e-04,-4.54220e-04,2.55733e-02);

  gStyle->SetPalette(1);
  TChain* convTree = new TChain("convTree");
  // set path to your trees here
  for (Int_t i=300;i<=399;i++)
    for (Int_t j=0;j<=27;j++)
      convTree->Add(Form("/d/data/BiBi09/conv.%04i.%02i.root",i,j));
  Int_t nEvents = convTree->GetEntries();
  printf("nEvents=%i\n",nEvents);
  SetBranchAddresses(convTree);
  TH1D* hB    = new TH1D("hB","",100,0,20);
  TH1D* hCent = new TH1D("hCent","",100,0,100);
  TH1D* hMult = new TH1D("hMult","",100,0,800);
  TH3D* hPi0CentPtY = new TH3D("hPi0CentPtY","",100,0,100,80,0,4,100,-5,5);
  TH3D* hEtaCentPtY = new TH3D("hEtaCentPtY","",100,0,100,80,0,4,100,-5,5);
  TH3D* hGenCentPtM = new TH3D("hGenCentPtM","",100,0,100,80,0,4,400,0,1);
  TH3D* hRecCentPtM = new TH3D("hRecCentPtM","",100,0,100,80,0,4,400,0,1);
  TH3D* hSigCentPtM = new TH3D("hSigCentPtM","",100,0,100,80,0,4,400,0,1);
  TH3D* hGenMultPtM = new TH3D("hGenMultPtM","",100,0,800,80,0,4,400,0,1);
  TH3D* hRecMultPtM = new TH3D("hRecMultPtM","",100,0,800,80,0,4,400,0,1);
  TH3D* hSigMultPtM = new TH3D("hSigMultPtM","",100,0,800,80,0,4,400,0,1);

  TH2D* hGaAllCentPt = new TH2D("hGaAllCentPt","",100,0,100,80,0,4);
  TH2D* hGaPi0CentPt = new TH2D("hGaPi0CentPt","",100,0,100,80,0,4);
  TH2D* hGaEtaCentPt = new TH2D("hGaEtaCentPt","",100,0,100,80,0,4);
  TH2D* hGaOthCentPt = new TH2D("hGaOthCentPt","",100,0,100,80,0,4);
  TH2D* hGaBreCentPt = new TH2D("hGaBreCentPt","",100,0,100,80,0,4); // gamma from bremstrahlung of electrons
  TH2D* hGaSigCentPt = new TH2D("hGaSigCentPt","",100,0,100,80,0,4); // gamma from sigma0
  TH2D* hGaBgdCentPt = new TH2D("hGaBgdCentPt","",100,0,100,80,0,4); // random background


  TH2D* hDedxAll   = new TH2D("hDedxAll","",60,0,1.5,100,0,10000);
  TH2D* hDedxSel   = new TH2D("hDedxSel","",60,0,1.5,100,0,10000);

  Int_t mcTrackId, iMc, vPdg[nRcMax], vMoId[nRcMax], iMoMc, vMoPdg[nRcMax], vGrId[nRcMax], iGrMc, vGrPdg[nRcMax];
  TLorentzVector vg[nPairMax];
  Int_t v1[nPairMax];
  Int_t v2[nPairMax];
  bool  vIsRcSelected[nRcMax];
  float vRcTheta[nRcMax];

  for (Int_t ev=0;ev<nEvents;ev++){
    if (ev%1000==0) printf("Event: %i\n",ev);
    convTree->GetEntry(ev);
    TVector3 vtx0(xVtxRc,yVtxRc,zVtxRc);
    if (fabs(zVtxRc)>50) continue;
    if (nMc<=2) continue;
    hB->Fill(bMc);
    double cent = 100*bMc*bMc/15./15.;
    TLorentzVector pMc;
    double mult = 0;
    Int_t nCharged = 0;
    for (Int_t i=0;i<nMc;i++){
      pMc.SetXYZM(vPxMc[i],vPyMc[i],vPzMc[i],0);
      if (fabs(pMc.Eta())<5 && (abs(vPdgCode[i])==2212 || abs(vPdgCode[i])==211 || abs(vPdgCode[i])==11 || abs(vPdgCode[i])==321)) nCharged++;
      if (fabs(vxMc[i])>1) continue;
      if (fabs(vyMc[i])>1) continue;
      if (fabs(vzMc[i])>50) continue;
      if (vPdgCode[i]!=22) continue;
      if (fabs(pMc.Eta())>1) continue;
      mult++;
    }
    if (nCharged>0) hCent->Fill(cent);
    hMult->Fill(mult);

    mcTrackId_to_iMc.clear();
    for (Int_t i=0;i<nMc;i++){
      mcTrackId_to_iMc[vMcTrackId[i]]=i;
      Int_t pdg = vPdgCode[i];
      if (pdg!=111 && pdg!=221) continue;
      pMc.SetXYZM(vPxMc[i],vPyMc[i],vPzMc[i],pdg==111 ? mPi0 : mEta);
      double m  = pMc.M();
      double pt = pMc.Pt();
      double y  = pMc.Rapidity();
      if (pdg==111 && fabs(vxMc[i])<1 && fabs(vyMc[i])<1) hPi0CentPtY->Fill(cent,pt,y);
      if (pdg==221 && fabs(vxMc[i])<1 && fabs(vyMc[i])<1) hEtaCentPtY->Fill(cent,pt,y);
      if (fabs(y)>yPiMax) continue;
      hGenCentPtM->Fill(cent,pt,m);
      hGenMultPtM->Fill(mult,pt,m);
    }

    for (Int_t i=0;i<nRc;i++){
      GetMcInfo(i, mcTrackId, iMc, vPdg[i], vMoId[i], iMoMc, vMoPdg[i], vGrId[i], iGrMc, vGrPdg[i]);
      if (vPdg[i]==-2) printf("%i\n",mcTrackId);
    }

    for (Int_t iRc=0;iRc<nRc;iRc++){
      TVector3 pe(vPxRcRefit[iRc],vPyRcRefit[iRc],vPzRcRefit[iRc]);
      // if (vTofTrueMatch[iRc]==0) vTofFlag[iRc]=0;
      vIsRcSelected[iRc] = IsElectronSelected(pe,vNhits[iRc],vDedx[iRc],vTofFlag[iRc],vMass2[iRc]);
      vRcTheta[iRc] = pe.Theta();
      if (vNhits[iRc]<nHitsMin) continue;
      if (pe.Pt()<ptMin) continue;
      if (fabs(pe.Eta())>etaMax) continue;
      if (fabs(vPdg[iRc])!=11 || vMoPdg[iRc]!=22) continue;
      hDedxAll->Fill(pe.Mag(),vDedx[iRc]);
      if (!vIsRcSelected[iRc]) continue;
      hDedxSel->Fill(pe.Mag(),vDedx[iRc]);
    }

    Int_t ng=0;
    // select photons
    for (Int_t iPair=0;iPair<nPair;iPair++){
      Int_t i1 = vPairIndex1[iPair];
      Int_t i2 = vPairIndex2[iPair];
      if (i1==i2) continue;
      if (!vIsRcSelected[i1]) continue;
      if (!vIsRcSelected[i2]) continue;
      TVector3 vtx(vPairX[iPair],vPairY[iPair],vPairZ[iPair]);
      TVector3 p12(vPairPx[iPair],vPairPy[iPair],vPairPz[iPair]);
      Float_t openingAngle = vPairOpening[iPair];
      Float_t openingAngleInPolarDirection = vRcTheta[i1] - vRcTheta[i2];
      Float_t psi = fabs(TMath::ASin(openingAngleInPolarDirection/openingAngle));
      Float_t r = vtx.Perp();
      vtx-=vtx0;
      Float_t pointingAngle = vtx.Angle(p12);
      if (!IsElectronPairSelected(vPairMass[iPair],vPairChi2[iPair],r, openingAngle, pointingAngle, psi, p12.Pt())) continue;
      vg[ng].SetVectM(p12,0);
      v1[ng]=i1;
      v2[ng]=i2;
      ng++;
    }

    // select pi0
    for (Int_t i=0;i<ng;i++){
      for (Int_t j=i+1;j<ng;j++){
        TLorentzVector p = vg[i]+vg[j];
        Double_t m  = p.M();
        Double_t pt = p.Pt();
        if (!IsPhotonPairSelected(p)) continue;
        if (v1[i]==v1[j] || v2[i]==v2[j]) continue;
        hRecCentPtM->Fill(cent,pt,m);
        hRecMultPtM->Fill(mult,pt,m);
        if (fabs(vPdg[v1[i]])!=11) continue;
        if (fabs(vPdg[v2[i]])!=11) continue;
        if (fabs(vPdg[v1[j]])!=11) continue;
        if (fabs(vPdg[v2[j]])!=11) continue;
        if (vMoId[v1[i]]!=vMoId[v2[i]]) continue;
        if (vMoId[v1[j]]!=vMoId[v2[j]]) continue;
        if (vMoPdg[v1[i]]!=22) continue;
        if (vMoPdg[v1[j]]!=22) continue;
        if (vGrId[v1[i]]!=vGrId[v1[j]]) continue;
        if (vGrPdg[v1[i]]!=111 && vGrPdg[v1[i]]!=221) continue;
        hSigCentPtM->Fill(cent,pt,m);
        hSigMultPtM->Fill(mult,pt,m);
      }
    }

    // select photons from pi0 and eta
    for (Int_t i=0;i<ng;i++){
      Double_t pt = vg[i].Pt();
      hGaAllCentPt->Fill(cent,pt);
      if (vPdg[v1[i]]==11 && vPdg[v2[i]]==-11 && vMoId[v1[i]]==vMoId[v2[i]] && vMoPdg[v1[i]]==22) {
        if      (vGrPdg[v1[i]]==111)       hGaPi0CentPt->Fill(cent,pt);
        else if (vGrPdg[v1[i]]==221)       hGaEtaCentPt->Fill(cent,pt);
        else if (abs(vGrPdg[v1[i]])==11)   hGaBreCentPt->Fill(cent,pt);
        else if (abs(vGrPdg[v1[i]])==3212) hGaSigCentPt->Fill(cent,pt);
        else {
          hGaOthCentPt->Fill(cent,pt);
        }
      } else {
        hGaBgdCentPt->Fill(cent,pt);
        // printf("%i %i %i %i %i %i %i %i\n",vPdg[v1[i]],vPdg[v2[i]],vMoId[v1[i]],vMoId[v2[i]],vMoPdg[v1[i]],vMoPdg[v2[i]],vGrId[v1[i]],vGrId[v2[i]]);
      }
    }
  }

  TFile* fout = new TFile("results.BiBi09.root","recreate");
  hB->Write();
  hCent->Write();
  hMult->Write();
  hGenCentPtM->Write();
  hRecCentPtM->Write();
  hSigCentPtM->Write();
  hGenMultPtM->Write();
  hRecMultPtM->Write();
  hSigMultPtM->Write();
  hPi0CentPtY->Write();
  hEtaCentPtY->Write();
  hGaAllCentPt->Write();
  hGaPi0CentPt->Write();
  hGaEtaCentPt->Write();
  hGaBreCentPt->Write();
  hGaSigCentPt->Write();
  hGaOthCentPt->Write();
  hGaBgdCentPt->Write();

  hDedxAll->Write();
  hDedxSel->Write();
  fout->Close();
}
