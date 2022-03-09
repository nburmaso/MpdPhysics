#include "TClonesArray.h"
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "map"
#include "TVector3.h"
#include "TF1.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TFile.h"

using namespace std;

void draw_expected_yields(Bool_t isPi0=1){
  gStyle->SetOptStat(0);

  float_t centMin =  0;
  float_t centMax =  100;

  TFile* f = new TFile("results.BiBi09.root");
  TH1D* hCent = (TH1D*) f->Get("hCent");
  TH3D* hGenCentPtM = (TH3D*) f->Get("hGenCentPtM");
  TH3D* hRecCentPtM = (TH3D*) f->Get("hRecCentPtM");
  TH3D* hSigCentPtM = (TH3D*) f->Get("hSigCentPtM");

  Double_t ptMin=0.0;
  Double_t ptMax=2.2;
  Int_t nBins = isPi0 ? 22 : 11;
  Double_t mPi0 = 0.1349766;
  Double_t mEta = 0.547862;

  hGenCentPtM->Rebin3D(1,nBins==22 ? 2 : 4,1);
  hRecCentPtM->Rebin3D(1,nBins==22 ? 2 : 4,1);
  hSigCentPtM->Rebin3D(1,nBins==22 ? 2 : 4,1);

  Int_t iCentMin = hRecCentPtM->GetXaxis()->FindFixBin(centMin+0.001);
  Int_t iCentMax = hRecCentPtM->GetXaxis()->FindFixBin(centMax-0.001);

  TH1D* hPt = new TH1D("hPt","",nBins,ptMin,ptMax);

  for (Int_t ip=1;ip<=nBins;ip++){
    TH1D* hGenM0 = hGenCentPtM->ProjectionZ(Form("hGenM%i",ip),iCentMin,iCentMax,ip,ip);
    TH1D* hRecM0 = hRecCentPtM->ProjectionZ(Form("hRecM%i",ip),iCentMin,iCentMax,ip,ip);
    TH1D* hSigM0 = hSigCentPtM->ProjectionZ(Form("hSigM%i",ip),iCentMin,iCentMax,ip,ip);
    Double_t pt = hGenCentPtM->GetYaxis()->GetBinCenter(ip);
    Double_t widthPi0 = (4.13954+ 0.488129*pt + 0.260571*pt*pt)*0.001;
    Double_t widthEta = (7.59045+ 0.909524*pt + 0.939862*pt*pt)*0.001;
    Double_t mPi0Min = mPi0 - 2.*widthPi0;
    Double_t mPi0Max = mPi0 + 2.*widthPi0;
    Double_t mEtaMin = mEta - 2.*widthEta;
    Double_t mEtaMax = mEta + 2.*widthEta;
    Int_t imPi0Min = hRecCentPtM->GetZaxis()->FindFixBin(mPi0Min);
    Int_t imPi0Max = hRecCentPtM->GetZaxis()->FindFixBin(mPi0Max);
    Int_t imEtaMin = hRecCentPtM->GetZaxis()->FindFixBin(mEtaMin);
    Int_t imEtaMax = hRecCentPtM->GetZaxis()->FindFixBin(mEtaMax);
    double_t allPi0 = hRecM0->Integral(imPi0Min,imPi0Max);
    double_t sigPi0 = hSigM0->Integral(imPi0Min,imPi0Max);
    double_t allEta = hRecM0->Integral(imEtaMin,imEtaMax);
    double_t sigEta = hSigM0->Integral(imEtaMin,imEtaMax);
    double_t sig = isPi0 ? sigPi0 : sigEta;
    double_t err = isPi0 ? sqrt(allPi0) : sqrt(allEta);
    hPt->SetBinContent(ip,sig);
    hPt->SetBinError(ip,err);
  }

  new TCanvas;
  hPt->Draw();

  Double_t nAnalysed = hCent->Integral(iCentMin,iCentMax);
  printf("%f\n",nAnalysed);
  Double_t nExpected = nAnalysed;
  Double_t scaleFactor = Double_t(nExpected)/nAnalysed;

  TH1D* hPtScaled = (TH1D*) hPt->Clone("hPtScaled");
  for (Int_t i=1;i<=hPt->GetNbinsX();i++){
    hPtScaled->SetBinContent(i,hPt->GetBinContent(i)*scaleFactor);
    hPtScaled->SetBinError(i,hPt->GetBinError(i)*sqrt(scaleFactor));
  }

  TCanvas* c = new TCanvas("c","",800,600);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.06);
  gPad->SetBottomMargin(0.12);
  hPtScaled->SetMarkerStyle(kFullCircle);
  hPtScaled->SetMarkerSize(1.);
  hPtScaled->SetMarkerColor(kBlue);
  hPtScaled->SetLineColor(kBlue);
  hPtScaled->SetLineWidth(2);
  hPtScaled->GetXaxis()->SetTitleSize(0.05);
  hPtScaled->GetYaxis()->SetTitleSize(0.05);
  hPtScaled->GetXaxis()->SetTitleOffset(1.1);
  hPtScaled->GetYaxis()->SetTitleOffset(1.5);
  hPtScaled->GetXaxis()->SetLabelSize(0.05);
  hPtScaled->GetYaxis()->SetLabelSize(0.05);
  hPtScaled->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hPtScaled->GetYaxis()->SetTitle(Form("Expected %s yield",isPi0 ? "#pi^{0}" : "#eta"));

  hPtScaled->SetMaximum(hPtScaled->GetMaximum()*1.3);
  hPtScaled->Draw();
  hPtScaled->SetMinimum(0);

  TLatex* l = new TLatex();
  l->SetNDC(1);
  l->SetTextAlign(12);
  l->SetTextFont(42);
  l->SetTextSize(0.043);
  l->DrawLatex(0.5,0.88,"Bi-Bi, #sqrt{#it{s}_{NN}} = 9.2 GeV");
  l->DrawLatex(0.5,0.82,"10^{7} min. bias URQMD events");
  gPad->Print("expected_yield.png");

  TFile* fyield = new TFile("yield.root","recreate");
  hPtScaled->Write();
  fyield->Close();
}
