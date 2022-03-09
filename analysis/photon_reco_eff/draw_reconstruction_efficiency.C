//R__ADD_INCLUDE_PATH($VMCWORKDIR/install/include)
#include "TSystem.h"
#include <TVector3.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "Rtypes.h"
#include "FairMCTracks.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TStyle.h"

void draw_reconstruction_efficiency() {
  gStyle->SetOptStat(0);
  
  TFile* f = new TFile("reconstruction_efficiency.root");
  TH1D* hPtConv  = (TH1D*) f->Get("hPtConv");
  TH1D* hPtAll   = (TH1D*) f->Get("hPtAll");
  TH1D* hPtAcc   = (TH1D*) f->Get("hPtAcc");
  TH1D* hPtRec   = (TH1D*) f->Get("hPtRec");
  TH1D* hPtRecGa = (TH1D*) f->Get("hPtRecGa");
  hPtConv->Sumw2();
  hPtAll->Sumw2();
  hPtAcc->Sumw2();
  hPtRec->Sumw2();
  hPtRecGa->Sumw2();
  hPtConv->Divide(hPtConv,hPtAll,1,1,"B");
  hPtAcc->Divide(hPtAcc,hPtAll,1,1,"B");
  hPtRec->Divide(hPtRec,hPtAll,1,1,"B");
  hPtRecGa->Divide(hPtRecGa,hPtAll,1,1,"B");

  TCanvas* c1 = new TCanvas("c1","",800,600);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.12);
  hPtConv->SetTitle("; #it{p}_{T} (GeV/#it{c}); Efficiency");
  hPtConv->SetLineWidth(2);
  hPtAcc->SetLineWidth(2);
  hPtRec->SetLineWidth(2);
  hPtRecGa->SetLineWidth(2);
  hPtConv->SetLineColor(kBlack);
  hPtAcc->SetLineColor(kRed);
  hPtRec->SetLineColor(kMagenta);
  hPtRecGa->SetLineColor(kBlue);
  hPtConv->SetMarkerColor(kBlack);
  hPtAcc->SetMarkerColor(kRed);
  hPtRec->SetMarkerColor(kMagenta);
  hPtRecGa->SetMarkerColor(kBlue);
  hPtConv->GetXaxis()->SetRangeUser(0,1.5-0.001);
  hPtConv->GetXaxis()->SetTitleSize(0.05);
  hPtConv->GetYaxis()->SetTitleSize(0.05);
  hPtConv->GetXaxis()->SetTitleOffset(1.1);
  hPtConv->GetYaxis()->SetTitleOffset(1.4);
  hPtConv->GetXaxis()->SetLabelSize(0.05);
  hPtConv->GetYaxis()->SetLabelSize(0.05);
  hPtConv->SetMinimum(0);
  hPtConv->SetMaximum(0.04);
  hPtConv->GetYaxis()->SetDecimals(1);
  hPtConv->Draw();
  hPtAcc->Draw("same");
  hPtRec->Draw("same");
  hPtRecGa->Draw("same");

  TLegend* l = new TLegend(0.41,0.15,0.98,0.37);
  l->SetBorderSize(0);
  l->AddEntry(hPtConv,"Converted #gamma #rightarrow e^{+}e^{-}, |#eta_{#gamma}| < 1");
  l->AddEntry(hPtAcc,"Accepted #gamma, |#eta_{e}| < 1, #it{p}^{e}_{T} > 50 MeV/c");
  l->AddEntry(hPtRec,"Reconstructed #gamma, electron cuts");
  l->AddEntry(hPtRecGa,"Reconstructed #gamma, pair cuts");
  l->Draw();

  TFile* feff = new TFile("eff_gamma.root","recreate");
  hPtRecGa->Write("hEffPtGa");
  feff->Close();
}
