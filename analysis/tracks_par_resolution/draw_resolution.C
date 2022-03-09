#include "TSystem.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

double maxPt = 2;
double minPt = 0;
int nBins = 40;
int start = 2;
int rebin_factor = 1;

double fitf(double *x, double *par)
{
    double arg = 0;
    if (par[2] != 0)
    {
      arg = (x[0] - par[1])/par[2];
    }
    double fitval = par[0]*exp(-0.5*arg*arg);
    return fitval;
}

void draw_resolution()
{
  gStyle->SetOptStat(0);

  auto *f = new TFile("resolution.root");
  auto *hPtRes_pip = (TH2D *)f->Get("hPtRes_pip");
  auto *hDCAZPt_pip = (TH2D *)f->Get("hDCAZPt_pip");
  auto *hDCARPt_pip = (TH2D *)f->Get("hDCARPt_pip");

  hPtRes_pip->Sumw2();
  hDCAZPt_pip->Sumw2();
  hDCARPt_pip->Sumw2();

  if (rebin_factor > 1)
  {
    hPtRes_pip->RebinX(rebin_factor);
    hDCAZPt_pip->RebinX(rebin_factor);
    hDCARPt_pip->RebinX(rebin_factor);
  }

  hDCARPt_pip->RebinY(4);

  // pt
  auto *hPtSigmas_pip = new TH1D("hPtSigmas_pip", "#pi^{+}; p_{t}; (p_{t}-p_{t,MC})/p_{t}", nBins, minPt, maxPt);
  // dca z
  auto *hDCAZSigmas_pip = new TH1D("hDCAZSigmas_pip", "#pi^{+}; p_{t}; DCA_{z}", nBins, minPt, maxPt);
  // dca r
  auto *hDCARSigmas_pip = new TH1D("hDCARSigmas_pip", "#pi^{+}; p_{t}; DCA_{r}", nBins, minPt, maxPt);

  // pt
  // -----------------------------------------------

  // pi
  for (int i = start; i < nBins; i++)
  {
    auto* slice_pt = hPtRes_pip->ProjectionY("slice_pt", i, i); // slice along pt bin number i
    auto slice_pt_fit = slice_pt->Fit("gaus", "S");
    float sigma = slice_pt_fit->Parameter(2);
    hPtSigmas_pip->SetBinContent(i, sigma);
    delete slice_pt;
  }

  // dca z
  // -----------------------------------------------

  // pi
  for (int i = start; i < nBins; i++)
  {
    auto* slice_pt = hDCAZPt_pip->ProjectionY("slice_pt", i, i); // slice along pt bin number i
    auto slice_pt_fit = slice_pt->Fit("gaus", "S");
    float sigma = slice_pt_fit->Parameter(2);
    hDCAZSigmas_pip->SetBinContent(i, sigma);
    delete slice_pt;
  }

  // dca
  // -----------------------------------------------

  // pi
  for (int i = start; i < nBins; i++)
  {
    auto* slice_pt = hDCARPt_pip->ProjectionY("slice_pt", i, i); // slice along pt bin number i
    // double h_sdev = slice_pt->GetStdDev();
    double h_mpv = slice_pt->GetMaximum() * (maxPt - minPt) / nBins + minPt;
    auto *fit_gaus = new TF1("fit_gaus", fitf, 0, 20, 3);
    fit_gaus->SetParNames("Constant","Mean_value","Sigma");
    fit_gaus->SetParameters(0, slice_pt->GetMean(), slice_pt->GetRMS());
    auto slice_pt_fit = slice_pt->Fit("fit_gaus", "S");
    float sigma = abs(slice_pt_fit->Parameter(2));
    hDCARSigmas_pip->SetBinContent(i, sigma);
    delete slice_pt;
  }

  // -----------------------------------------------
  maxPt = 1.5 - 1e-3;

  hPtSigmas_pip->GetXaxis()->SetRangeUser(minPt, maxPt);
  hDCAZSigmas_pip->GetXaxis()->SetRangeUser(minPt, maxPt);
  hDCARSigmas_pip->GetXaxis()->SetRangeUser(minPt, maxPt);


  TCanvas* c1 = new TCanvas("c1","",1200,1000);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.12);
  hPtSigmas_pip->GetYaxis()->SetRangeUser(0.001, 0.059);
  hPtSigmas_pip->SetTitle("; #it{p}_{T} (GeV/#it{c}); #it{p}_{T} resolution");
  hPtSigmas_pip->SetLineWidth(2);
  hPtSigmas_pip->SetLineColor(kBlack);
  hPtSigmas_pip->SetMarkerColor(kBlack);
  hPtSigmas_pip->GetXaxis()->SetTitleSize(0.04);
  hPtSigmas_pip->GetYaxis()->SetTitleSize(0.04);
  hPtSigmas_pip->GetXaxis()->SetTitleFont(132);
  hPtSigmas_pip->GetYaxis()->SetTitleFont(132);
  hPtSigmas_pip->GetXaxis()->SetTitleOffset(1.1);
  hPtSigmas_pip->GetYaxis()->SetTitleOffset(1.4);
  hPtSigmas_pip->GetXaxis()->SetLabelSize(0.04);
  hPtSigmas_pip->GetYaxis()->SetLabelSize(0.04);
  hPtSigmas_pip->GetXaxis()->SetLabelFont(132);
  hPtSigmas_pip->GetYaxis()->SetLabelFont(132);
  hPtSigmas_pip->SetMinimum(0);
  hPtSigmas_pip->GetYaxis()->SetDecimals(1);
  hPtSigmas_pip->Draw();

  TCanvas* c2 = new TCanvas("c2","",1200,1000);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.12);
  hDCAZSigmas_pip->GetYaxis()->SetRangeUser(0, 2);
  hDCAZSigmas_pip->SetTitle("; #it{p}_{T} (GeV/#it{c}); DCA_{z} resolution (cm)");
  hDCAZSigmas_pip->SetLineWidth(2);
  hDCAZSigmas_pip->SetLineColor(kBlack);
  hDCAZSigmas_pip->SetMarkerColor(kBlack);
  hDCAZSigmas_pip->GetXaxis()->SetTitleSize(0.04);
  hDCAZSigmas_pip->GetYaxis()->SetTitleSize(0.04);
  hDCAZSigmas_pip->GetXaxis()->SetTitleFont(132);
  hDCAZSigmas_pip->GetYaxis()->SetTitleFont(132);
  hDCAZSigmas_pip->GetXaxis()->SetTitleOffset(1.1);
  hDCAZSigmas_pip->GetYaxis()->SetTitleOffset(1.4);
  hDCAZSigmas_pip->GetXaxis()->SetLabelSize(0.04);
  hDCAZSigmas_pip->GetYaxis()->SetLabelSize(0.04);
  hDCAZSigmas_pip->GetXaxis()->SetLabelFont(132);
  hDCAZSigmas_pip->GetYaxis()->SetLabelFont(132);
  hDCAZSigmas_pip->SetMinimum(0);
  hDCAZSigmas_pip->GetYaxis()->SetDecimals(1);
  hDCAZSigmas_pip->Draw();

  TCanvas* c3 = new TCanvas("c3","",1200,1000);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.12);
  hDCARSigmas_pip->GetYaxis()->SetRangeUser(0, 2);
  hDCARSigmas_pip->SetTitle("; #it{p}_{T} (GeV/#it{c}); DCA_{r} resolution (cm)");
  hDCARSigmas_pip->SetLineWidth(2);
  hDCARSigmas_pip->SetLineColor(kBlack);
  hDCARSigmas_pip->SetMarkerColor(kBlack);
  hDCARSigmas_pip->GetXaxis()->SetTitleSize(0.04);
  hDCARSigmas_pip->GetYaxis()->SetTitleSize(0.04);
  hDCARSigmas_pip->GetXaxis()->SetTitleFont(132);
  hDCARSigmas_pip->GetYaxis()->SetTitleFont(132);
  hDCARSigmas_pip->GetXaxis()->SetTitleOffset(1.1);
  hDCARSigmas_pip->GetYaxis()->SetTitleOffset(1.4);
  hDCARSigmas_pip->GetXaxis()->SetLabelSize(0.04);
  hDCARSigmas_pip->GetYaxis()->SetLabelSize(0.04);
  hDCARSigmas_pip->GetXaxis()->SetLabelFont(132);
  hDCARSigmas_pip->GetYaxis()->SetLabelFont(132);
  hDCARSigmas_pip->SetMinimum(0);
  hDCARSigmas_pip->GetYaxis()->SetDecimals(1);
  hDCARSigmas_pip->Draw();

  c1->Print("pt_resolution.png");
  c2->Print("dcaz_resolution.png");
  c3->Print("dcar_resolution.png");
}
