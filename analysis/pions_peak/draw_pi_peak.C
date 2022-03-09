#include "TCanvas.h"

void draw_pi_peak(double centMin = 0., double centMax = 100, double ptMin = 0.8, double ptMax = 2, Bool_t isPi0 = 1){
  gStyle->SetOptStat(0);
  double mPi0 = 0.1349766;
  double mEta = 0.547862;
  double mMes = isPi0 ? mPi0 : mEta;
  double mMinBgd   = isPi0 ? 0.15:  0.5;
  double mMaxBgd   = isPi0 ? 0.25:  0.6;
  double mMinSig   = isPi0 ? 0.125: 0.525;
  double mMaxSig   = isPi0 ? 0.16 : 0.59;
  double mMinRec   = isPi0 ? 0.05 : 0.4;
  double mMaxRec   = isPi0 ? 0.25 : 0.7;
  double mMinRange = isPi0 ? 0.00 : 0.4;
  double mMaxRange = isPi0 ? 0.26 : 0.7;

  //-------------------------------
  auto* f00 = new TFile("results.BiBi09.root");

  //-------------------------------
  auto* hRecCentPtM00 = (TH3D*) f00->Get("hRecCentPtM");
  auto* hSigCentPtM00 = (TH3D*) f00->Get("hSigCentPtM");
  hRecCentPtM00->Rebin3D(1,1,1);
  hSigCentPtM00->Rebin3D(1,1,1);
  int iCentMin00 =  hRecCentPtM00->GetXaxis()->FindFixBin(centMin+0.001);
  int iCentMax00 =  hRecCentPtM00->GetXaxis()->FindFixBin(centMax-0.001);
  int iPtMin00 = hRecCentPtM00->GetYaxis()->FindFixBin(ptMin+0.001);
  int iPtMax00 = hRecCentPtM00->GetYaxis()->FindFixBin(ptMax-0.001);
  auto* hRecM00 = hRecCentPtM00->ProjectionZ("hRecM00",iCentMin00,iCentMax00,iPtMin00,iPtMax00);
  auto* hSigM00 = hSigCentPtM00->ProjectionZ("hSigM00",iCentMin00,iCentMax00,iPtMin00,iPtMax00);


  //-------------------------------
  auto* fgaus00 = new TF1("fgaus00","gausn",0,0.2);
  fgaus00->SetParameter(0,hSigM00->Integral()*hSigM00->GetBinWidth(1));
  fgaus00->FixParameter(1,mMes);
  fgaus00->SetParameter(2,0.005);


  //-------------------------------
  auto* fbg = new TF1("fbg","pol2",0.15,0.25);
  fbg->SetParameter(0,1);
  fbg->SetParameter(1,1);
  fbg->SetParameter(2,1);
  fbg->SetLineColor(kBlack);

  //-------------------------------
  auto* fsum00 = new TF1("fsum00","gausn+pol2(3)",0.1,0.25);
  fsum00->SetParameter(0,fgaus00->GetParameter(0));
  fsum00->SetParameter(1,fgaus00->GetParameter(1));
  fsum00->SetParameter(2,fgaus00->GetParameter(2));
  fsum00->SetParameter(3,fbg->GetParameter(0));
  fsum00->SetParameter(4,fbg->GetParameter(1));
  fsum00->SetParameter(5,fbg->GetParameter(2));

  //-------------------------------
  hRecM00->SetTitle(Form(";m_{#gamma#gamma} (GeV/#it{c}^{2});Entries per %.1f MeV/#it{c}^{2}",hRecM00->GetBinWidth(1)*1000.));
  hRecM00->SetLineColor(kBlue);

  //-------------------------------
  auto* cM = new TCanvas("cM","cM",1200,1000);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.12);
  hRecM00->GetXaxis()->SetTitleSize(0.04);
  hRecM00->GetYaxis()->SetTitleSize(0.04);
  hRecM00->GetXaxis()->SetTitleFont(132);
  hRecM00->GetYaxis()->SetTitleFont(132);
  hRecM00->GetXaxis()->SetTitleOffset(1.1);
  hRecM00->GetYaxis()->SetTitleOffset(1.4);
  hRecM00->GetXaxis()->SetLabelSize(0.04);
  hRecM00->GetYaxis()->SetLabelSize(0.04);
  hRecM00->GetXaxis()->SetLabelFont(132);
  hRecM00->GetYaxis()->SetLabelFont(132);

  hRecM00->SetLineWidth(4);
  hRecM00->SetMinimum(0);
  hRecM00->GetXaxis()->SetRangeUser(mMinRange,mMaxRange);
  hRecM00->Sumw2();
  hRecM00->GetYaxis()->SetRangeUser(0.9,1.5*hRecM00->GetMaximum());

  hRecM00->Draw("same");
  hRecM00->Fit(fbg,"LQN","",mMinBgd,mMaxBgd);
  hSigM00->Fit(fgaus00,"QN","",mMinSig,mMaxSig);
  fsum00->SetParameter(0,fgaus00->GetParameter(0));
  fsum00->SetParameter(1,fgaus00->GetParameter(1));
  fsum00->SetParameter(2,fgaus00->GetParameter(2));
  fsum00->SetParameter(3,fbg->GetParameter(0));
  fsum00->SetParameter(4,fbg->GetParameter(1));
  fsum00->SetParameter(5,fbg->GetParameter(2));
  hRecM00->Fit(fsum00,"Q","",mMinRec,mMaxRec);

  //-------------------------------
  float p0_val00 = fgaus00->GetParameter(0);
  float p0_err00 = fsum00->GetParError(0);
  float w_val00  = fgaus00->GetParameter(2);
  float w_err00  = fgaus00->GetParError(2);
  double mMin00 = mMes-2*w_val00;
  double mMax00 = mMes+2*w_val00;
  int iMMin00 = hSigM00->GetXaxis()->FindFixBin(mMin00+0.0001);
  int iMMax00 = hSigM00->GetXaxis()->FindFixBin(mMax00-0.0001);
  double sigPi0_val00 = hSigM00->Integral(iMMin00,iMMax00);
  double bgdPi0_val00 = hRecM00->Integral(iMMin00,iMMax00)-sigPi0_val00;
  float n_val00  = sigPi0_val00;
  float n_err00  = n_val00/p0_val00*p0_err00;
  float signif00 = sigPi0_val00 / sqrt(sigPi0_val00 + bgdPi0_val00);

  printf("signal=%f\n", sigPi0_val00);

  TLegend* l = new TLegend(0.16,0.84,0.46,0.69);
  l->SetBorderSize(0);
  l->SetTextSize(0.04);
  l->SetTextFont(132);
  l->AddEntry(hRecM00,Form("#sigma_{#pi} = %.1f #pm %.1f MeV/#it{c}^{2}, S/#sqrt{S+B} = %.2f", 1000*w_val00, 1000*w_err00, signif00));
  l->Draw();

  auto* tex = new TLatex();
  tex->SetNDC(1);
  tex->SetTextAlign(12);
  tex->SetTextFont(132);
  tex->SetTextSize(0.04);
  tex->DrawLatex(0.18,0.93,"Bi#minusBi, #sqrt{#it{s}_{NN}} = 9.2 GeV, 10M URQMD events");
  tex->DrawLatex(0.18,0.87,Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}^{2}, |y| < 1",ptMin,ptMax));
}
