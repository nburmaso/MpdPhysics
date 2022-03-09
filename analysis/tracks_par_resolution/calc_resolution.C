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

int nHitsMin=20;
float ptMin=0.05;
float etaMax=1.0;
float chi2Max=10;
float psiMax=10;
float rMin=10;
float yPiMax=1.;

int n_bins_pt = 50;
double pt_low = 0;
double pt_hi = 2;
int n_bins_pres = 1000;
double pres_low = -1;
double pres_hi = 1;
int n_bins_rres = 4000;
double rres_low = 0;
double rres_hi = 20;
int n_bins_zres = 1000;
double zres_low = -20;
double zres_hi = 20;

void calc_resolution()
{
  TChain *tree = new TChain("convTree");
  tree->Add("conv.root");
  SetBranchAddresses(tree);
  int nEvents = tree->GetEntries();
  printf("nEvents = %i\n", nEvents);

  // e+
  auto *hPtRes_ep  = new TH2D("hPtRes_ep", "",  n_bins_pt, pt_low, pt_hi, n_bins_pres, pres_low, pres_hi);
  auto *hDCARPt_ep = new TH2D("hDCARPt_ep", "", n_bins_pt, pt_low, pt_hi, n_bins_rres, rres_low, rres_hi);
  auto *hDCAZPt_ep = new TH2D("hDCAZPt_ep", "", n_bins_pt, pt_low, pt_hi, n_bins_zres, zres_low, zres_hi);

  // e-
  auto *hPtRes_em  = new TH2D("hPtRes_em", "",  n_bins_pt, pt_low, pt_hi, n_bins_pres, pres_low, pres_hi);
  auto *hDCARPt_em = new TH2D("hDCARPt_em", "", n_bins_pt, pt_low, pt_hi, n_bins_rres, rres_low, rres_hi);
  auto *hDCAZPt_em = new TH2D("hDCAZPt_em", "", n_bins_pt, pt_low, pt_hi, n_bins_zres, zres_low, zres_hi);

  // pi+
  auto *hPtRes_pip  = new TH2D("hPtRes_pip", "",  n_bins_pt, pt_low, pt_hi, n_bins_pres, pres_low, pres_hi);
  auto *hDCARPt_pip = new TH2D("hDCARPt_pip", "", n_bins_pt, pt_low, pt_hi, n_bins_rres, rres_low, rres_hi);
  auto *hDCAZPt_pip = new TH2D("hDCAZPt_pip", "", n_bins_pt, pt_low, pt_hi, n_bins_zres, zres_low, zres_hi);

  // pi-
  auto *hPtRes_pim  = new TH2D("hPtRes_pim", "",  n_bins_pt, pt_low, pt_hi, n_bins_pres, pres_low, pres_hi);
  auto *hDCARPt_pim = new TH2D("hDCARPt_pim", "", n_bins_pt, pt_low, pt_hi, n_bins_rres, rres_low, rres_hi);
  auto *hDCAZPt_pim = new TH2D("hDCAZPt_pim", "", n_bins_pt, pt_low, pt_hi, n_bins_zres, zres_low, zres_hi);

  // p+
  auto *hPtRes_pp  = new TH2D("hPtRes_pp", "",  n_bins_pt, pt_low, pt_hi, n_bins_pres, pres_low, pres_hi);
  auto *hDCARPt_pp = new TH2D("hDCARPt_pp", "", n_bins_pt, pt_low, pt_hi, n_bins_rres, rres_low, rres_hi);
  auto *hDCAZPt_pp = new TH2D("hDCAZPt_pp", "", n_bins_pt, pt_low, pt_hi, n_bins_zres, zres_low, zres_hi);

  // p-
  auto *hPtRes_pm  = new TH2D("hPtRes_pm", "",  n_bins_pt, pt_low, pt_hi, n_bins_pres, pres_low, pres_hi);
  auto *hDCARPt_pm = new TH2D("hDCARPt_pm", "", n_bins_pt, pt_low, pt_hi, n_bins_rres, rres_low, rres_hi);
  auto *hDCAZPt_pm = new TH2D("hDCAZPt_pm", "", n_bins_pt, pt_low, pt_hi, n_bins_zres, zres_low, zres_hi);

  map<int, int> mcTrackId_to_iMc;

  for (int iEv = 0; iEv < nEvents; iEv++)
  {
    if (iEv % 10000 == 0)
      printf("Event: %i\n", iEv);
    tree->GetEntry(iEv);

    mcTrackId_to_iMc.clear();
    for (int i = 0; i < nMc; i++)
      mcTrackId_to_iMc[vMcTrackId[i]] = i;

    for (int iRc = 0; iRc < nRc; iRc++)
    {
      int iMc = mcTrackId_to_iMc[vMcTrackRef[iRc]];
      int pdg = vPdgCode[iMc];
      int moId = vMotherId[iMc];
      TVector3 p(vPxRcRefit[iRc], vPyRcRefit[iRc], vPzRcRefit[iRc]);
      if (vNhits[iRc] < nHitsMin)
        continue;
      if (p.Pt() < ptMin)
        continue;
      if (fabs(p.Eta()) > etaMax)
        continue;
      float pt_mc = std::hypot(vPxMc[iMc], vPyMc[iMc]);
      float pt_res = (p.Pt() - pt_mc) / pt_mc;
      float dca_r = std::hypot(vDcaX[iRc], vDcaY[iRc]);
      float dca_z = vDcaZ[iRc];
      if (pdg == 11)
      {
        hPtRes_ep->Fill(p.Pt(), pt_res);
        hDCARPt_ep->Fill(p.Pt(), dca_r);
        hDCAZPt_ep->Fill(p.Pt(), dca_z);
      }
      if (pdg == -11)
      {
        hPtRes_em->Fill(p.Pt(), pt_res);
        hDCARPt_em->Fill(p.Pt(), dca_r);
        hDCAZPt_em->Fill(p.Pt(), dca_z);
      }
      if (moId == -1 && pdg == 211)
      {
        hPtRes_pip->Fill(p.Pt(), pt_res);
        hDCARPt_pip->Fill(p.Pt(), dca_r);
        hDCAZPt_pip->Fill(p.Pt(), dca_z);
      }
      if (moId == -1 && pdg == -211)
      {
        hPtRes_pim->Fill(p.Pt(), pt_res);
        hDCARPt_pim->Fill(p.Pt(), dca_r);
        hDCAZPt_pim->Fill(p.Pt(), dca_z);
      }
      if (moId == -1 && pdg == 2212)
      {
        hPtRes_pp->Fill(p.Pt(), pt_res);
        hDCARPt_pp->Fill(p.Pt(), dca_r);
        hDCAZPt_pp->Fill(p.Pt(), dca_z);
      }
      if (moId == -1 && pdg == -2212)
      {
        hPtRes_pm->Fill(p.Pt(), pt_res);
        hDCARPt_pm->Fill(p.Pt(), dca_r);
        hDCAZPt_pm->Fill(p.Pt(), dca_z);
      }
    }
  }

  auto *f = new TFile("resolution.root", "recreate");
  hPtRes_ep->Write();
  hDCARPt_ep->Write();
  hDCAZPt_ep->Write();
  hPtRes_em->Write();
  hDCARPt_em->Write();
  hDCAZPt_em->Write();
  hPtRes_pip->Write();
  hDCARPt_pip->Write();
  hDCAZPt_pip->Write();
  hPtRes_pim->Write();
  hDCARPt_pim->Write();
  hDCAZPt_pim->Write();
  hPtRes_pp->Write();
  hDCARPt_pp->Write();
  hDCAZPt_pp->Write();
  hPtRes_pm->Write();
  hDCARPt_pm->Write();
  hDCAZPt_pm->Write();
  f->Close();
}
