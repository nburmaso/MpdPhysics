#if !defined(__CINT__) && !defined(__CLING__)
#include "FairRunSim.h"
#include "FairPrimaryGenerator.h"
#include "FairBoxGenerator.h"
#include "FairRuntimeDb.h"
#endif
R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "TVirtualMC.h"
#include "macro/mpd/mpdloadlibs.C"
#include "geometry.C"

void sim(TString inFile = "urqmd.f14", TString outFile="sim.root", Int_t nEvents = 1000) {
  gDebug = 0;

  FairRunSim *fRun = new FairRunSim();
  fRun->SetPythiaDecayer(1);
  fRun->SetName("TGeant4");
  fRun->SetOutputFile(outFile);
  geometry(fRun);

  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();
  primGen->SetBeam(0.0,0.0,0.1,0.1);
  primGen->SetTarget(0.0,24.0);
  primGen->SmearGausVertexZ(kTRUE);
  primGen->SmearVertexXY(kTRUE);
  fRun->SetGenerator(primGen);
  fRun->SetUserDecay(kTRUE);

  MpdUrqmdGenerator* urqmdGen = new MpdUrqmdGenerator(inFile);
  urqmdGen->SetEventPlane(0,TMath::TwoPi());
  primGen->AddGenerator(urqmdGen);

  MpdMultiField* fField = new MpdMultiField();
  MpdConstField* fMagField = new MpdConstField();
  fMagField->SetField(0., 0., 5.);
  fMagField->SetFieldRegion(-230, 230, -230, 230, -375, 375);
  fField->AddField(fMagField);
  fRun->SetField(fField);

  fRun->Init();

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* output = new FairParRootFileIo(kTRUE);
  output->open(gFile);
  rtdb->setOutput(output);
  MpdMultiFieldPar* Par = (MpdMultiFieldPar*) rtdb->getContainer("MpdMultiFieldPar");
  Par->SetParameters(fField);
  Par->setInputVersion(fRun->GetRunId(), 1);
  Par->setChanged();
  rtdb->saveOutput();
  rtdb->print();

  fRun->Run(nEvents);
}
