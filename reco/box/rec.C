#if !defined(__CINT__) && !defined(__CLING__)
// ROOT includes
#include "TString.h"
#include "TStopwatch.h"
#include "TSystem.h"
// Fair includes
#include "FairRunAna.h"
#include "FairFileSource.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTask.h"
// MPD includes
#include "MpdTpcDigitizerAZ.h"
#include "MpdTpcClusterFinderMlem.h"
#include "MpdKalmanFilter.h"
#include "MpdVertexZfinder.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdKfPrimaryVertexFinder.h"
#include "MpdTofHitProducer.h"
#include "MpdTofMatching.h"
#include "MpdZdcDigiProducer.h"
#include "MpdFillDstTask.h"
#endif

void rec(TString inFile = "sim.root", TString outFile = "rec.root", Int_t nEvents = 200000) {
  TStopwatch timer;
  timer.Start();
  FairRunAna* fRun = new FairRunAna();

  FairSource* fFileSource = new FairFileSource(inFile);
  fRun->SetSource(fFileSource);
  fRun->SetOutputFile(outFile);
  fRun->SetGenerateRunInfo(false);
  fRun->SetUseFairLinks(true);

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parIo1 = new FairParRootFileIo();
  parIo1->open(inFile.Data());
  rtdb->setFirstInput(parIo1);
  rtdb->setOutput(parIo1);
  rtdb->saveOutput();

  MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
  fRun->AddTask(kalman);

  MpdTpcDigitizerAZ* tpcDigitizer = new MpdTpcDigitizerAZ();
  tpcDigitizer->SetPersistence(kFALSE);
  fRun->AddTask(tpcDigitizer);

  MpdTpcClusterFinderMlem *tpcClusAZ = new MpdTpcClusterFinderMlem();
  fRun->AddTask(tpcClusAZ);

  FairTask* vertZ = new MpdVertexZfinder();
  fRun->AddTask(vertZ);

  MpdTpcKalmanFilter* recoKF = new MpdTpcKalmanFilter("Kalman filter");
  recoKF->UseTpcHit(kFALSE);
  fRun->AddTask(recoKF);

  FairTask* findVtx = new MpdKfPrimaryVertexFinder("Vertex finder");
  fRun->AddTask(findVtx);

  MpdTofHitProducer* tofHit = new MpdTofHitProducer("Hit producer");
  fRun->AddTask(tofHit);

  MpdTofMatching* tofMatch = new MpdTofMatching("TOF matching");
  fRun->AddTask(tofMatch);

  MpdFillDstTask* fillDST = new MpdFillDstTask("MpdDst task");
  fRun->AddTask(fillDST);

  fRun->Init();

  fRun->Run(0,nEvents);

  timer.Stop();
  printf("Reco finished successfully: %.1f seconds\n",timer.RealTime());
}
