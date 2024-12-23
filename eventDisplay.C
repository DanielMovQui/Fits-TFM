void eventDisplay()
{
  //-----User Settings:-----------------------------------------------
  TString  InputDataFile     ="heliossim.root";
  TString  ParFile       ="heliospar.root";
  TString  OutputDataFile	 ="heliostest.root";


  // -----   Reconstruction run   -------------------------------------------
  FairRunAna *fRun= new FairRunAna();
  FairRootFileSink *sink = new FairRootFileSink(OutputDataFile);
  FairFileSource *source = new FairFileSource(InputDataFile);
  fRun->SetSource(source);
  fRun->SetSink(sink);

  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(ParFile.Data());
  rtdb->setFirstInput(parInput1);

  FairEventManager *fMan= new FairEventManager();

  //----------------------Traks and points -------------------------------------
  FairMCPointDraw *AtSiArrayPoints = new FairMCPointDraw("AtSiArrayPoint", kBlue, kFullSquare);

  fMan->AddTask(AtSiArrayPoints);

  fMan->Init();

}
