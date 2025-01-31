void B10_dp_sim_alpha(Int_t nEvents = 10000, TString mcEngine = "TGeant4")
{

  TString dir = getenv("VMCWORKDIR");

  // Output file name
  TString outFile ="heliossim.root";

  // Parameter file name
  TString parFile="heliospar.root";

  // -----   Timer   --------------------------------------------------------
  TStopwatch timer;
  timer.Start();
  // ------------------------------------------------------------------------
  
  // -----   Create simulation run   ----------------------------------------
  FairRunSim* run = new FairRunSim();
  run->SetName(mcEngine);              // Transport engine
  run->SetOutputFile(outFile);          // Output file
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  // ------------------------------------------------------------------------


  // -----   Create media   -------------------------------------------------
  run->SetMaterials("media.geo");       // Materials
  // ------------------------------------------------------------------------

  // -----   Create geometry   ----------------------------------------------

  FairModule* cave= new AtCave("CAVE");
  cave->SetGeometryFileName("cave.geo");
  run->AddModule(cave);
  

/*
 FairDetector *ATTPC = new AtTpc("ATTPC", kTRUE);
   ATTPC->SetGeometryFileName("ATTPC_He600torr_v2.root");
   // ATTPC->SetModifyGeometry(kTRUE);
   run->AddModule(ATTPC);
*/
   

  FairDetector* HELIOS = new AtSiArray("HELIOS", kTRUE);
  HELIOS->SetGeometryFileName("HELIOS_SiArray_v1.0.root");
  //HELIOS->SetModifyGeometry(kTRUE);
  run->AddModule(HELIOS);

 // ------------------------------------------------------------------------

    // -----   Magnetic field   -------------------------------------------
    // Constant Field
    AtConstField  *fMagField = new AtConstField();
    fMagField->SetField(0., 0. ,-22.5 ); // values are in kG  
    fMagField->SetFieldRegion(-300, 300,-300, 300, -500, 500); // values are in cm 
    fMagField->Print();
                          //  (xmin,xmax,ymin,ymax,zmin,zmax)
    run->SetField(fMagField);
    // --------------------------------------------------------------------



  // -----   Create PrimaryGenerator   --------------------------------------
  FairPrimaryGenerator* primGen = new FairPrimaryGenerator();




                  // Beam Information
                Int_t z = 5;  // Atomic number
	        Int_t a = 10; // Mass number
	        Int_t q = 0;   // Charge State
	        Int_t m = 1;   // Multiplicity  NOTE: Due the limitation of the TGenPhaseSpace accepting only pointers/arrays the maximum multiplicity has been set to 10 particles. 
	        Double_t px = 0.000/a;  // X-Momentum / per nucleon!!!!!!
	        Double_t py = 0.000/a;  // Y-Momentum / per nucleon!!!!!! 
	        Double_t pz = 1.370/a;  // Z-Momentum / per nucleon!!!!!!  
  	        Double_t BExcEner = 0.0; 
                Double_t Bmass = 10.012937; //Mass in GeV
                Double_t NomEnergy = 1.2; //Used to force the beam to stop within a certain energy range.  

	          // AtTPCIonGenerator* ionGen = new AtTPCIonGenerator("Ion",z,a,q,m,px,py,pz,BExcEner,Bmass,NomEnergy);
	          // ionGen->SetSpotRadius(0,-100,0);
	          // add the ion generator

	          // primGen->AddGenerator(ionGen);

  		  //primGen->SetBeam(1,1,0,0); //These parameters change the position of the vertex of every track added to the Primary Generator
		  // primGen->SetTarget(30,0);

		 // Variables for 2-Body kinematics reaction
                  std::vector<Int_t> Zp; // Zp
		  std::vector<Int_t> Ap; // Ap
                  std::vector<Int_t> Qp;//Electric charge
                  Int_t mult;  //Number of particles
 		  std::vector<Double_t> Pxp; //Px momentum X
		  std::vector<Double_t> Pyp; //Py momentum Y
		  std::vector<Double_t> Pzp; //Pz momentum Z
                  std::vector<Double_t> Mass; // Masses
		  std::vector<Double_t> ExE; // Excitation energy
 		  Double_t ResEner; // Energy of the beam (Useless for the moment)


		  // Note: Momentum will be calculated from the phase Space according to the residual energy of the beam


	          mult = 4; //Number of Nuclei involved in the reaction (Should be always 4) THIS DEFINITION IS MANDATORY (and the number of particles must be the same)
                  ResEner = 140.0; // For fixed target mode (Si Array) in MeV 

                  // ---- Beam ----
                  Zp.push_back(z); // 40Ar TRACKID=0
            	  Ap.push_back(a); //
            	  Qp.push_back(q);
                  Pxp.push_back(px);
            	  Pyp.push_back(py);
            	  Pzp.push_back(pz);
            	  Mass.push_back(Bmass);
            	  ExE.push_back(BExcEner);

                  // ---- Target ----
                 Zp.push_back(1); // 
		 Ap.push_back(2); //
		 Qp.push_back(0); //
		 Pxp.push_back(0.0);
                 Pyp.push_back(0.0);
		 Pzp.push_back(0.0);
                 Mass.push_back(2.0135532);
		 ExE.push_back(0.0);//In MeV

                  //--- Scattered -----
                Zp.push_back(5); //
                Ap.push_back(11); //
                Qp.push_back(0);
          	Pxp.push_back(0.0);
          	Pyp.push_back(0.0);
          	Pzp.push_back(0.0);
          	Mass.push_back(11.009305);
          	ExE.push_back(9.2); //


                 // ---- Recoil -----
		 Zp.push_back(1); // p  
		 Ap.push_back(1); //
		 Qp.push_back(0); //
		 Pxp.push_back(0.0);
                 Pyp.push_back(0.0);
		 Pzp.push_back(0.0);
                 Mass.push_back(1.00783);
		 ExE.push_back(0.0);

                 Double_t ThetaMinCMS = 10.0;
                 Double_t ThetaMaxCMS = 45.0;

                 AtTPC2Body *TwoBody = new AtTPC2Body("TwoBody", &Zp, &Ap, &Qp, mult, &Pxp, &Pyp, &Pzp, &Mass, &ExE,
                                                      ResEner, ThetaMinCMS, ThetaMaxCMS);
                 TwoBody->SetFixedTargetPosition(0.0, 0.0, 0.0);
                 TwoBody->SetFixedBeamMomentum(0.0, 0.0, pz * a);
                 primGen->AddGenerator(TwoBody);

                 // Setting decay
   // Set the parameters of the decay generator


   std::vector<std::vector<Int_t>> zDecay;
   std::vector<std::vector<Int_t>> aDecay;
   std::vector<std::vector<Int_t>> qDecay;
   std::vector<std::vector<Double_t>> massDecay;

   Int_t zB;
   Int_t aB;
   Double_t massDecayB;
   Double_t massTarget;
   Double_t exEnergy;
   std::vector<Double_t> SepEne;

   Int_t TotDecayCases = 1; // the number of decay channel (case) to be considered

   zDecay.resize(TotDecayCases);
   aDecay.resize(TotDecayCases);
   qDecay.resize(TotDecayCases);
   massDecay.resize(TotDecayCases);

   zB = 5; // 11B
   aB = 11;
   massDecayB = 11.009305;
   massTarget = 0.0;
   exEnergy = 0.0; // NB: Set to zero for sequential decay
// Energía de separación de \( ^{11}B \) para el decaimiento en 7Li + α (8.88 MeV)
SepEne.push_back(8.664); // Energía de separación en MeV

// 7Li
zDecay.at(0).push_back(3);
aDecay.at(0).push_back(7);
qDecay.at(0).push_back(0);
massDecay.at(0).push_back(6.94174064); // masa de 7Li en u

// alpha
zDecay.at(0).push_back(2);
aDecay.at(0).push_back(4);
qDecay.at(0).push_back(0);
massDecay.at(0).push_back(4.00260325415);      // masa de alpha en u

AtTPCIonDecay *decay =
      new AtTPCIonDecay(&zDecay, &aDecay, &qDecay, &massDecay, zB, aB, massDecayB, massTarget, exEnergy, &SepEne);
   decay->SetSequentialDecay(kTRUE);
   primGen->AddGenerator(decay);


                 run->SetGenerator(primGen);

                 // ------------------------------------------------------------------------

                 //---Store the visualiztion info of the tracks, this make the output file very large!!
                 //--- Use it only to display but not for production!
                 //run->SetStoreTraj(kTRUE);

                 // -----   Initialize simulation run   ------------------------------------
                 run->Init();
                 // ------------------------------------------------------------------------

                 // Trajectory filters

                 // -----   Runtime database   ---------------------------------------------

                 Bool_t kParameterMerged = kTRUE;
                 FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
                 parOut->open(parFile.Data());
                 rtdb->setOutput(parOut);
                 rtdb->saveOutput();
                 rtdb->print();
                 // ------------------------------------------------------------------------

                 // -----   Start run   ----------------------------------------------------
                 run->Run(nEvents);

                 // You can export your ROOT geometry ot a separate file
                 //run->CreateGeometryFile("geofile_helios_full.root");
                 // ------------------------------------------------------------------------

                 // -----   Finish   -------------------------------------------------------
                 timer.Stop();
                 Double_t rtime = timer.RealTime();
                 Double_t ctime = timer.CpuTime();
                 cout << endl << endl;
                 cout << "Macro finished succesfully." << endl;
                 cout << "Output file is " << outFile << endl;
                 cout << "Parameter file is " << parFile << endl;
                 cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
                 // ------------------------------------------------------------------------
}
