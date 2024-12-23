// Generic analysis for HELIOS MC simulation
// Y. Ayyad ayyadlim@frib.msu.edu 12/30/2020

void sim_analysis(Int_t num_ev=100000)
{

    TString mcFileNameHead = "heliossim";
    TString mcFileNameTail = ".root";
    TString mcFileName     = mcFileNameHead + mcFileNameTail;
    TString outFileNameHead = "heliosana";
    TString outFileNameTail = ".root";
    TString outFileName     = outFileNameHead + outFileNameTail;

    AtSiPoint* point = new AtSiPoint();
    TClonesArray *pointArray=0;

    TFile* file = new TFile(mcFileName.Data(),"READ");
    TTree* tree = (TTree*) file -> Get("cbmsim");


    tree = (TTree*) file -> Get("cbmsim");
    tree -> SetBranchAddress("AtSiArrayPoint", &pointArray);
    Int_t nEvents = tree -> GetEntriesFast();

    if(nEvents>num_ev) nEvents=num_ev;

    //Histograms
    TH2D *Eloss_vs_Z = new TH2D("Eloss_vs_Z","ELoss_vs_Z",1000,0,-1000,1000,0,10);
    TH2D *ElossTarget_vs_Z = new TH2D("ElossTarget_vs_Z","ELossTarget_vs_Z",1000,0,-1000,1000,0,10);
    TH2D *Recoil_XY = new TH2D("Recoil_XY","Recoil_XY",1000,-50.0,50.0,1000,-50.0,50.0);
    TH2D *ThetaLabSim_vs_Z = new TH2D("ThetaLabSim_vs_Z","ThetaLabSim_vs_Z",1000,0,-1000,1000,-180,180);

      TCanvas *c1 = new TCanvas();
    c1->Divide(2,2);
    c1->Draw();

    for(Int_t iEvent=0; iEvent<nEvents; iEvent++)
    {

	tree->GetEvent(iEvent);
        // tree -> GetEntry(iEvent);
        Int_t n = pointArray -> GetEntries();
        std::cout<<" ---- Event Number : "<<iEvent<<std::endl;
        double zpos_sca = 0.0;
        double energyLoss_sca = 0.0;
        double energyLossTarget_sca = 0.0;
        double zpos_rec = 0.0;
        double xpos_rec = 0.0;
        double ypos_rec = 0.0;
        double energyLoss_rec = 0.0;
	double angle_sca = 0.0;
	double energy_sca = 0.0;

           for(Int_t i=0; i<n; i++) {

            point = (AtSiPoint*) pointArray -> At(i);
            auto VolName=point->GetVolName();
            auto trackID = point -> GetTrackID();
            //std::cout<<" Volume Name : "<<VolName<<std::endl;
            //std::cout<<" Track ID : "<<trackID<<std::endl;

             if(trackID==1) //Scattered proton
	     {
              if(VolName!="siliconTarget"){
	       std::cout<<" Volume Name : "<<VolName<<std::endl;
               std::cout<<" Track ID : "<<trackID<<std::endl;
               std::cout<<" Point number : "<<i<<std::endl;
	       energyLoss_sca+=( point -> GetEnergyLoss() )*1000;//MeV
	       zpos_sca=point->GetZ()*10;
               angle_sca=point->GetAIni();
	       std::cout<<" Energy Loss "<<energyLoss_sca<<"  Z pos "<<zpos_sca<<" Point ID : "<<i<<"\n";
               }else if(VolName=="siliconTarget"){
                 
                 energyLossTarget_sca+=( point -> GetEnergyLoss() )*1000;//MeV;
               }

	     }

             if(trackID==0) //Heavy recoil
	     {
	      if(VolName!="siliconTarget"){
	       std::cout<<" Volume Name : "<<VolName<<std::endl;
               std::cout<<" Track ID : "<<trackID<<std::endl;
               std::cout<<" Point number : "<<i<<std::endl;
	       energyLoss_rec+=( point -> GetEnergyLoss() )*1000;//MeV
	       zpos_rec=point->GetZ()*10;
               xpos_rec=point->GetX(zpos_rec)*10;
               ypos_rec=point->GetY(zpos_rec)*10;
	       //std::cout<<" Energy Loss "<<energyLoss_rec<<"  Z pos "<<zpos_rec<<" Point ID : "<<i<<"\n";
              }
	     }
				

           }

	Eloss_vs_Z->Fill(zpos_sca,energyLoss_sca);
        Recoil_XY->Fill(xpos_rec,ypos_rec);
        ThetaLabSim_vs_Z->Fill(zpos_sca,angle_sca);

        ElossTarget_vs_Z->Fill(zpos_sca,energyLossTarget_sca);

    }
    
        std::map<TString, TH1D*> detectorHistograms;
    
    for(Int_t iEvent=0; iEvent<nEvents; iEvent++)
    {
        tree->GetEvent(iEvent);
        Int_t n = pointArray -> GetEntries();
        std::cout<<" ---- Event Number : "<<iEvent<<std::endl;

        // Map to store energy deposition for each detector in this event
        std::map<TString, Double_t> eventEnergyDeposition;

        for(Int_t i=0; i<n; i++) {
            point = (AtSiPoint*) pointArray -> At(i);
            auto VolName = point->GetVolName();
            auto trackID = point -> GetTrackID();
            Double_t energyLoss = (point -> GetEnergyLoss()) * 1000; // MeV

            // Accumulate energy deposition for each detector
            if (VolName != "siliconTarget") {
                eventEnergyDeposition[VolName] += energyLoss;

                // Create histogram for this detector if it doesn't exist
                if (detectorHistograms.find(VolName) == detectorHistograms.end()) {
                    TString histName = "EnergyDeposition_" + VolName;
                    detectorHistograms[VolName] = new TH1D(histName, histName, 1000, 0, 40);
                }
            }

            // ... [rest of the loop remains unchanged] ...
        }

        // Fill histograms with total energy deposition for each detector in this event
        for (const auto& pair : eventEnergyDeposition) {
            detectorHistograms[pair.first]->Fill(pair.second);
        }

        // ... [rest of the event loop remains unchanged] ...
    }

    int canvasCount = (detectorHistograms.size() + 3) / 4; // Round up division
    std::vector<TCanvas*> canvases;

    auto histIt = detectorHistograms.begin();
    for (int i = 0; i < canvasCount; ++i) {
        TString canvasName = TString::Format("c_EnergyDeposition_%d", i+1);
        TCanvas* canvas = new TCanvas(canvasName, canvasName, 800, 600);
        canvas->Divide(2, 2);
        canvases.push_back(canvas);

        for (int j = 0; j < 4 && histIt != detectorHistograms.end(); ++j, ++histIt) {
            canvas->cd(j+1);
            histIt->second->Draw();
        }
    }

    c1->cd(1);
    Eloss_vs_Z->Draw("zcol");
    c1->cd(2);
    Recoil_XY->Draw("zcol");
    c1->cd(3);
    ElossTarget_vs_Z->Draw("zcol");

}
