void sim_analysis_alpha_PID(Int_t num_ev=10000)

{
    TString mcFileNameHead = "heliossim";
    TString mcFileNameTail = ".root";
    TString mcFileName     = mcFileNameHead + mcFileNameTail;
    TString outFileNameHead = "heliosana";
    TString outFileNameTail = ".root";
    TString outFileName     = outFileNameHead + outFileNameTail;

    AtSiPoint* point = new AtSiPoint();
    TClonesArray *pointArray = nullptr;

    TFile* file = new TFile(mcFileName.Data(), "READ");
    TTree* tree = (TTree*) file->Get("cbmsim");
    tree->SetBranchAddress("AtSiArrayPoint", &pointArray);
    Int_t nEvents = tree->GetEntriesFast();

    if (nEvents > num_ev) nEvents = num_ev;

    Int_t totalEvents = 0;
    Int_t protonDetectedEvents = 0;
    Int_t lithiumDetectedEvents = 0;
    Int_t alphaDetectedEvents = 0;
    Int_t coincidenceEvents = 0;

    // Crear histogramas PID combinados para cada sector del detector QQQ
    TH2F* hPID_Combined[4];
    for (int i = 0; i < 4; i++) {
        hPID_Combined[i] = new TH2F(Form("hPID_Combined_QQQ%d", i+1), 
                                    Form("PID Recoil %d;E_{res} (MeV);#DeltaE (MeV)", i+1), 
                                    1000, 0, 50, 1000, 0, 50);
    }

    for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
        tree->GetEvent(iEvent);
        Int_t n = pointArray->GetEntries();

        bool protonDetected = false;
        bool lithiumDetected = false;
        bool alphaDetected = false;

        Double_t dE_Li[4] = {0}, E_res_Li[4] = {0};
        Double_t dE_Alpha[4] = {0}, E_res_Alpha[4] = {0};

        for (Int_t i = 0; i < n; i++) {
            point = (AtSiPoint*) pointArray->At(i);
            auto VolName = point->GetVolName();
            auto trackID = point->GetTrackID();
            auto energy = point->GetEnergyLoss(); // Convertir a MeV

            if (trackID == 1) { // Supongamos que trackID 1 es protón dispersado
                if (VolName.Contains("siliconYZL") || VolName.Contains("siliconYZR") || 
                    VolName.Contains("siliconXZT") || VolName.Contains("siliconXZB")) {
                    protonDetected = true;
                }
            }

            int sector = -1;
            if (VolName.Contains("Sector1")) sector = 0;
            else if (VolName.Contains("Sector2")) sector = 1;
            else if (VolName.Contains("Sector3")) sector = 2;
            else if (VolName.Contains("Sector4")) sector = 3;

            if (sector >= 0) {
                if (VolName.Contains("QQQ")) {
                    if (trackID == 0) { // 7Li
                        dE_Li[sector] += energy*1000 ; // aqui hay que modificar
                        lithiumDetected = true;
                    } else if (trackID == 2) { // Alpha
                        dE_Alpha[sector] += energy *1000; // y aqui
                        alphaDetected = true;
                    }
                } else if (VolName.Contains("siliconE")) {
                    if (trackID == 0) { // 7Li
                        E_res_Li[sector] += energy*1000 ; // y aqui
                    } else if (trackID == 2) { // Alpha
                        E_res_Alpha[sector] += energy*1000 ; // y aqui
                    }
                }
            }
        }

        // Llenar histogramas PID combinados
        for (int i = 0; i < 4; i++) {
            if (dE_Li[i] > 0 && E_res_Li[i] > 0) {
                hPID_Combined[i]->Fill(E_res_Li[i], dE_Li[i], 1); // Color rojo para litio
            }
            if (dE_Alpha[i] > 0 && E_res_Alpha[i] > 0) {
                hPID_Combined[i]->Fill(E_res_Alpha[i], dE_Alpha[i], 2); // Color azul para alpha
            }
        }

        totalEvents++;
        if (protonDetected) protonDetectedEvents++;
        if (lithiumDetected) lithiumDetectedEvents++;
        if (alphaDetected) alphaDetectedEvents++;
        if (protonDetected && (lithiumDetected || alphaDetected)) coincidenceEvents++;
    }

    // Cálculo de eficiencias
    double protonEfficiency = static_cast<double>(protonDetectedEvents) / totalEvents * 100;
    double lithiumEfficiency = static_cast<double>(lithiumDetectedEvents) / totalEvents * 100;
    double alphaEfficiency = static_cast<double>(alphaDetectedEvents) / totalEvents * 100;
    double coincidenceEfficiency = static_cast<double>(coincidenceEvents) / totalEvents * 100;

    // Imprimir resultados
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Total eventos procesados: " << totalEvents << std::endl;
    std::cout << "Eventos con detección de protones: " << protonDetectedEvents 
              << " (" << protonEfficiency << "%)" << std::endl;
    std::cout << "Eventos con detección de 7Li: " << lithiumDetectedEvents 
              << " (" << lithiumEfficiency << "%)" << std::endl;
    std::cout << "Eventos con detección de partículas alfa: " << alphaDetectedEvents 
              << " (" << alphaEfficiency << "%)" << std::endl;
    std::cout << "Eventos en coincidencia (protón y (7Li o alfa)): " << coincidenceEvents 
              << " (" << coincidenceEfficiency << "%)" << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    // Dibujar los plots PID combinados
    TCanvas* c_combined = new TCanvas("c_combined", "Combined PID Plots", 800, 800);
    c_combined->Divide(2, 2);
    for (int i = 0; i < 4; i++) {
        c_combined->cd(i+1);
        hPID_Combined[i]->Draw("colz");
        gPad->SetLogz(); // Escala logarítmica para mejor visualización
    }
}


