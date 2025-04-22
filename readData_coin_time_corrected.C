#include <TCanvas.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TTree.h>
#include <iostream>
#include <cmath>

Float_t Ex = 0;
Int_t detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x[24];
Float_t thetaCM = 0;

void readData_coin_time_corrected() {
    TFile *f = new TFile("h082_10BDP_trace_run013_015-019_025-041.root");
    TTree *tree = (TTree*)f->Get("tree");

    tree->SetBranchAddress("Ex", &Ex);
    tree->SetBranchAddress("e", e);
    tree->SetBranchAddress("rdt", rdt);
    tree->SetBranchAddress("detID", &detID);
    tree->SetBranchAddress("coinTime", &coinTime);
    tree->SetBranchAddress("x", x);
    tree->SetBranchAddress("thetaCM", &thetaCM);

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("Ex", 1);
    tree->SetBranchStatus("x", 1);
    tree->SetBranchStatus("rdt", 1);
    tree->SetBranchStatus("thetaCM", 1);
    tree->SetBranchStatus("coinTime", 1);
    tree->SetBranchStatus("detID", 1);
    tree->SetBranchStatus("e", 1);

    // Cuts
    auto cutProtonRecoil1 = new TCutG("CUTPROTONRECOIL1", 6);
    cutProtonRecoil1->SetPoint(0, 3220.905, 1444.03);
    cutProtonRecoil1->SetPoint(1, 3311.691, 1170.896);
    cutProtonRecoil1->SetPoint(2, 3262.632, 1016.418);
    cutProtonRecoil1->SetPoint(3, 3164.516, 1311.94);
    cutProtonRecoil1->SetPoint(4, 3219.777, 1448.507);
    cutProtonRecoil1->SetPoint(5, 3220.905, 1444.03);

    auto cutProtonRecoil2 = new TCutG("CUTPROTONRECOIL2", 6);
    cutProtonRecoil2->SetPoint(0, 3119.599, 1380.306);
    cutProtonRecoil2->SetPoint(1, 3046.031, 1344.799);
    cutProtonRecoil2->SetPoint(2, 3112.911, 1077.659);
    cutProtonRecoil2->SetPoint(3, 3197.625, 1119.928);
    cutProtonRecoil2->SetPoint(4, 3120.713, 1381.996);
    cutProtonRecoil2->SetPoint(5, 3119.599, 1380.306);

    auto cutProtonRecoil3 = new TCutG("CUTPROTONRECOIL3", 6);
    cutProtonRecoil3->SetPoint(0, 3261.741, 1338.105);
    cutProtonRecoil3->SetPoint(1, 3168.548, 1306.482);
    cutProtonRecoil3->SetPoint(2, 3249.519, 1100.931);
    cutProtonRecoil3->SetPoint(3, 3348.823, 1144.413);
    cutProtonRecoil3->SetPoint(4, 3264.797, 1342.058);
    cutProtonRecoil3->SetPoint(5, 3261.741, 1338.105);

    auto cutProtonRecoil4 = new TCutG("CUTPROTONRECOIL4", 6);
    cutProtonRecoil4->SetPoint(0, 3196.594, 1393.424);
    cutProtonRecoil4->SetPoint(1, 3121.597, 1296.035);
    cutProtonRecoil4->SetPoint(2, 3229.732, 1090.438);
    cutProtonRecoil4->SetPoint(3, 3311.706, 1178.358);
    cutProtonRecoil4->SetPoint(4, 3194.85, 1389.366);
    cutProtonRecoil4->SetPoint(5, 3196.594, 1393.424);

    auto cutBoronRecoil1 = new TCutG("cutBoronRecoil1", 8);
    cutBoronRecoil1->SetPoint(0, 65.059, 4708.98);
    cutBoronRecoil1->SetPoint(1, 410.621, 3944.01);
    cutBoronRecoil1->SetPoint(2, 2157.11, 2218.75);
    cutBoronRecoil1->SetPoint(3, 3333.89, 1795.57);
    cutBoronRecoil1->SetPoint(4, 2913.61, 2414.06);
    cutBoronRecoil1->SetPoint(5, 46.38, 5083.33);
    cutBoronRecoil1->SetPoint(6, 93.0775, 4660.16);
    cutBoronRecoil1->SetPoint(7, 65.059, 4708.98);

    auto cutBoronRecoil2 = new TCutG("CUTBORONRECOIL2", 8);
    cutBoronRecoil2->SetPoint(0, 32.3707, 4660.16);
    cutBoronRecoil2->SetPoint(1, 1302.54, 2739.58);
    cutBoronRecoil2->SetPoint(2, 3067.71, 1567.71);
    cutBoronRecoil2->SetPoint(3, 3086.39, 1990.89);
    cutBoronRecoil2->SetPoint(4, 2768.85, 2479.17);
    cutBoronRecoil2->SetPoint(5, 23.0312, 5115.89);
    cutBoronRecoil2->SetPoint(6, 41.7102, 4578.78);
    cutBoronRecoil2->SetPoint(7, 32.3707, 4660.16);

    auto cutBoronRecoil3 = new TCutG("CUTBORONRECOIL3", 9);
    cutBoronRecoil3->SetPoint(0, 74.3985, 4529.95);
    cutBoronRecoil3->SetPoint(1, 2054.37, 2121.09);
    cutBoronRecoil3->SetPoint(2, 3361.91, 1567.71);
    cutBoronRecoil3->SetPoint(3, 3305.87, 1909.51);
    cutBoronRecoil3->SetPoint(4, 2866.91, 2446.61);
    cutBoronRecoil3->SetPoint(5, 149.115, 4936.85);
    cutBoronRecoil3->SetPoint(6, 93.0775, 4481.12);
    cutBoronRecoil3->SetPoint(7, 93.0775, 4481.12);
    cutBoronRecoil3->SetPoint(8, 74.3985, 4529.95);

    auto cutBoronRecoil4 = new TCutG("CUTBORONRECOIL4", 9);
    cutBoronRecoil4->SetPoint(0, 41.7102, 4692.71);
    cutBoronRecoil4->SetPoint(1, 1106.41, 3097.66);
    cutBoronRecoil4->SetPoint(2, 2563.38, 1828.12);
    cutBoronRecoil4->SetPoint(3, 3329.22, 1665.36);
    cutBoronRecoil4->SetPoint(4, 3198.46, 1974.61);
    cutBoronRecoil4->SetPoint(5, 2815.54, 2511.72);
    cutBoronRecoil4->SetPoint(6, 41.7102, 5115.89);
    cutBoronRecoil4->SetPoint(7, 41.7102, 4660.16);
    cutBoronRecoil4->SetPoint(8, 41.7102, 4692.71);

    auto cutLitium1 = new TCutG("CUTLITIUM1", 10);
    cutLitium1->SetPoint(0, 82.5365, 3299.11);
    cutLitium1->SetPoint(1, 386.778, 1441.96);
    cutLitium1->SetPoint(2, 1321.24, 584.821);
    cutLitium1->SetPoint(3, 3461.8, 263.393);
    cutLitium1->SetPoint(4, 4037.68, 424.107);
    cutLitium1->SetPoint(5, 3374.87, 888.393);
    cutLitium1->SetPoint(6, 1245.18, 1316.96);
    cutLitium1->SetPoint(7, 93.4023, 3299.11);
    cutLitium1->SetPoint(8, 82.5365, 3227.68);
    cutLitium1->SetPoint(9, 82.5365, 3299.11);

    auto cutLitium3 = new TCutG("CUTLITIUM3", 8);
    cutLitium3->SetPoint(0, 267.255, 2129.46);
    cutLitium3->SetPoint(1, 3689.98, 540.179);
    cutLitium3->SetPoint(2, 2798.98, 236.607);
    cutLitium3->SetPoint(3, 1038.73, 629.464);
    cutLitium3->SetPoint(4, 256.389, 2165.18);
    cutLitium3->SetPoint(5, 441.107, 2004.46);
    cutLitium3->SetPoint(6, 419.376, 2040.18);
    cutLitium3->SetPoint(7, 267.255, 2129.46);

    auto cutLitium2 = new TCutG("CUTLITIUM2", 7);
    cutLitium2->SetPoint(0, 202.06, 2084.82);
    cutLitium2->SetPoint(1, 3461.8, 531.25);
    cutLitium2->SetPoint(2, 2549.07, 334.821);
    cutLitium2->SetPoint(3, 930.068, 727.679);
    cutLitium2->SetPoint(4, 136.865, 1977.68);
    cutLitium2->SetPoint(5, 245.523, 2102.68);
    cutLitium2->SetPoint(6, 202.06, 2084.82);

    auto cutLitium4 = new TCutG("CUTLITIUM4", 6);
    cutLitium4->SetPoint(0, 147.731, 2040.18);
    cutLitium4->SetPoint(1, 3776.9, 575.893);
    cutLitium4->SetPoint(2, 3092.36, 272.321);
    cutLitium4->SetPoint(3, 886.605, 736.607);
    cutLitium4->SetPoint(4, 158.597, 2093.75);
    cutLitium4->SetPoint(5, 147.731, 2040.18);

    auto cutalpha1 = new TCutG("CUTALPHA1", 7);
    cutalpha1->SetPoint(0, 39.0733, 1120.54);
    cutalpha1->SetPoint(1, 39.0733, 102.679);
    cutalpha1->SetPoint(2, 2896.78, 138.393);
    cutalpha1->SetPoint(3, 1527.69, 513.393);
    cutalpha1->SetPoint(4, 28.2075, 1120.54);
    cutalpha1->SetPoint(5, 39.0733, 1120.54);
    cutalpha1->SetPoint(6, 39.0733, 1120.54);

    auto cutalpha3 = new TCutG("CUTALPHA3", 8);
    cutalpha3->SetPoint(0, 2766.39, 254.464);
    cutalpha3->SetPoint(1, 701.886, 736.607);
    cutalpha3->SetPoint(2, 104.268, 1522.32);
    cutalpha3->SetPoint(3, 17.3418, 1058.04);
    cutalpha3->SetPoint(4, 202.06, 165.179);
    cutalpha3->SetPoint(5, 2733.79, 129.464);
    cutalpha3->SetPoint(6, 2733.79, 290.179);
    cutalpha3->SetPoint(7, 2766.39, 254.464);

    auto cutalpha2 = new TCutG("CUTALPHA2", 8);
    cutalpha2->SetPoint(0, 115.134, 1227.68);
    cutalpha2->SetPoint(1, 1799.33, 424.107);
    cutalpha2->SetPoint(2, 2657.73, 245.536);
    cutalpha2->SetPoint(3, 2049.24, 120.536);
    cutalpha2->SetPoint(4, 343.315, 227.679);
    cutalpha2->SetPoint(5, 28.2075, 1138.39);
    cutalpha2->SetPoint(6, 158.597, 1227.68);
    cutalpha2->SetPoint(7, 115.134, 1227.68);

    auto cutalpha4 = new TCutG("CUTALPHA4", 7);
    cutalpha4->SetPoint(0, 49.9391, 1111.61);
    cutalpha4->SetPoint(1, 332.45, 200.893);
    cutalpha4->SetPoint(2, 2538.2, 58.0357);
    cutalpha4->SetPoint(3, 2646.86, 308.036);
    cutalpha4->SetPoint(4, 538.899, 772.321);
    cutalpha4->SetPoint(5, 39.0733, 1111.61);
    cutalpha4->SetPoint(6, 49.9391, 1111.61);

    TH2F* coinTimeVsX[24];
    for (int i = 0; i < 24; ++i) {
        coinTimeVsX[i] = new TH2F(Form("coinTimeVsX_%d", i),
                                  Form("Coincidence Time vs X - Detector eH[%d]", i),
                                  10000, -2, 2, 10000, -100, 100);
        coinTimeVsX[i]->GetXaxis()->SetTitle("X (Longitudinal Position)");
        coinTimeVsX[i]->GetYaxis()->SetTitle("Coincidence Time (ns)");
    }

    Long64_t nentries = tree->GetEntries();
    std::cout << "Number of entries : " << nentries << "\n";

    // Coeficientes polinómicos para la corrección (pol4)
    float pol4[24][5] = {
        {22.2632, -0.443764, -95.8043, 1.28447, 80.6196},
        {25.6633, 0.939019, -86.6471, -1.27647, 70.5687},
        {39.4859, -1.42573, -115.269, 2.50988, 89.4895},
        {26.2832, -4.68976, -88.0606, 12.0906, 71.8652},
        {30.0488, 4.98238, -80.8835, -8.06761, 59.2426},
        {36.7207, 3.9082, -95.6975, -5.3619, 71.3784},
        {20.5705, 3.91639, -27.3175, -4.67329, 15.0882},
        {25.4729, 2.61548, -32.0358, -0.697981, 17.4898},
        {20.7718, 3.45079, -30.4761, -5.86632, 20.0606},
        {47.5046, 2.57078, -158.759, -6.83235, 135.341},
        {40.6197, -0.42643, -96.5613, 0.173386, 72.1999},
        {0, 0, 0, 0, 0},
        {21.3337, 3.09734, -31.0887, -4.00537, 17.6514},
        {23.8342, 0.712722, -24.1388, 0.185379, 7.58745},
        {18.5734, 2.8507, -23.9463, -4.71869, 13.1247},
        {25.1945, 3.23653, -27.8029, -4.67166, 12.2644},
        {24.0988, -7.92557, -74.3601, 10.7257, 59.0229},
        {30.4273, -1.36079, -113.128, 1.84176, 100.789},
        {19.3952, 4.53279, -71.1824, -6.23276, 58.6395},
        {32.0866, -0.850581, -104.966, 3.29582, 82.7438},
        {26.6236, 6.76763, -79.3186, -10.5532, 64.7973},
        {26.5402, 4.55678, -51.5956, -5.42489, 35.7765},
        {0, 0, 0, 0, 0},
        {28.9171, -1.07817, -82.8819, 2.867, 74.3833}
    };

    TH2F* correctedHeatmap[24];
      for (int i = 0; i < 24; ++i) {
        correctedHeatmap[i] = new TH2F(Form("correctedHeatmap_%d", i),
                                  Form("Corrected Coincidence Time vs X - Detector eH[%d]", i),
                                  1000, -2, 2, 1000, -100, 100);
        correctedHeatmap[i]->GetXaxis()->SetTitle("X (Longitudinal Position)");
        correctedHeatmap[i]->GetYaxis()->SetTitle("Coincidence Time Corrected (ns)");
    }

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        if (x[detID] < -0.95 || x[detID] > 0.95 || thetaCM < 10 || e[detID] < 1)
            continue;

        //if (!cutProtonRecoil1->IsInside(rdt[0], rdt[1]) && !cutProtonRecoil2->IsInside(rdt[2], rdt[3]) &&            !cutProtonRecoil3->IsInside(rdt[4], rdt[5]) && !cutProtonRecoil4->IsInside(rdt[6], rdt[7]))
           // continue;

        // Calculate polynomial correction
        double poly_correction = pol4[detID][0]
                                + pol4[detID][1] * x[detID]
                                + pol4[detID][2] * x[detID] * x[detID]
                                + pol4[detID][3] * x[detID] * x[detID] * x[detID]
                                + pol4[detID][4] * x[detID] * x[detID] * x[detID] * x[detID];

        // Apply the correction to coinTime
        double correctedCoinTime = coinTime - poly_correction;

        // Fill the corrected heatmap
        correctedHeatmap[detID]->Fill(x[detID], correctedCoinTime);
    }

    TCanvas* canvas = new TCanvas("canvas", "Corrected Heatmaps", 1800, 1200);
    canvas->Divide(6, 4);

    // Draw the corrected heatmaps
    for (auto detID = 0; detID < 24; detID++) {
      // Crear un canvas grande con 6x4 pads
        canvas->cd(detID + 1);
        correctedHeatmap[detID]->Draw("COLZ");
        correctedHeatmap[detID]->GetXaxis()->SetTitle("X");
        correctedHeatmap[detID]->GetXaxis()->SetRangeUser(-2, 2);
        correctedHeatmap[detID]->GetYaxis()->SetTitle("Coincidence Time (ns)");
        correctedHeatmap[detID]->SetTitle(Form("Coincidence Time vs X [%d]", detID));

    }
        

}





