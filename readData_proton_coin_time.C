Float_t Ex = 0;
Int_t detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x[24];
Float_t thetaCM = 0;

void readData_proton_coin_time()
{
TFile *f = new TFile("h082_10BDP_trace_run013_015-019_025-041.root");
TTree *tree = (TTree*)f->Get("tree");

tree->SetBranchAddress("Ex",&Ex);
tree->SetBranchAddress("e", e);
tree->SetBranchAddress("rdt", rdt);
tree->SetBranchAddress("detID", &detID);
tree->SetBranchAddress("coinTime", &coinTime);
tree->SetBranchAddress("x", x);
tree->SetBranchAddress("thetaCM", &thetaCM);

tree->SetBranchStatus("*", 0);
tree->SetBranchStatus("Ex", 1);
tree->SetBranchStatus("x", 1);
tree->SetBranchStatus("z", 1);
tree->SetBranchStatus("rdt", 1);
tree->SetBranchStatus("thetaCM", 1);
tree->SetBranchStatus("coinTime", 1);
tree->SetBranchStatus("detID", 1);
tree->SetBranchStatus("e", 1);
   
// Cuts blob proton de la izq.
auto cutProtonRecoil1 = new TCutG("CUTPROTONRECOIL1",6);
    cutProtonRecoil1->SetVarX("rdtH[0]");
    cutProtonRecoil1->SetVarY("");
    cutProtonRecoil1->SetTitle("Graph");
    cutProtonRecoil1->SetFillStyle(1000);
    cutProtonRecoil1->SetPoint(0,3220.905,1444.03);
    cutProtonRecoil1->SetPoint(1,3311.691,1170.896);
    cutProtonRecoil1->SetPoint(2,3262.632,1016.418);
    cutProtonRecoil1->SetPoint(3,3164.516,1311.94);
    cutProtonRecoil1->SetPoint(4,3219.777,1448.507);
    cutProtonRecoil1->SetPoint(5,3220.905,1444.03);

auto cutProtonRecoil2 = new TCutG("CUTPROTONRECOIL2",6);
   cutProtonRecoil2->SetVarX("rdtH[1]");
   cutProtonRecoil2->SetVarY("");
   cutProtonRecoil2->SetTitle("Graph");
   cutProtonRecoil2->SetFillStyle(1000);
   cutProtonRecoil2->SetPoint(0,3119.599,1380.306);
   cutProtonRecoil2->SetPoint(1,3046.031,1344.799);
   cutProtonRecoil2->SetPoint(2,3112.911,1077.659);
   cutProtonRecoil2->SetPoint(3,3197.625,1119.928);
   cutProtonRecoil2->SetPoint(4,3120.713,1381.996);
   cutProtonRecoil2->SetPoint(5,3119.599,1380.306);

auto cutProtonRecoil3 = new TCutG("CUTPROTONRECOIL3",6);
   cutProtonRecoil3->SetVarX("rdtH[2]");
   cutProtonRecoil3->SetVarY("");
   cutProtonRecoil3->SetTitle("Graph");
   cutProtonRecoil3->SetFillStyle(1000);
   cutProtonRecoil3->SetPoint(0,3261.741,1338.105);
   cutProtonRecoil3->SetPoint(1,3168.548,1306.482);
   cutProtonRecoil3->SetPoint(2,3249.519,1100.931);
   cutProtonRecoil3->SetPoint(3,3348.823,1144.413);
   cutProtonRecoil3->SetPoint(4,3264.797,1342.058);
   cutProtonRecoil3->SetPoint(5,3261.741,1338.105);

auto cutProtonRecoil4 = new TCutG("CUTPROTONRECOIL4",6);
   cutProtonRecoil4->SetVarX("rdtH[3]");
   cutProtonRecoil4->SetVarY("");
   cutProtonRecoil4->SetTitle("Graph");
   cutProtonRecoil4->SetFillStyle(1000);
   cutProtonRecoil4->SetPoint(0,3196.594,1393.424);
   cutProtonRecoil4->SetPoint(1,3121.597,1296.035);
   cutProtonRecoil4->SetPoint(2,3229.732,1090.438);
   cutProtonRecoil4->SetPoint(3,3311.706,1178.358);
   cutProtonRecoil4->SetPoint(4,3194.85,1389.366);
   cutProtonRecoil4->SetPoint(5,3196.594,1393.424);

   auto cutBoronRecoil1 = new TCutG("cutBoronRecoil1",8);
   cutBoronRecoil1->SetVarX("rdtH[0]");
   cutBoronRecoil1->SetVarY("");
   cutBoronRecoil1->SetTitle("Graph");
   cutBoronRecoil1->SetFillStyle(1000);
   cutBoronRecoil1->SetPoint(0,65.059,4708.98);
   cutBoronRecoil1->SetPoint(1,410.621,3944.01);
   cutBoronRecoil1->SetPoint(2,2157.11,2218.75);
   cutBoronRecoil1->SetPoint(3,3333.89,1795.57);
   cutBoronRecoil1->SetPoint(4,2913.61,2414.06);
   cutBoronRecoil1->SetPoint(5,46.38,5083.33);
   cutBoronRecoil1->SetPoint(6,93.0775,4660.16);
   cutBoronRecoil1->SetPoint(7,65.059,4708.98);

auto cutBoronRecoil2 = new TCutG("CUTBORONRECOIL2",8);
   cutBoronRecoil2->SetVarX("rdtH[1]");
   cutBoronRecoil2->SetVarY("");
   cutBoronRecoil2->SetTitle("Graph");
   cutBoronRecoil2->SetFillStyle(1000);
   cutBoronRecoil2->SetPoint(0,32.3707,4660.16);
   cutBoronRecoil2->SetPoint(1,1302.54,2739.58);
   cutBoronRecoil2->SetPoint(2,3067.71,1567.71);
   cutBoronRecoil2->SetPoint(3,3086.39,1990.89);
   cutBoronRecoil2->SetPoint(4,2768.85,2479.17);
   cutBoronRecoil2->SetPoint(5,23.0312,5115.89);
   cutBoronRecoil2->SetPoint(6,41.7102,4578.78);
   cutBoronRecoil2->SetPoint(7,32.3707,4660.16);

auto cutBoronRecoil3 = new TCutG("CUTBORONRECOIL3",9);
   cutBoronRecoil3->SetVarX("rdtH[2]");
   cutBoronRecoil3->SetVarY("");
   cutBoronRecoil3->SetTitle("Graph");
   cutBoronRecoil3->SetFillStyle(1000);
   cutBoronRecoil3->SetPoint(0,74.3985,4529.95);
   cutBoronRecoil3->SetPoint(1,2054.37,2121.09);
   cutBoronRecoil3->SetPoint(2,3361.91,1567.71);
   cutBoronRecoil3->SetPoint(3,3305.87,1909.51);
   cutBoronRecoil3->SetPoint(4,2866.91,2446.61);
   cutBoronRecoil3->SetPoint(5,149.115,4936.85);
   cutBoronRecoil3->SetPoint(6,93.0775,4481.12);
   cutBoronRecoil3->SetPoint(7,93.0775,4481.12);
   cutBoronRecoil3->SetPoint(8,74.3985,4529.95);

auto cutBoronRecoil4 = new TCutG("CUTBORONRECOIL4",9);
   cutBoronRecoil4->SetVarX("rdtH[3]");
   cutBoronRecoil4->SetVarY("");
   cutBoronRecoil4->SetTitle("Graph");
   cutBoronRecoil4->SetFillStyle(1000);
   cutBoronRecoil4->SetPoint(0,41.7102,4692.71);
   cutBoronRecoil4->SetPoint(1,1106.41,3097.66);
   cutBoronRecoil4->SetPoint(2,2563.38,1828.12);
   cutBoronRecoil4->SetPoint(3,3329.22,1665.36);
   cutBoronRecoil4->SetPoint(4,3198.46,1974.61);
   cutBoronRecoil4->SetPoint(5,2815.54,2511.72);
   cutBoronRecoil4->SetPoint(6,41.7102,5115.89);
   cutBoronRecoil4->SetPoint(7,41.7102,4660.16);
   cutBoronRecoil4->SetPoint(8,41.7102,4692.71);

   auto cutLitium1 = new TCutG("CUTLITIUM1",10);
   cutLitium1->SetVarX("rdtH[0]");
   cutLitium1->SetVarY("");
   cutLitium1->SetTitle("Graph");
   cutLitium1->SetFillStyle(1000);
   cutLitium1->SetPoint(0,82.5365,3299.11);
   cutLitium1->SetPoint(1,386.778,1441.96);
   cutLitium1->SetPoint(2,1321.24,584.821);
   cutLitium1->SetPoint(3,3461.8,263.393);
   cutLitium1->SetPoint(4,4037.68,424.107);
   cutLitium1->SetPoint(5,3374.87,888.393);
   cutLitium1->SetPoint(6,1245.18,1316.96);
   cutLitium1->SetPoint(7,93.4023,3299.11);
   cutLitium1->SetPoint(8,82.5365,3227.68);
   cutLitium1->SetPoint(9,82.5365,3299.11);

auto cutLitium3 = new TCutG("CUTLITIUM3",8);
   cutLitium3->SetVarX("rdtH[2]");
   cutLitium3->SetVarY("");
   cutLitium3->SetTitle("Graph");
   cutLitium3->SetFillStyle(1000);
   cutLitium3->SetPoint(0,267.255,2129.46);
   cutLitium3->SetPoint(1,3689.98,540.179);
   cutLitium3->SetPoint(2,2798.98,236.607);
   cutLitium3->SetPoint(3,1038.73,629.464);
   cutLitium3->SetPoint(4,256.389,2165.18);
   cutLitium3->SetPoint(5,441.107,2004.46);
   cutLitium3->SetPoint(6,419.376,2040.18);
   cutLitium3->SetPoint(7,267.255,2129.46);

auto cutLitium2 = new TCutG("CUTLITIUM2",7);
   cutLitium2->SetVarX("rdtH[1]");
   cutLitium2->SetVarY("");
   cutLitium2->SetTitle("Graph");
   cutLitium2->SetFillStyle(1000);
   cutLitium2->SetPoint(0,202.06,2084.82);
   cutLitium2->SetPoint(1,3461.8,531.25);
   cutLitium2->SetPoint(2,2549.07,334.821);
   cutLitium2->SetPoint(3,930.068,727.679);
   cutLitium2->SetPoint(4,136.865,1977.68);
   cutLitium2->SetPoint(5,245.523,2102.68);
   cutLitium2->SetPoint(6,202.06,2084.82);

auto cutLitium4 = new TCutG("CUTLITIUM4",6);
   cutLitium4->SetVarX("rdtH[3]");
   cutLitium4->SetVarY("");
   cutLitium4->SetTitle("Graph");
   cutLitium4->SetFillStyle(1000);
   cutLitium4->SetPoint(0,147.731,2040.18);
   cutLitium4->SetPoint(1,3776.9,575.893);
   cutLitium4->SetPoint(2,3092.36,272.321);
   cutLitium4->SetPoint(3,886.605,736.607);
   cutLitium4->SetPoint(4,158.597,2093.75);
   cutLitium4->SetPoint(5,147.731,2040.18);

auto cutalpha1 = new TCutG("CUTALPHA1",7);
   cutalpha1->SetVarX("rdtH[0]");
   cutalpha1->SetVarY("");
   cutalpha1->SetTitle("Graph");
   cutalpha1->SetFillStyle(1000);
   cutalpha1->SetPoint(0,39.0733,1120.54);
   cutalpha1->SetPoint(1,39.0733,102.679);
   cutalpha1->SetPoint(2,2896.78,138.393);
   cutalpha1->SetPoint(3,1527.69,513.393);
   cutalpha1->SetPoint(4,28.2075,1120.54);
   cutalpha1->SetPoint(5,39.0733,1120.54);
   cutalpha1->SetPoint(6,39.0733,1120.54);

auto cutalpha3 = new TCutG("CUTALPHA3",8);
   cutalpha3->SetVarX("rdtH[2]");
   cutalpha3->SetVarY("");
   cutalpha3->SetTitle("Graph");
   cutalpha3->SetFillStyle(1000);
   cutalpha3->SetPoint(0,2766.39,254.464);
   cutalpha3->SetPoint(1,701.886,736.607);
   cutalpha3->SetPoint(2,104.268,1522.32);
   cutalpha3->SetPoint(3,17.3418,1058.04);
   cutalpha3->SetPoint(4,202.06,165.179);
   cutalpha3->SetPoint(5,2733.79,129.464);
   cutalpha3->SetPoint(6,2733.79,290.179);
   cutalpha3->SetPoint(7,2766.39,254.464);

auto cutalpha2 = new TCutG("CUTALPHA2",8);
   cutalpha2->SetVarX("rdtH[1]");
   cutalpha2->SetVarY("");
   cutalpha2->SetTitle("Graph");
   cutalpha2->SetFillStyle(1000);
   cutalpha2->SetPoint(0,115.134,1227.68);
   cutalpha2->SetPoint(1,1799.33,424.107);
   cutalpha2->SetPoint(2,2657.73,245.536);
   cutalpha2->SetPoint(3,2049.24,120.536);
   cutalpha2->SetPoint(4,343.315,227.679);
   cutalpha2->SetPoint(5,28.2075,1138.39);
   cutalpha2->SetPoint(6,158.597,1227.68);
   cutalpha2->SetPoint(7,115.134,1227.68);

auto cutalpha4 = new TCutG("CUTALPHA4",7);
   cutalpha4->SetVarX("rdtH[3]");
   cutalpha4->SetVarY("");
   cutalpha4->SetTitle("Graph");
   cutalpha4->SetFillStyle(1000);
   cutalpha4->SetPoint(0,49.9391,1111.61);
   cutalpha4->SetPoint(1,332.45,200.893);
   cutalpha4->SetPoint(2,2538.2,58.0357);
   cutalpha4->SetPoint(3,2646.86,308.036);
   cutalpha4->SetPoint(4,538.899,772.321);
   cutalpha4->SetPoint(5,39.0733,1111.61);
   cutalpha4->SetPoint(6,49.9391,1111.61);


//Histograms
TH1F* coinTimeH = new TH1F("coinTimeH","coinTimeH",1000,-1000,1000);

//TH2F* xh[24];

TH2F* coinTimeVsX[24];
for (int i = 0; i < 24; ++i) {
    coinTimeVsX[i] = new TH2F(Form("coinTimeVsX_%d", i), 
                              Form("Coincidence Time vs X - Detector eH[%d]", i), 
                              100, -2, 2, 100, -100, 100);
                              coinTimeVsX[i]->GetXaxis()->SetTitle("X (Longitudinal Position)");
        coinTimeVsX[i]->GetYaxis()->SetTitle("Coincidence Time (ns)");
}

TH1F* eH[24];
for (auto i = 0; i < 24; ++i)
      eH[i] = new TH1F(Form("eH[%i]", i), Form("eH[%i]", i), 1000, -2, 18); 

TH2F* rdtH[4];
for (auto i = 0; i < 4; ++i)
      rdtH[i] = new TH2F(Form("rdtH[%i]", i), Form("rdtH[%i]", i), 1000, 0, 6000, 1000, 0, 6000);

TH1F* exH[24];
for (auto i = 0; i < 24; i++)
      exH[i] = new TH1F(Form("exH[%i]", i), Form("exH[%i]", i), 1000, -2, 18); 

TH1F* exTotalH = new TH1F("exTotalH","exTotalH",200,0,20);  
/*
TH2F* coinTimeVsExTotalH = new TH2F("coinTimeVsExTotalH", "Coincidence Time vs Total Excitation Energy", 
   100, 9, 15, 100, -100, 150);
coinTimeVsExTotalH->GetXaxis()->SetTitle("ExTotalH (MeV)");
coinTimeVsExTotalH->GetYaxis()->SetTitle("Coincidence Time (ns)");
*/

 Long64_t nentries = tree->GetEntries();
 std::cout<<" Number of entries : "<<nentries<<"\n";
   for (Long64_t i=0;i<nentries;i++) {
     tree->GetEntry(i);
/*
    if(i % 1000000 == 0){
         std::cout<<" Entry number : "<<i<<"\n";
         //if(!std::isnan(e[2])) std::cout<<" Energy index 0 "<<e[2]<<"\n"; 
    }   
*/
     // if(coinTime<-10 || coinTime>20) continue; 

      coinTimeH->Fill(coinTime);

      if (x[detID] < -0.95 || x[detID] > 0.95 || thetaCM < 10 || e[detID] < 1)
       continue;

      // if (!cutProtonRecoil1->IsInside(rdt[0],rdt[1])  && !cutProtonRecoil2->IsInside(rdt[2],rdt[3]) && !cutProtonRecoil3->IsInside(rdt[4],rdt[5]) && !cutProtonRecoil4->IsInside(rdt[6],rdt[7]))
      //continue; 

     // if (!cutBoronRecoil1->IsInside(rdt[0],rdt[1])  && !cutBoronRecoil2->IsInside(rdt[2],rdt[3]) && !cutBoronRecoil3->IsInside(rdt[4],rdt[5]) && !cutBoronRecoil4->IsInside(rdt[6],rdt[7]))
      //continue; 

     // bool passLitiumCut = cutLitium1->IsInside(rdt[0], rdt[1]) || cutLitium2->IsInside(rdt[2], rdt[3]) || cutLitium3->IsInside(rdt[4], rdt[5]) || cutLitium4->IsInside(rdt[6], rdt[7]);
     // bool passAlphaCut = cutalpha1->IsInside(rdt[0], rdt[1]) || cutalpha2->IsInside(rdt[2], rdt[3]) || cutalpha3->IsInside(rdt[4], rdt[5]) || cutalpha4->IsInside(rdt[6], rdt[7]);

       // if (!passLitiumCut && !passAlphaCut)
          //  continue;

      exTotalH->Fill(Ex);
      //Fill the new histogram
    //coinTimeVsExTotalH->Fill(Ex, coinTime);
/*
      for(auto i=0;i<24;++i){
        xh[detID]->Fill(x[detID], coinTime);
      }
      
*/
/*
std::cout << "Event " << i << ": detID = " << detID
          << ", x = " << x[detID] << ", coinTime = " << coinTime << std::endl;

          */
if (detID >= 0 && detID < 24) {  // Asegurar que detID es válido
    coinTimeVsX[detID]->Fill(x[detID], coinTime);
}        
     for(auto i=0;i<24;++i){
        eH[i]->Fill(e[i]);
        if(detID==i) exH[i]->Fill(Ex);
    }
/*
    Float_t totalEnergy_e = 0;
    for (int j = 0; j < 24; j++) {
        if(detID==j) totalEnergy_e += Ex;
    }

    Float_t energyDifference; // Declare energyDifference outside the loop

    for(auto i=0; i<4; ++i) {
        rdtH[i]->Fill(rdt[i*2], rdt[i*2+1]);
        Float_t pairRdtEnergy = (rdt[i * 2] + rdt[i * 2 + 1])/1000;
        energyDifference = 100.0 - totalEnergy_e - pairRdtEnergy; // Calculate it inside the loop
        energyDifferenceH->Fill(energyDifference); // Fill the histogram inside the loop
    }

*/
        // Crear un nuevo canvas para cada rdtH

  }//events

/*
// Crear un nuevo canvas para todos los rdtH
TCanvas *c3 = new TCanvas("c3", "PIDs with Cuts", 800, 600);
c3->Divide(2,2);

for (int i = 0; i < 4; ++i) {
    c3->cd(i+1);
    
    // Dibujar el histograma
    rdtH[i]->Draw("colz");
    
    rdtH[i]->GetXaxis()->SetTitle("E_{res} (keV)");
    rdtH[i]->GetYaxis()->SetTitle("DE (keV)");
    rdtH[i]->SetStats(0);
    
    // Dibujar los cuts correspondientes
    switch(i) {
        case 0:
            cutProtonRecoil1->Draw("same");
            break;
        case 1:
            cutProtonRecoil2->Draw("same");
            break;
        case 2:
            cutProtonRecoil3->Draw("same");
            break;
        case 3:
            cutProtonRecoil4->Draw("same");
            break;
    }
    
    // Actualizar el pad
    gPad->Update();
}

// Actualizar el canvas completo
c3->Update();
*/
gROOT->SetBatch(kFALSE);
/*
TCanvas *c4 = new TCanvas();
coinTimeH->Draw();

//Draw the new histogram
TCanvas *c8 = new TCanvas("c8", "Coincidence Time vs ExTotalH", 800, 600);
coinTimeVsExTotalH->SetMarkerSize(1.5);  // Ajusta este valor según lo necesites
coinTimeVsExTotalH->Draw("scat");

TCanvas *c1 = new TCanvas();
c1->Divide(6,4);
  for(auto i=0;i<24;++i){
    c1->cd(i+1);
    exH[i]->Draw();
  }


TCanvas *c2 = new TCanvas();
c2->Divide(6,4);
  for(auto i=0;i<24;++i){
    c2->cd(i+1);
    eH[i]->Draw();
  }
  */
/*
TCanvas *c3 = new TCanvas();
c3->Divide(2,2);
  for(auto i=0;i<4;++i){
    c3->cd(i+1);
    rdtH[i]->Draw("colz");
  } 
   */ 
  
    // Para hacer los profiles a mano 1 ventana en 1 ventana.
/*
   for (int i = 0; i < 24; ++i) {
      // Crear una nueva ventana (canvas) para cada histograma
      TCanvas* canvas = new TCanvas(Form("canvas_%d", i), Form("Canvas %d", i), 800, 600);
  
      // Dibujar el histograma en la ventana creada
      coinTimeVsX[i]->Draw("colz");
      coinTimeVsX[i]->GetXaxis()->SetTitle("X");
      coinTimeVsX[i]->GetXaxis()->SetRangeUser(-0.95, 0.95);
      coinTimeVsX[i]->GetYaxis()->SetTitle("Coincidence Time (ns)");
      coinTimeVsX[i]->SetTitle(Form("Coincidence Time vs X [%d]", i));
  }
  */

  TCanvas* canvas = new TCanvas("canvas", " CoinTimeVsX", 1800, 1200);
    canvas->Divide(6, 4);

    // Draw the corrected heatmaps
    for (auto detID = 0; detID < 24; detID++) {
      // Crear un canvas grande con 6x4 pads
        canvas->cd(detID + 1);
        coinTimeVsX[detID]->Draw("COLZ");
        coinTimeVsX[detID]->GetXaxis()->SetTitle("X");
        coinTimeVsX[detID]->GetXaxis()->SetRangeUser(-2, 2);
        coinTimeVsX[detID]->GetYaxis()->SetTitle("Coincidence Time (ns)");
        coinTimeVsX[detID]->SetTitle(Form("Coincidence Time vs X [%d]", detID));

    }

   
/*
TCanvas* c9 = new TCanvas("c", "Histograma combinado", 800, 600);
hSum->Draw("COLZ");


// Dibujar el nuevo histograma
TCanvas *c5 = new TCanvas("c5", "Energy Difference", 800, 600);
energyDifferenceH->Draw();
*/
// Mostrar el Canvas
gPad->Update();


}
