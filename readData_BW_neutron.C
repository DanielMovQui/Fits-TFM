Float_t Ex = 0;
Int_t detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x = 0;
Float_t thetaCM = 0;


void readData_BW_neutron()
{
TFile *f = new TFile("h082_10BDP_trace_run013_015-019_025-041.root");
TTree *tree = (TTree*)f->Get("tree");

tree->SetBranchAddress("Ex",&Ex);
tree->SetBranchAddress("e", e);
tree->SetBranchAddress("rdt", rdt);
tree->SetBranchAddress("detID", &detID);
tree->SetBranchAddress("coinTime", &coinTime);
tree->SetBranchAddress("x", &x);
tree->SetBranchAddress("thetaCM", &thetaCM);

tree->SetBranchStatus("*", 0);
tree->SetBranchStatus("Ex", 1);
tree->SetBranchStatus("x", 1);
tree->SetBranchStatus("rdt", 1);
tree->SetBranchStatus("thetaCM", 1);
tree->SetBranchStatus("coinTime", 1);
tree->SetBranchStatus("detID", 1);
tree->SetBranchStatus("e", 1);

//Cuts
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
   cutProtonRecoil2->Draw("l");

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

auto cutBoronRecoil1 = new TCutG("CUTBORONRECOIL1",19);
   cutBoronRecoil1->SetVarX("rdtH[0]");
   cutBoronRecoil1->SetVarY("");
   cutBoronRecoil1->SetTitle("Graph");
   cutBoronRecoil1->SetFillStyle(1000);
   cutBoronRecoil1->SetPoint(0,128.5845,5519.572);
   cutBoronRecoil1->SetPoint(1,437.6618,4705.767);
   cutBoronRecoil1->SetPoint(2,976.4588,3987.703);
   cutBoronRecoil1->SetPoint(3,1749.152,3221.769);
   cutBoronRecoil1->SetPoint(4,2434.134,2663.276);
   cutBoronRecoil1->SetPoint(5,2881.043,2455.835);
   cutBoronRecoil1->SetPoint(6,3014.698,2535.62);
   cutBoronRecoil1->SetPoint(7,3378.073,2415.943);
   cutBoronRecoil1->SetPoint(8,3453.254,2128.718);
   cutBoronRecoil1->SetPoint(9,3373.896,1801.6);
   cutBoronRecoil1->SetPoint(10,3236.065,1761.708);
   cutBoronRecoil1->SetPoint(11,2642.97,2009.04);
   cutBoronRecoil1->SetPoint(12,1711.562,2671.254);
   cutBoronRecoil1->SetPoint(13,851.1572,3397.296);
   cutBoronRecoil1->SetPoint(14,178.7051,4370.67);
   cutBoronRecoil1->SetPoint(15,45.05008,4681.831);
   cutBoronRecoil1->SetPoint(16,15.81304,5288.196);
   cutBoronRecoil1->SetPoint(17,86.81729,5551.485);
   cutBoronRecoil1->SetPoint(18,128.5845,5519.572); 

auto cutBoronRecoil2 = new TCutG("CUTBORONRECOIL2",16);
   cutBoronRecoil2->SetVarX("rdtH[1]");
   cutBoronRecoil2->SetVarY("");
   cutBoronRecoil2->SetTitle("Graph");
   cutBoronRecoil2->SetFillStyle(1000);
   cutBoronRecoil2->SetPoint(0,193.6272,5099.459);
   cutBoronRecoil2->SetPoint(1,737.4532,4212.577);
   cutBoronRecoil2->SetPoint(2,1289.646,3605.343);
   cutBoronRecoil2->SetPoint(3,2004.986,2950.169);
   cutBoronRecoil2->SetPoint(4,2540.446,2614.592);
   cutBoronRecoil2->SetPoint(5,2967.14,2510.723);
   cutBoronRecoil2->SetPoint(6,3301.802,2406.854);
   cutBoronRecoil2->SetPoint(7,3276.703,1951.428);
   cutBoronRecoil2->SetPoint(8,2653.394,1943.438);
   cutBoronRecoil2->SetPoint(9,1661.958,2662.531);
   cutBoronRecoil2->SetPoint(10,616.1382,3637.302);
   cutBoronRecoil2->SetPoint(11,105.7783,4508.204);
   cutBoronRecoil2->SetPoint(12,26.29606,4955.64);
   cutBoronRecoil2->SetPoint(13,201.9937,5115.439);
   cutBoronRecoil2->SetPoint(14,201.9937,5115.439);
   cutBoronRecoil2->SetPoint(15,193.6272,5099.459); 

auto cutBoronRecoil3 = new TCutG("CUTBORONRECOIL3",16);
   cutBoronRecoil3->SetVarX("rdtH[2]");
   cutBoronRecoil3->SetVarY("");
   cutBoronRecoil3->SetTitle("Graph");
   cutBoronRecoil3->SetFillStyle(1000);
   cutBoronRecoil3->SetPoint(0,142.6864,5002.586);
   cutBoronRecoil3->SetPoint(1,544.0122,4525.102);
   cutBoronRecoil3->SetPoint(2,1048.773,3821.015);
   cutBoronRecoil3->SetPoint(3,2045.881,3011.72);
   cutBoronRecoil3->SetPoint(4,2790.609,2550.422);
   cutBoronRecoil3->SetPoint(5,3138.149,2518.05);
   cutBoronRecoil3->SetPoint(6,3452.59,2323.82);
   cutBoronRecoil3->SetPoint(7,3460.865,1789.685);
   cutBoronRecoil3->SetPoint(8,2956.104,1708.756);
   cutBoronRecoil3->SetPoint(9,2401.696,2032.474);
   cutBoronRecoil3->SetPoint(10,345.418,3982.874);
   cutBoronRecoil3->SetPoint(11,97.17527,4427.986);
   cutBoronRecoil3->SetPoint(12,64.07623,4929.749);
   cutBoronRecoil3->SetPoint(13,134.4117,5018.772);
   cutBoronRecoil3->SetPoint(14,134.4117,5018.772);
   cutBoronRecoil3->SetPoint(15,142.6864,5002.586);

auto cutBoronRecoil4 = new TCutG("CUTBORONRECOIL4",17);
   cutBoronRecoil4->SetVarX("rdtH[3]");
   cutBoronRecoil4->SetVarY("");
   cutBoronRecoil4->SetTitle("Graph");
   cutBoronRecoil4->SetFillStyle(1000);
   cutBoronRecoil4->SetPoint(0,180.9049,5017.467);
   cutBoronRecoil4->SetPoint(1,805.8984,4237.362);
   cutBoronRecoil4->SetPoint(2,1377.788,3565.144);
   cutBoronRecoil4->SetPoint(3,2174.348,2876.327);
   cutBoronRecoil4->SetPoint(4,2766.662,2486.275);
   cutBoronRecoil4->SetPoint(5,3085.286,2436.481);
   cutBoronRecoil4->SetPoint(6,3346.721,2411.584);
   cutBoronRecoil4->SetPoint(7,3407.995,2129.418);
   cutBoronRecoil4->SetPoint(8,3277.278,1764.262);
   cutBoronRecoil4->SetPoint(9,2774.832,1855.551);
   cutBoronRecoil4->SetPoint(10,1937.422,2270.501);
   cutBoronRecoil4->SetPoint(11,585.3125,3706.227);
   cutBoronRecoil4->SetPoint(12,160.4803,4353.548);
   cutBoronRecoil4->SetPoint(13,13.42304,4834.89);
   cutBoronRecoil4->SetPoint(14,58.35721,5025.766);
   cutBoronRecoil4->SetPoint(15,184.9899,5025.766);
   cutBoronRecoil4->SetPoint(16,180.9049,5017.467);


//Histograms
TH1F* coinTimeH = new TH1F("coinTimeH","coinTimeH",1000,-1000,1000);

TH1F* eH[24];
for (auto i = 0; i < 24; ++i)
      eH[i] = new TH1F(Form("eH[%i]", i), Form("eH[%i]", i), 1000, -2, 18); 

TH2F* rdtH[4];
for (auto i = 0; i < 4; ++i)
      rdtH[i] = new TH2F(Form("rdtH[%i]", i), Form("rdtH[%i]", i), 1000, 0, 6000, 1000, 0, 6000);

TH1F* exH[24];
for (auto i = 0; i < 24; i++)
      exH[i] = new TH1F(Form("exH[%i]", i), Form("exH[%i]", i), 1000, -2, 18); 

TH1F* exTotalH = new TH1F("exTotalH","exTotalH",300,-2,18);      


 Long64_t nentries = tree->GetEntries();
 std::cout<<" Number of entries : "<<nentries<<"\n";
   for (Long64_t i=0;i<nentries;i++) {
     tree->GetEntry(i);

    if(i % 1000000 == 0){
         std::cout<<" Entry number : "<<i<<"\n";
         //if(!std::isnan(e[2])) std::cout<<" Energy index 0 "<<e[2]<<"\n"; 
    }   

     if(coinTime<-10 || coinTime>20) continue; 

      coinTimeH->Fill(coinTime);

      if (x < -0.95 || x > 0.95 || thetaCM < 10 || e[detID] < 1)
        continue;

    //if (!cutProtonRecoil1->IsInside(rdt[0],rdt[1])  && !cutProtonRecoil2->IsInside(rdt[2],rdt[3]) && !cutProtonRecoil3->IsInside(rdt[4],rdt[5]) && !cutProtonRecoil4->IsInside(rdt[6],rdt[7]))
      //continue; 

     if (!cutBoronRecoil1->IsInside(rdt[0],rdt[1])  && !cutBoronRecoil2->IsInside(rdt[2],rdt[3]) && !cutBoronRecoil3->IsInside(rdt[4],rdt[5]) && !cutBoronRecoil4->IsInside(rdt[6],rdt[7]))
      continue; 
    
      exTotalH->Fill(Ex);

     for(auto i=0;i<24;++i){
        eH[i]->Fill(e[i]);
        if(detID==i) exH[i]->Fill(Ex);
    }

     for(auto i=0;i<4;++i)
        rdtH[i]->Fill(rdt[i*2],rdt[i*2+1]);

  }//events

gROOT->SetBatch(kFALSE);

/*
TCanvas *c1 = new TCanvas();
c1->Divide(2, 4);
int selectedIndices[] = {3, 4, 8, 9, 10, 14, 15, 21};

for (int i = 0; i < sizeof(selectedIndices) / sizeof(selectedIndices[0]); ++i) {
    int index = selectedIndices[i];
    c1->cd(i + 1);
    exH[index]->Draw();
}
*/
/*
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

  TCanvas *c3 = new TCanvas();
  c3->Divide(2,2);
    for(auto i=0;i<4;++i){
      c3->cd(i+1);
      rdtH[i]->Draw("colz");
    } 

  c3->cd(1);
  cutProtonRecoil1->Draw("same"); 
  cutBoronRecoil1->Draw("same"); 
  c3->cd(2);
  cutProtonRecoil2->Draw("same"); 
  cutBoronRecoil2->Draw("same"); 
  c3->cd(3);
  cutProtonRecoil3->Draw("same");
  cutBoronRecoil3->Draw("same"); 
  c3->cd(4);
  cutProtonRecoil4->Draw("same");
  cutBoronRecoil4->Draw("same"); 

  TCanvas *c4 = new TCanvas();
  coinTimeH->Draw();


  TCanvas *c5 = new TCanvas();
  exTotalH->Draw();
*/

TF1 *bw1 = new TF1("m1", "[0] / ((x * x - [1] * [1]) * (x * x - [1] * [1]) + [1] * [1] * [2] * [2])", 11.2, 12.0);
TF1 *bw2 = new TF1("m2", "[3] / ((x * x - [4] * [4]) * (x * x - [4] * [4]) + [4] * [4] * [5] * [5])", 11.5, 12.5);
TF1 *bw3 = new TF1("m3", "[6] / ((x * x - [7] * [7]) * (x * x - [7] * [7]) + [7] * [7] * [8] * [8])", 12.0, 13.0);
TF1 *bw4 = new TF1("m4", "[9] / ((x * x - [10] * [10]) * (x * x - [10] * [10]) + [10] * [10] * [11] * [11])", 12.5, 14);
TF1 *bw5 = new TF1("m5", "[12] / ((x * x - [13] * [13]) * (x * x - [13] * [13]) + [13] * [13] * [14] * [14])", 12.8, 14.4);
TF1 *bw6 = new TF1("m6", "[15] / ((x * x - [16] * [16]) * (x * x - [16] * [16]) + [16] * [16] * [17] * [17])", 13.7, 14.5);

// Definir parámetros iniciales para cada Breit-Wigner
bw1->SetParameters(300, 11.6, 0.180);  // Parámetros iniciales: amplitud, media, anchura (ancho a media altura)
bw2->SetParameters(1000, 11.9, 0.194);
bw3->SetParameters(2500, 12.5, 0.20);
bw4->SetParameters(2500, 13.3, 0.30);
bw5->SetParameters(2500, 13.6, 0.3);
bw6->SetParameters(1000, 14, 0.5);

// Definir total dada por la suma de los 8 picos
TF1 *total = new TF1("mstotal", "[0] / ((x * x - [1] * [1]) * (x * x - [1] * [1]) + [1] * [1] * [2] * [2]) + [3] / ((x * x - [4] * [4]) * (x * x - [4] * [4]) + [4] * [4] * [5] * [5]) + [6] / ((x * x - [7] * [7]) * (x * x - [7] * [7]) + [7] * [7] * [8] * [8]) + [9] / ((x * x - [10] * [10]) * (x * x - [10] * [10]) + [10] * [10] * [11] * [11]) + [12] / ((x * x - [13] * [13]) * (x * x - [13] * [13]) + [13] * [13] * [14] * [14]) + [15] / ((x * x - [16] * [16]) * (x * x - [16] * [16]) + [16] * [16] * [17] * [17])", 11, 15);

// Ajustar cada función a los datos teniendo en cuenta la anterior
exTotalH->Fit(bw1, "R");
exTotalH->Fit(bw2, "R+");
exTotalH->Fit(bw3, "R+");
exTotalH->Fit(bw4, "R+");
exTotalH->Fit(bw5, "R+");
exTotalH->Fit(bw6, "R+");

// Obtener los parámetros del fit
Double_t par[18];
bw1->GetParameters(&par[0]);
bw2->GetParameters(&par[3]);
bw3->GetParameters(&par[6]);
bw4->GetParameters(&par[9]);
bw5->GetParameters(&par[12]);
bw6->GetParameters(&par[15]);

total->SetParameters(par);

// Configurar nombres de parámetros uno por uno
total->SetParName(0, "Amp1");
total->SetParName(1, "Mean1");
total->SetParName(2, "Width1");
total->SetParName(3, "Amp2");
total->SetParName(4, "Mean2");
total->SetParName(5, "Width2");
total->SetParName(6, "Amp3");
total->SetParName(7, "Mean3");
total->SetParName(8, "Width3");
total->SetParName(9, "Amp4");
total->SetParName(10, "Mean4");
total->SetParName(11, "Width4");
total->SetParName(12, "Amp5");
total->SetParName(13, "Mean5");
total->SetParName(14, "Width5");
total->SetParName(15, "Amp6");
total->SetParName(16, "Mean6");
total->SetParName(17, "Width6");

// Establecer límites de parámetros (si es necesario)
total->SetParLimits(1, 11.50, 11.7);
total->SetParLimits(2, 0.1, 0.4);
total->SetParLimits(4, 11.80, 12);
total->SetParLimits(5, 0.1, 0.4);
total->SetParLimits(7, 12.4, 12.6);
total->SetParLimits(8, 0.180, 1);
total->SetParLimits(10, 13.25, 13.35);
total->SetParLimits(11, 0.1, 0.4);
total->SetParLimits(13, 13.55, 13.65);
total->SetParLimits(14, 0.1, 0.4);
total->SetParLimits(16, 14.05, 14.15);
total->SetParLimits(17, 0.1, 0.4);

// Ajustar solo la función total y desactivar la visualización de la línea de ajuste resultante
exTotalH->Fit(total, "R+");

// Configurar opciones de visualización para la línea de ajuste resultante
total->SetLineColor(kRed);
total->SetLineWidth(2);
exTotalH->SetTitle("Energy distribution ^{10}B"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");
// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[18];
total->GetParameters(par_total);

// Crear nuevas funciones Breit-Wigner con los parámetros ajustados
TF1 *new_bw1 = new TF1("new_m1", "[0] / ((x * x - [1] * [1]) * (x * x - [1] * [1]) + [1] * [1] * [2] * [2])", 10, 13);
TF1 *new_bw2 = new TF1("new_m2", "[3] / ((x * x - [4] * [4]) * (x * x - [4] * [4]) + [4] * [4] * [5] * [5])", 10, 13);
TF1 *new_bw3 = new TF1("new_m3", "[6] / ((x * x - [7] * [7]) * (x * x - [7] * [7]) + [7] * [7] * [8] * [8])", 11, 14);
TF1 *new_bw4 = new TF1("new_m4", "[9] / ((x * x - [10] * [10]) * (x * x - [10] * [10]) + [10] * [10] * [11] * [11])", 11, 15);
TF1 *new_bw5 = new TF1("new_m5", "[12] / ((x * x - [13] * [13]) * (x * x - [13] * [13]) + [13] * [13] * [14] * [14])", 11, 15);
TF1 *new_bw6 = new TF1("new_m6", "[15] / ((x * x - [16] * [16]) * (x * x - [16] * [16]) + [16] * [16] * [17] * [17])", 11, 15);

new_bw1->SetParameters(par_total);
new_bw2->SetParameters(par_total);
new_bw3->SetParameters(par_total);
new_bw4->SetParameters(par_total);
new_bw5->SetParameters(par_total);
new_bw6->SetParameters(par_total);

// Configurar nombres de parámetros
new_bw1->SetParName(0, "Amp1");
new_bw1->SetParName(1, "Mean1");
new_bw1->SetParName(2, "Width1");

new_bw2->SetParName(3, "Amp2");
new_bw2->SetParName(4, "Mean2");
new_bw2->SetParName(5, "Width2");

new_bw3->SetParName(6, "Amp3");
new_bw3->SetParName(7, "Mean3");
new_bw3->SetParName(8, "Width3");

new_bw4->SetParName(9, "Amp4");
new_bw4->SetParName(10, "Mean4");
new_bw4->SetParName(11, "Width4");

new_bw5->SetParName(12, "Amp5");
new_bw5->SetParName(13, "Mean5");
new_bw5->SetParName(14, "Width5");

new_bw6->SetParName(15, "Amp6");
new_bw6->SetParName(16, "Mean6");
new_bw6->SetParName(17, "Width6");

// Dibujar las nuevas funciones en el mismo Canvas
new_bw1->SetLineColor(kBlack);
new_bw2->SetLineColor(kBlack);
new_bw3->SetLineColor(kBlack);
new_bw4->SetLineColor(kBlack);
new_bw5->SetLineColor(kBlack);
new_bw6->SetLineColor(kBlack);

new_bw1->Draw("SAME");
new_bw2->Draw("SAME");
new_bw3->Draw("SAME");
new_bw4->Draw("SAME");
new_bw5->Draw("SAME");
new_bw6->Draw("SAME");

// Mostrar el Canvas
gPad->Update();



}

