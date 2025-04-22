#include <TMath.h> 
#include <iostream>
#include <cmath> 
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TF1.h>
#include <TROOT.h>
#include <TPad.h>
#include <math.h>
#include <TLegend.h>

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
   cutBoronRecoil1->Draw(""); 

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
   cutBoronRecoil2->Draw("");

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
   cutBoronRecoil3->Draw("");

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
   cutBoronRecoil4->Draw("");


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

TH1F* exTotalH = new TH1F("exTotalH","exTotalH",100,11,15);      


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

      if (x < -0.95 || x > 0.95 || thetaCM < 10|| e[detID] < 1)
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
/*
  double binEfficiency10k[60] = {
    0.3828, 0.3859, 0.3891, 0.3910, 0.3934, 0.3937, 0.3945, 0.3905, 0.3888, 0.3890, 0.3855, 0.3835, 0.3804, 0.3803, 0.3770, 0.3745, 0.3695, 0.3637, 0.3645, 0.3659, 0.3652, 0.3591, 0.3518, 0.3489, 0.3463, 0.3433, 0.3427, 0.3434, 0.3481, 0.3483, 0.3486, 0.3416, 0.3417, 0.3371, 0.3314, 0.3336, 0.3315, 0.3287, 0.3246, 0.3205, 
    0.3210, 0.3158, 0.3154, 0.3097, 0.3091, 0.3002, 0.2977, 0.2960, 0.2984, 0.3032, 0.3053, 0.3014, 0.2991, 0.2934, 0.2894, 0.2862, 0.2864, 0.2820, 0.2778, 0.2773
    };

for (int i = 1; i <= 60 ; i++) {  
    int efficiencyIndex = i - 1 ;
    double efficiency = binEfficiency10k[efficiencyIndex];
      double content = exTotalH->GetBinContent(i);
      double error = exTotalH->GetBinError(i);
      double content_real = content / efficiency;
      exTotalH->SetBinContent(i, content_real);    
    std::cout << "Content of the bin:"<< content<<"\n";
      std::cout << "Content of the bin efficiency:"<< content_real<<"\n";
      std::cout << "Efficiency of the bin:"<< efficiency<<"\n";

}

*/

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
/*
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
total->SetNpx(5000);
*/
exTotalH->SetTitle("Energy distribution ^{10}B + n, ThetaCM = (40, 45)"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");
// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
//total->Draw("SAME");
/*
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

new_bw1->SetNpx(5000);
new_bw2->SetNpx(5000);
new_bw3->SetNpx(5000);
new_bw4->SetNpx(5000);
new_bw5->SetNpx(5000);
new_bw6->SetNpx(5000);

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
*/
// Configurar la leyenda
TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
//legend->AddEntry(new_bw1, "BW individual fits"); // Agregar entrada para la primera función BW ajustada
//legend->AddEntry(total, "Total Fit"); // Agregar entrada para el ajuste total
legend->SetBorderSize(0); // Sin borde
legend->Draw(); // Dibujar la leyenda

// Mostrar el Canvas
gPad->Update();
/*
// Calcular las integrales de cada función sobre su dominio
double integral1 = new_bw1->Integral(10, 15);
double integral2 = new_bw2->Integral(10, 15);
double integral3 = new_bw3->Integral(10, 15);
double integral4 = new_bw4->Integral(10, 15);
double integral5 = new_bw5->Integral(10, 15);
double integral6 = new_bw6->Integral(10, 15);

// Imprimir los resultados
std::cout << "Integral de new_bw1 en [10, 13]: " << integral1 << std::endl;
std::cout << "Integral de new_bw2 en [10, 13]: " << integral2 << std::endl;
std::cout << "Integral de new_bw3 en [11, 14]: " << integral3 << std::endl;
std::cout << "Integral de new_bw4 en [11, 15]: " << integral4 << std::endl;
std::cout << "Integral de new_bw5 en [11, 15]: " << integral5 << std::endl;
std::cout << "Integral de new_bw6 en [11, 15]: " << integral6 << std::endl;

*/
}

