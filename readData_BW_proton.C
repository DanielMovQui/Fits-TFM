
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


// Constantes de masa y conversiones
const double u_to_kg = 1.66053906660e-27; // Unidad de masa atómica a kg
const double MeV_to_J = 1.60218e-13; // Conversión de MeV a J
const double hbar = 1.0545718e-34; // Planck constante reducida en J·s
const double mass_B10_u = 10.012937; // Masa de boro-10 en unidades de masa atómica
const double mass_neutron_u = 1.008665; // Masa de neutrón en unidades de masa atómica

double binContents[] = {
    0.286418,
    0.340845,
    0.415375,
    0.462497,
    0.50954,
    0.521469,
    0.59722,
    0.556928,
    0.577389,
    0.64948,
    0.619734,
    0.663233,
    0.717629,
    0.783121,
    0.77146,
    0.729515,
    0.828862,
    0.807528,
    0.869546,
    0.79778,
    0.77868,
    0.75586,
    0.801451,
    0.863992,
    0.84627,
    0.842383,
    0.84106,
    0.856531,
    0.913934,
    0.80581,
    0.834675,
    0.856308,
    0.806832,
    0.848633,
    0.842008
};




double BW(double *x, double *par, int bw_index) {

      // Definir los rangos de energía para cada Breit-Wigner (E_min y E_max)
        static const double E_min_values[8] = {10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9};  // Energías mínimas
        static const double E_max_values[8] = {13, 13, 13, 14.4, 14.4, 14.4, 14.4, 14.4};  // Energías máximas
      
      // Determinar el bin correspondiente a la energía E
    

    // Calcular la energía central del bin

    // Definir los rangos de energía para cada Breit-Wigner (E_min y E_max)
    int num_bins = 10000; 
    double bin_width = (E_max_values[bw_index] - E_min_values[bw_index]) / num_bins;  // Ancho de cada bin

    double E = x[0];                 // Energía actual en el ajuste (energía de los datos)
    double Amp = par[0];             // Amplitud
    double E0 = par[1];              // Energía central
    double Gamma0 = par[2];          // Anchura en E0
    double sigma = par[3];           // Resolución (desviación estándar de la Gaussiana)             // Radio nuclear 
    int bin_index = (E - E_min_values[bw_index]) / bin_width;
    double E_bin = E_min_values[bw_index] + (bin_index + 0.5) * bin_width;

    // Calculando el denominador de la función Breit-Wigner usando la energía central del bin
    double denominator = (E_bin * E_bin - E0 * E0) * (E_bin * E_bin - E0 * E0) + Gamma0 * Gamma0 * E0 * E0;

    double gaussian = exp(-0.5 * pow((E - E_bin) / sigma, 2)) / (sigma * sqrt(2 * M_PI));  // Distribución Gaussiana

    return Amp * gaussian / denominator  ;
}

double TotalFunction(double *x, double *par) {
    double E = x[0];
    double total = 0.0;

    // Número de funciones BW
    const int numBW = 8;

    for (int i = 0; i < numBW; i++) {
        total += BW(&E, par + i * 4, i);
    }

    int bin_index = static_cast<int>((E - 10.9) / (14.4 - 10.9) * 35);
    if (bin_index >= 0 && bin_index < 35) {
        total += binContents[bin_index];
        
    }


    return total;
}


void readData_BW_proton()
{
    TFile *f = new TFile("h082_10BDP_trace_run013_015-019_025-041.root");
    TTree *tree = (TTree*)f->Get("tree");

    tree->SetBranchAddress("Ex", &Ex);
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

    // Definición de cortes de protones

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

TH1F* exTotalH = new TH1F("exTotalH","exTotalH",50,9,14.4);      


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

    if (!cutProtonRecoil1->IsInside(rdt[0],rdt[1])  && !cutProtonRecoil2->IsInside(rdt[2],rdt[3]) && !cutProtonRecoil3->IsInside(rdt[4],rdt[5]) && !cutProtonRecoil4->IsInside(rdt[6],rdt[7]))
      continue; 

    // if (!cutBoronRecoil1->IsInside(rdt[0],rdt[1])  && !cutBoronRecoil2->IsInside(rdt[2],rdt[3]) && !cutBoronRecoil3->IsInside(rdt[4],rdt[5]) && !cutBoronRecoil4->IsInside(rdt[6],rdt[7]))
      //continue; 
    
      exTotalH->Fill(Ex);



   


for(auto i=0;i<24;++i){
        eH[i]->Fill(e[i]);
        if(detID==i) exH[i]->Fill(Ex);
    }

     for(auto i=0;i<4;++i)
        rdtH[i]->Fill(rdt[i*2],rdt[i*2+1]);

  }//events

  // Eficiencia de que 10Be y proton se detecten en el mismo sector del detector.
double binEfficiency10kcoincsector[] = {0.0593, 0.0282, 0.0156, 0.0121, 0.0085, 0.0062, 0.0050, 0.0048, 0.0049,
0.0043, 0.0038, 0.0036, 0.0032, 0.0028, 0.0028, 0.0031, 0.0035, 0.0029, 0.0027, 0.0029, 0.0030, 0.003, 0.0034, 
0.0037, 0.0036, 0.0033, 0.0042, 0.0036, 0.0031, 0.0031, 0.0036, 0.0037};

// Eficiencia de detectar el primer proton y 10Be, menos las coincidencias de 10Be y segundo proton en sector del recoil.
double binEfficiency10k[] = {1,1,1,1, 0.1297, 0.1546, 0.1623, 0.1641, 0.1622, 0.1539, 0.1495, 0.1461, 
0.1452, 0.1410, 0.1390, 0.1317, 0.1256, 0.1258, 0.1216, 0.1171, 0.1167, 0.1117, 0.1106, 0.1074, 0.1062, 0.1009, 0.1007, 
0.0951, 0.0931, 0.0900, 0.0861, 0.0821, 0.0792, 0.0765, 0.0763, 0.0801};

/*
  double binEfficiency10k[] = { 0.3776, 0.3828, 0.3909, 0.3890, 0.3906, 0.3883, 0.3943, 0.3899, 0.3877, 0.3833, 0.3774, 0.3727, 0.3706, 0.3664, 0.3659,
    0.3584, 0.3516, 0.3439, 0.3432, 0.3454, 0.3499, 0.3449, 0.3429, 0.3352, 0.3328, 0.3302, 
    0.3253, 0.3212, 0.3218, 0.3125, 0.3059, 0.3031, 0.3006, 0.2979, 0.3053, 0.3021
                                
    };

*/

/*
  // Después de llenar el histograma exTotalH
for (int i = 1; i <= exTotalH->GetNbinsX(); i++) {
  double binCenter = exTotalH->GetBinCenter(i);
  if (binCenter >= 10.9 && binCenter <= 14.4) {
    int efficiencyIndex = static_cast<int>((binCenter - 10.9) / (14.4 - 10.9) * 35);
    if (efficiencyIndex >= 0 && efficiencyIndex < 35) {
      double efficiency = binEfficiency10k[efficiencyIndex];
      if (efficiency > 0) {
        double content = exTotalH->GetBinContent(i);
        double error = exTotalH->GetBinError(i);
        exTotalH->SetBinContent(i, content / efficiency);
        exTotalH->SetBinError(i, error / efficiency);
      }
    }
  }
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
      // Dibujar el corte correspondiente en cada RDT
    if (i == 0) cutProtonRecoil1->Draw("same");
    if (i == 1) cutProtonRecoil3->Draw("same");
    if (i == 2) cutProtonRecoil2->Draw("same");
    if (i == 3) cutProtonRecoil4->Draw("same");
    } 

  TCanvas *c4 = new TCanvas();
  coinTimeH->Draw();

*/
/*
gROOT->SetBatch(kFALSE);

TH1F *h = new TH1F("", "", 35, 10.9, 14.4);

for (int i = 1; i <= 35; i++) {
    h->SetBinContent(i, binContents[i-1]);
}

TF1 *bw1 = new TF1("bw1", [=](double *x, double *par) { 
    return BW(x, par, 0); 
}, 10.9, 13, 4);

TF1 *bw2 = new TF1("bw2", [=](double *x, double *par) { 
    return BW(x, par, 1); 
}, 10.9, 13, 4);

TF1 *bw3 = new TF1("bw3", [=](double *x, double *par) { 
    return BW(x, par, 2); 
}, 10.9, 13, 4);

TF1 *bw4 = new TF1("bw4", [=](double *x, double *par) { 
    return BW(x, par, 3); 
}, 10.9, 14.4, 4);

TF1 *bw5 = new TF1("bw5", [=](double *x, double *par) { 
    return BW(x, par, 4); 
}, 10.9, 14.4, 4);

TF1 *bw6 = new TF1("bw6", [=](double *x, double *par) { 
    return BW(x, par, 5); 
}, 10.9, 14.4, 4);

TF1 *bw7 = new TF1("bw7", [=](double *x, double *par) { 
    return BW(x, par, 6); 
}, 10.9, 14.4, 4);

TF1 *bw8 = new TF1("bw8", [=](double *x, double *par) { 
    return BW(x, par, 7); 
}, 10.9, 14.4, 4);


bw1->SetParameters(300, 11.2, 0.2);  // Parámetros iniciales: amplitud, media, anchura (ancho a media altura)
bw2->SetParameters(300, 11.5, 0.2); 
bw3->SetParameters(300, 11.9, 0.2); 
bw4->SetParameters(300, 12.1, 0.2); 
bw5->SetParameters(300, 12.5, 0.2); 
bw6->SetParameters(300, 13.2, 0.2); 
bw7->SetParameters(300, 13.6, 0.2); 
bw8->SetParameters(300, 14, 0.2); 

// Definir total dada por la suma de los 16 picos
TF1 *total = new TF1("mstotal", TotalFunction, 10.5, 14.4, 32);

// Ajustar cada función a los datos teniendo en cuenta la anterior
exTotalH->Fit(bw1, "R");
exTotalH->Fit(bw2, "R+");
exTotalH->Fit(bw3, "R+");
exTotalH->Fit(bw4, "R+");
exTotalH->Fit(bw5, "R+");
exTotalH->Fit(bw6, "R+");
exTotalH->Fit(bw7, "R+");
exTotalH->Fit(bw8, "R+");

Double_t par[32]; // 

bw1->GetParameters(&par[0]);
bw2->GetParameters(&par[4]);
bw3->GetParameters(&par[8]);
bw4->GetParameters(&par[12]);
bw5->GetParameters(&par[16]);
bw6->GetParameters(&par[20]);
bw7->GetParameters(&par[24]);
bw8->GetParameters(&par[28]);

// Configurar nombres de parámetros
for (int i = 0; i < 8; i++) {
    total->SetParName(i * 4 + 0, Form("Amp%d", i + 1));
    total->SetParName(i * 4 + 1, Form("Mean%d", i + 1));
    total->SetParName(i * 4 + 2, Form("Width%d", i + 1));
    total->SetParName(i * 4 + 3, Form("Sigma%d", i + 1));
}

// Establecer límites de parámetros (si es necesario)
total->SetParLimits(1, 11.25, 11.38);
total->SetParLimits(2, 0.165, 0.177);
total->FixParameter(3, 0.0657);
total->SetParLimits(5, 11.45, 11.55);
total->SetParLimits(6, 0.150, 0.17);
total->FixParameter(7, 0.0657);
total->SetParLimits(9, 11.8, 11.9);
total->SetParLimits(10, 0.227, 0.229);
total->FixParameter(11, 0.0657);
total->SetParLimits(13, 12, 12.12);
total->SetParLimits(14, 0.148, 0.296);
total->FixParameter(15, 0.0657);
total->SetParLimits(17, 12.4, 12.6);
total->SetParLimits(18, 0.403, 0.629);
total->FixParameter(19, 0.0657);
total->SetParLimits(21, 13.25, 13.35);
total->SetParLimits(22, 0.340, 0.420);
total->FixParameter(23, 0.0657);
total->SetParLimits(25, 13.55, 13.65);
total->SetParLimits(26, 0.222, 0.242);
total->FixParameter(27, 0.0657);
total->SetParLimits(29, 13.97, 14.03);
total->SetParLimits(30, 0.224, 0.288);
total->FixParameter(31, 0.0657);

// Ajustar solo la función total y desactivar la visualización de la línea de ajuste resultante
exTotalH->Fit(total, "R+");

// Configurar opciones de visualización para la línea de ajuste resultante
total->SetLineColor(kRed);
total->SetLineWidth(2);
total->SetNpx(10000);

exTotalH->SetTitle("Energy distribution ^{10}Be + p"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");
exTotalH->SetStats(0);

// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
h->SetFillColor(kGreen);
h->Draw("HIST SAME");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[32];
total->GetParameters(par_total);

TF1 *new_bw1 = new TF1("new_bw1", [=](double *x, double *par_total) { 
    return BW(x, par_total, 0); 
}, 10.9, 13, 4);

TF1 *new_bw2 = new TF1("new_bw2", [=](double *x, double *par_total) { 
    return BW(x, par_total, 1); 
}, 10.9, 13, 4);

TF1 *new_bw3 = new TF1("new_bw3", [=](double *x, double *par_total) { 
    return BW(x, par_total, 2); 
}, 10.9, 13, 4);

TF1 *new_bw4 = new TF1("new_bw4", [=](double *x, double *par_total) { 
    return BW(x, par_total, 3); 
}, 10.9, 14.4, 4);

TF1 *new_bw5 = new TF1("new_bw5", [=](double *x, double *par_total) { 
    return BW(x, par_total, 4); 
}, 10.9, 14.4, 4);

TF1 *new_bw6 = new TF1("new_bw6", [=](double *x, double *par_total) { 
    return BW(x, par_total, 5); 
}, 10.9, 14.4, 4);

TF1 *new_bw7 = new TF1("new_bw7", [=](double *x, double *par_total) { 
    return BW(x, par_total, 6); 
}, 10.9, 14.4, 4);

TF1 *new_bw8 = new TF1("new_bw8", [=](double *x, double *par_total) { 
    return BW(x, par_total, 7); 
}, 10.9, 14.4, 4);

new_bw1->SetNpx(10000);
new_bw2->SetNpx(10000);
new_bw3->SetNpx(10000);
new_bw4->SetNpx(10000);
new_bw5->SetNpx(10000);
new_bw6->SetNpx(10000);
new_bw7->SetNpx(10000);
new_bw8->SetNpx(10000);


new_bw1->SetParameters(par_total);
new_bw2->SetParameters(par_total+4);
new_bw3->SetParameters(par_total+8);
new_bw4->SetParameters(par_total+12);
new_bw5->SetParameters(par_total+16);
new_bw6->SetParameters(par_total+20);
new_bw7->SetParameters(par_total+24);
new_bw8->SetParameters(par_total+28);

new_bw1->SetLineColor(kBlack);
new_bw2->SetLineColor(kBlack);
new_bw3->SetLineColor(kBlack);
new_bw4->SetLineColor(kBlack);
new_bw5->SetLineColor(kBlack);
new_bw6->SetLineColor(kBlack);
new_bw7->SetLineColor(kBlack);
new_bw8->SetLineColor(kBlack);

new_bw1->Draw("SAME");
new_bw2->Draw("SAME");
new_bw3->Draw("SAME");
new_bw4->Draw("SAME");
new_bw5->Draw("SAME");
new_bw6->Draw("SAME");
new_bw7->Draw("SAME");
new_bw8->Draw("SAME");
*/
// Configurar la leyenda
exTotalH->SetTitle("Energy distribution ^{10}Be + p"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");
exTotalH->SetStats(0);

// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
//legend->AddEntry(new_bw1, "BW individual fits"); // Agregar entrada para la primera función BW ajustada
//legend->AddEntry(total, "Total Fit"); 
//legend->AddEntry(h, "Phase Space"); 
legend->SetBorderSize(0); // Sin borde
legend->Draw(); // Dibujar la leyenda

// Mostrar el Canvas
gPad->Update();
/*
// Calcular la integral de cada función Breit-Wigner ajustada
Double_t integral_bw1 = new_bw1->Integral(10.9, 14.4) * 10;
Double_t integral_bw2 = new_bw2->Integral(10.9, 14.4) * 10;
Double_t integral_bw3 = new_bw3->Integral(10.9, 14.4) * 10;
Double_t integral_bw4 = new_bw4->Integral(10.9, 14.4) * 10;
Double_t integral_bw5 = new_bw5->Integral(10.9, 14.4) * 10;
Double_t integral_bw6 = new_bw6->Integral(10.9, 14.4) * 10;
Double_t integral_bw7 = new_bw7->Integral(10.9, 14.4) * 10;
Double_t integral_bw8 = new_bw8->Integral(10.9, 14.4) * 10;


// Imprimir los resultados
cout << "Integral of new_bw1: " << integral_bw1 << endl;
cout << "Integral of new_bw2: " << integral_bw2 << endl;
cout << "Integral of new_bw3: " << integral_bw3 << endl;
cout << "Integral of new_bw4: " << integral_bw4 << endl;
cout << "Integral of new_bw5: " << integral_bw5 << endl;
cout << "Integral of new_bw6: " << integral_bw6 << endl;
cout << "Integral of new_bw7: " << integral_bw7 << endl;
cout << "Integral of new_bw8: " << integral_bw8 << endl;
*/


}

