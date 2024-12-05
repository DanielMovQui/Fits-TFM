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
const double r = 85 * 1.0e-12;

// Función para calcular la masa reducida
double calcularMasaReducida() {
    double mass_B10_kg = mass_B10_u * u_to_kg;
    double mass_neutron_kg = mass_neutron_u * u_to_kg;
    double masa_reducida = (mass_B10_kg * mass_neutron_kg) / (mass_B10_kg + mass_neutron_kg);
    return masa_reducida;
}

double BWModificada(double *x, double *par, int l, int bw_index) {
    double mu = calcularMasaReducida();

    // Definir los rangos de energía para cada Breit-Wigner (E_min y E_max)
    static const double E_min_values[6] = {11, 11, 11, 12, 12, 12};  // Energías mínimas
    static const double E_max_values[6] = {14.5, 14.5, 14.5, 14.5, 14.5, 14.5};  // Energías máximas

    double E = x[0];                 // Energía actual en el ajuste (energía de los datos)
    double Amp = par[0];             // Amplitud
    double E0 = par[1];              // Energía central
    double Gamma0 = par[2];          // Anchura en E0
    double r = par[3];               // Radio nuclear 

    // Asegurarse de que bw_index está dentro del rango válido
    if (bw_index < 0 || bw_index >= 6) {
        return 0.0; // Retornar 0 si el índice es inválido
    }

    // Número total de bines en el histograma
    int num_bins = 60; 
    double bin_width = (E_max_values[bw_index] - E_min_values[bw_index]) / num_bins;  // Ancho de cada bin

    // Determinar el bin correspondiente a la energía E
    int bin_index = (E - E_min_values[bw_index]) / bin_width;
    if (bin_index < 0 || bin_index >= num_bins) {
        // Si la energía está fuera del rango del histograma
        return 0.0;
    }

    // Calcular la energía central del bin
    double E_bin = E_min_values[bw_index] + (bin_index + 0.5) * bin_width;


    // Calculando las penetrabilidades para diferentes valores de l
    double Gamma_pen_l0 = Gamma0 * sqrt(E_bin / E0);
    double Gamma_pen_l1 = Gamma0 * pow(sqrt(E_bin / E0), 3.0 / 2) * (2 * E0 / (E_bin + E0)) * 
        ((1 + (2 * MeV_to_J * mu * E0 * pow(r, 2) / pow(hbar, 2))) / 
         (1 + (2 * MeV_to_J * mu * E_bin * pow(r, 2) / pow(hbar, 2))));
    double Gamma_pen_l2 = Gamma0 * pow(sqrt(E_bin / E0), 5.0 / 2) * (2 * E0 / (E_bin + E0)) * 
        ((9 + (6 * mu * MeV_to_J * E0 * pow(r, 2) / pow(hbar, 2)) + 
         (2 * mu * MeV_to_J * MeV_to_J * E0 * pow(r, 2) / pow(hbar, 2)) * (2)) / 
         (9 + (6 * mu * MeV_to_J * E0 * pow(r, 2) / pow(hbar, 2)) + 
          (2 * mu * MeV_to_J * MeV_to_J * E_bin * pow(r, 2) / pow(hbar, 2)) * (2)));

    // Selección de la Gamma efectiva según l
    double Gamma_eff;
    if (l == 0) {
        Gamma_eff = Gamma_pen_l0;
    } else if (l == 1) {
        Gamma_eff = Gamma_pen_l1;
    } else if (l == 2) {
        Gamma_eff = Gamma_pen_l2;
    } else {
        Gamma_eff = Gamma0;  // Para otros valores de l, usar Gamma0 por defecto
    }

    // Calculando el denominador de la función Breit-Wigner usando la energía central del bin
    double denominator = (E_bin * E_bin - E0 * E0) * (E_bin * E_bin - E0 * E0) + Gamma_eff * Gamma_eff * E0 * E0;

    return Amp / denominator;
}


// Definición de la función total
double TotalFunction(double *x, double *par) {
    double E = x[0];
    double total = 0.0;

    // Número de funciones BW
    const int numBW = 6;
    int l_values[6] = {0, 1, 2, 0, 0, 5}; 

    for (int i = 0; i < numBW; i++) {
        // Cada BW tiene 4 parámetros (Amplitud, Energía central, Gamma, Radio)
        double Amp = par[i * 4];         // Amplitud
        double E0 = par[i * 4 + 1];      // Energía central
        double Gamma0 = par[i * 4 + 2];   // Anchura
        double R = par[i * 4 + 3];        // Radio

        // Sumar la contribución de la función BW
        total += BWModificada(&E, par + i * 4, l_values[i],i);
    }
    return total;
}

void readData_BW_neutron_penetrability4()
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

TH1F* exTotalH = new TH1F("exTotalH","exTotalH",60,11,15);      


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



// Crear las funciones BW con los parámetros adecuados
TF1 *bw1 = new TF1("bw1", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 0); 
}, 11, 14.5, 4);  // E_min = 11.2, E_max = 12.0, número de parámetros = 4

TF1 *bw2 = new TF1("bw2", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 1); 
}, 11, 14.5, 4);  // E_min = 11.5, E_max = 12.5, número de parámetros = 4

TF1 *bw3 = new TF1("bw3", [=](double *x, double *par) { 
    return BWModificada(x, par, 2, 2); 
}, 11, 14.5, 4);  // E_min = 12.0, E_max = 13.0, número de parámetros = 4

TF1 *bw4 = new TF1("bw4", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 3); 
}, 12, 14.5, 4);  // E_min = 12.5, E_max = 14.0, número de parámetros = 4

TF1 *bw5 = new TF1("bw5", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 4); 
}, 12, 14.5, 4);  // E_min = 12.8, E_max = 14.4, número de parámetros = 4

TF1 *bw6 = new TF1("bw6", [=](double *x, double *par) { 
    return BWModificada(x, par, 5, 5); 
}, 12, 14.5, 4);  // E_min = 13.7, E_max = 14.5, número de parámetros = 4


// Configurar los parámetros iniciales para cada función BW
bw1->SetParameters(500, 11.6, 0.180, r);
bw2->SetParameters(1500, 11.9, 0.294, r);
bw3->SetParameters(2800, 12.5, 0.50, r);
bw4->SetParameters(2700, 13.3, 0.40, r);
bw5->SetParameters(2000, 13.6, 0.3, r);
bw6->SetParameters(1000, 14.0, 0.3, r);

// Crear la función total (sumatoria de las funciones Breit-Wigner)
TF1 *total = new TF1("mstotal", TotalFunction, 11, 15, 24);

// Ajustar las funciones BW a los datos experimentales
exTotalH->Fit(bw1, "R");
exTotalH->Fit(bw2, "R+");
exTotalH->Fit(bw3, "R+");
exTotalH->Fit(bw4, "R+");
exTotalH->Fit(bw5, "R+");
exTotalH->Fit(bw6, "R+");

// Obtener los parámetros ajustados
Double_t par[24];
bw1->GetParameters(&par[0]); // 0-4: bw1
bw2->GetParameters(&par[4]); // 5-9: bw2
bw3->GetParameters(&par[8]); // 10-14: bw3
bw4->GetParameters(&par[12]); // 15-19: bw4
bw5->GetParameters(&par[16]); // 20-24: bw5
bw6->GetParameters(&par[20]); // 25-29: bw6

// Configurar los parámetros para la función total
total->SetParameters(par);

// Configurar nombres de parámetros
for (int i = 0; i < 6; i++) {
    total->SetParName(i * 4 + 0, Form("Amp%d", i + 1));
    total->SetParName(i * 4 + 1, Form("Mean%d", i + 1));
    total->SetParName(i * 4 + 2, Form("Width%d", i + 1));
    total->SetParName(i * 4 + 3, Form("R%d", i + 1));
}


// Establecer límites para cada parámetro de las funciones BW
total->SetParLimits(1, 11.65, 11.75); // Mean1
total->SetParLimits(2, 0.01, 0.4);   // Width1

total->SetParLimits(5, 11.9, 12); // Mean2
total->SetParLimits(6, 0.1, 0.4);    // Width2
total->FixParameter(7, 8.5e-11);
total->SetParLimits(9, 12.35, 12.6); // Mean3
total->SetParLimits(10, 0.1, 0.8);   // Width3
total->FixParameter(11, 8.5e-11);
total->SetParLimits(13, 13.25, 13.35); // Mean4
total->SetParLimits(14, 0.1, 0.4);   // Width4

total->SetParLimits(17, 13.55, 13.65); // Mean5
total->SetParLimits(18, 0.1, 0.4);   // Width5

total->SetParLimits(21, 14.05, 14.15); // Mean6
total->SetParLimits(22, 0.1, 0.4);   // Width6

// Ajustar la función total y visualizar el ajuste
exTotalH->Fit(total, "R+");

// Configurar opciones de visualización
total->SetLineColor(kRed);
total->SetLineWidth(2);
total->SetNpx(10000);
exTotalH->SetTitle("Energy distribution ^{10}B"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");
exTotalH->Draw("HIST");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[24];
total->GetParameters(par_total);

// Crear nuevas funciones BW con los parámetros ajustados (newBW)
TF1 *new_bw1 = new TF1("new_bw1", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 0); 
}, 11, 14.5, 4);

TF1 *new_bw2 = new TF1("new_bw2", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 1); 
}, 11, 14.5, 4);

TF1 *new_bw3 = new TF1("new_bw3", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 2, 2); 
}, 11, 14.5, 4);

TF1 *new_bw4 = new TF1("new_bw4", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 3); 
}, 12, 14.5, 4);

TF1 *new_bw5 = new TF1("new_bw5", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 4); 
}, 12, 14.5, 4);

TF1 *new_bw6 = new TF1("new_bw6", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 5, 5); 
}, 12, 14.5, 4);


new_bw1->SetNpx(10000);
new_bw2->SetNpx(10000);
new_bw3->SetNpx(10000);
new_bw4->SetNpx(10000);
new_bw5->SetNpx(10000);
new_bw6->SetNpx(10000);

new_bw1->SetParameters(par_total);
new_bw2->SetParameters(par_total + 4);
new_bw3->SetParameters(par_total + 8);
new_bw4->SetParameters(par_total + 12);
new_bw5->SetParameters(par_total + 16);
new_bw6->SetParameters(par_total + 20);

new_bw1->SetLineColor(kBlack);
new_bw1->SetLineWidth(2);  // Establecer un grosor de línea mayor
new_bw1->Draw("SAME");

new_bw2->SetLineColor(kBlack);
new_bw2->SetLineWidth(2);
new_bw2->Draw("SAME");

new_bw3->SetLineColor(kBlack);
new_bw3->SetLineWidth(2);
new_bw3->Draw("SAME");

new_bw4->SetLineColor(kBlack);
new_bw4->SetLineWidth(2);
new_bw4->Draw("SAME");

new_bw5->SetLineColor(kBlack);
new_bw5->SetLineWidth(2);
new_bw5->Draw("SAME");

new_bw6->SetLineColor(kBlack);
new_bw6->SetLineWidth(2);
new_bw6->Draw("SAME");

// Configurar la leyenda
TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
legend->AddEntry(new_bw1, "BW individual fits"); // Agregar entrada para la primera función BW ajustada
legend->AddEntry(total, "Total Fit"); // Agregar entrada para el ajuste total
legend->SetBorderSize(0); // Sin borde
legend->Draw(); // Dibujar la leyenda

// Mostrar el Canvas
gPad->Update();

}