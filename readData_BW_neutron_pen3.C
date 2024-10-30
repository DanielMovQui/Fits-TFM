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
#include <gsl/gsl_sf_bessel.h> // Biblioteca de funciones de Bessel


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
    std::cout << "Masa reducida calculada: " << masa_reducida << std::endl;
    return masa_reducida;
}

// Función para calcular el número de onda k dado E
double k(double E) {
    double E_J = E * MeV_to_J;
    double mu = calcularMasaReducida();
    if (mu <= 0) {
        std::cerr << "Error: Masa reducida no válida." << std::endl;
        return NAN; // Retornar NaN si mu no es válido
    }
    double k_value = sqrt(2 * mu * E_J) / hbar;
    std::cout << "k(E) calculado: " << k_value << std::endl;
    return k_value;
}

// Función para calcular el número de onda k0 dado E0
double k0(double E0) {
    if (E0<0){
        std::cerr << "Error: El argumento p debe ser mayor que cero para l > 0." << std::endl;
    return NAN; // Retornar NaN en caso de error
    }
    double E0_J = E0 * MeV_to_J;
    double mu = calcularMasaReducida();
    double k0_value = sqrt(2 * mu * E0_J) / hbar;
    std::cout << "k0(E0) calculado: " << k0_value << std::endl;
    return k0_value;
}

// Función para calcular la función de Bessel J_l(p) usando GSL
double besselJ_sferica(int l, double p) {
   
    return gsl_sf_bessel_jl(l, p);
}

// Función para calcular la función de Bessel Y_l(p) usando GSL
double besselY_sferica(int l, double p) {
    
    return gsl_sf_bessel_yl(l, p);
}

// Función para calcular la penetrabilidad s_l
double calcularPenetrabilidad(int l, double k_value, double r) {
    double p = k_value * r;
    double J = besselJ_sferica(l, p);
    double Y = besselY_sferica(l, p);
    double penetrabilidad = J * J + Y * Y;
    std::cout << "Penetrabilidad calculada: " << penetrabilidad << std::endl;
    return penetrabilidad;
}

// Función para calcular la anchura Γ(E)
double calcularAnchura(double E, double E0, int l, double r, double gamma_0) {
    double k_value = k(E);
    double k0_value = k0(E0);
    double P_l = k0_value * calcularPenetrabilidad(l, k_value, r);
    double P_l0 = k_value * calcularPenetrabilidad(l, k0_value, r);
    // Comprobación de valores válidos
    if (std::isnan(P_l) || std::isnan(P_l0) || P_l == 0) {
        std::cerr << "Error en el cálculo de la penetrabilidad." << std::endl;
        return NAN;
    }
    std::cout << "gamma_0: " << gamma_0 << std::endl;
    std::cout << "P_l: " << P_l << std::endl;
    std::cout << "P_l0: " << P_l0 << std::endl;
    return gamma_0 * (P_l0 / P_l);
}

double BWModificada(double *x, double *par) {
    double E = x[0];                 // Energía actual en el ajuste
    double Amp = par[0];             // Amplitud
    double E0 = par[1];              // Energía central
    double Gamma0 = par[2];          // Anchura en E0
    int l = static_cast<int>(par[3]);// Momento angular l
    double R = par[4];               // Radio nuclear

    // Cálculo de la anchura dependiente de la energía actual E
    double GammaE = calcularAnchura(E, E0, l, R, Gamma0);

    // Comprueba que GammaE no sea NaN o un valor negativo que pueda causar problemas
    if (GammaE < 0) {
        GammaE = 0;  // Evitar resultados negativos
    }

    double denominator = (E * E - E0 * E0) * (E * E - E0 * E0) + E0 * E0 * GammaE * GammaE;
if (denominator == 0) {
    return 0; // o un valor apropiado
}
return Amp / denominator;
}

// Definición de la función total
double TotalFunction(double *x, double *par) {
    double E = x[0];
    double total = 0.0;

    // Número de funciones BW
    const int numBW = 6;

    for (int i = 0; i < numBW; i++) {
        // Cada BW tiene 5 parámetros, así que se multiplica por 5
        double Amp = par[i * 5];         // Amplitud
        double E0 = par[i * 5 + 1];      // Energía central
        double Gamma = par[i * 5 + 2];   // Anchura
        int l = static_cast<int>(par[i * 5 + 3]); // Momento angular l
        double R = par[i * 5 + 4];        // Radio

        // Sumar la contribución de la función BW
        total += BWModificada(&E, par + i * 5);
    }
    return total;
}


void readData_BW_neutron_pen3()
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

// Configurar las funciones Breit-Wigner modificadas
    TF1 *bw1 = new TF1("bw1", BWModificada, 11.2, 12.0, 5);
    TF1 *bw2 = new TF1("bw2", BWModificada, 11.5, 12.5, 5);
    TF1 *bw3 = new TF1("bw3", BWModificada, 12.0, 13.0, 5);
    TF1 *bw4 = new TF1("bw4", BWModificada, 12.5, 14.0, 5);
    TF1 *bw5 = new TF1("bw5", BWModificada, 12.8, 14.4, 5);
    TF1 *bw6 = new TF1("bw6", BWModificada, 13.7, 14.5, 5);

    // Parámetros iniciales para cada BW: Amplitud, Energía Central, Gamma0, Momento angular l, Radio r
    bw1->SetParameters(500, 11.6, 0.180, 1, r);
    bw2->SetParameters(1500, 11.9, 0.294, 2, r);
    bw3->SetParameters(2800, 12.5, 0.50, 1, r);
    bw4->SetParameters(2700, 13.3, 0.40, 0, r);
    bw5->SetParameters(2000, 13.6, 0.3, 0, r);
    bw6->SetParameters(1000, 14.0, 0.3, 5, r);


  // Crear la función total
TF1 *total = new TF1("mstotal", TotalFunction, 11, 15, 6 * 5);
  
// Ajustar cada función a los datos teniendo en cuenta la anterior
exTotalH->Fit(bw1, "R");
exTotalH->Fit(bw2, "R+");
exTotalH->Fit(bw3, "R+");
exTotalH->Fit(bw4, "R+");
exTotalH->Fit(bw5, "R+");
exTotalH->Fit(bw6, "R+");

// Obtener los parámetros del fit
Double_t par[30];

bw1->GetParameters(&par[0]); // 0-4: bw1
bw2->GetParameters(&par[5]); // 5-9: bw2
bw3->GetParameters(&par[10]); // 10-14: bw3
bw4->GetParameters(&par[15]); // 15-19: bw4
bw5->GetParameters(&par[20]); // 20-24: bw5
bw6->GetParameters(&par[25]); // 25-29: bw6


total->SetParameters(par);

// Configurar nombres de parámetros
for (int i = 0; i < 6; i++) {
    total->SetParName(i * 5 + 0, Form("Amp%d", i + 1));
    total->SetParName(i * 5 + 1, Form("Mean%d", i + 1));
    total->SetParName(i * 5 + 2, Form("Width%d", i + 1));
    total->SetParName(i * 5 + 3, Form("l%d", i + 1));
    total->SetParName(i * 5 + 4, Form("R%d", i + 1));
}


// Establecer límites de parámetros 
total->SetParLimits(1, 11.50, 11.7); // Mean1
total->SetParLimits(2, 0.1, 0.4);    // Width1
total->SetParLimits(6, 11.85, 12.05); // Mean2
total->SetParLimits(7, 0.1, 0.4);    // Width2
total->SetParLimits(11, 12.4, 12.6); // Mean3
total->SetParLimits(12, 0.180, 1);   // Width3
total->SetParLimits(16, 13.25, 13.35); // Mean4
total->SetParLimits(17, 0.1, 0.4);   // Width4
total->SetParLimits(21, 13.55, 13.65); // Mean5
total->SetParLimits(22, 0.1, 0.4);   // Width5
total->SetParLimits(26, 14.05, 14.15); // Mean6
total->SetParLimits(27, 0.1, 0.4);   // Width6

// Ajustar solo la función total y desactivar la visualización de la línea de ajuste resultante
exTotalH->Fit(total, "R+");

// Configurar opciones de visualización para la línea de ajuste resultante
total->SetLineColor(kRed);
total->SetLineWidth(2);
total->SetNpx(5000);
exTotalH->SetTitle("Energy distribution ^{10}B"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");
// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[30];

total->GetParameters(par_total);

// Al definir las nuevas funciones, asegurarse de que utilizan BWModificada
TF1 *new_bw1 = new TF1("new_bw1", BWModificada, 10, 13);
TF1 *new_bw2 = new TF1("new_bw2", BWModificada, 10, 13);
TF1 *new_bw3 = new TF1("new_bw3", BWModificada, 11, 14);
TF1 *new_bw4 = new TF1("new_bw4", BWModificada, 11, 15);
TF1 *new_bw5 = new TF1("new_bw5", BWModificada, 11, 15);
TF1 *new_bw6 = new TF1("new_bw6", BWModificada, 11, 15);

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

// Configurar la leyenda
TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
legend->AddEntry(new_bw1, "BW individual fits"); // Agregar entrada para la primera función BW ajustada
legend->AddEntry(total, "Total Fit"); // Agregar entrada para el ajuste total
legend->SetBorderSize(0); // Sin borde
legend->Draw(); // Dibujar la leyenda

// Mostrar el Canvas
gPad->Update();

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
}


