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
#include "binContents.C"

Float_t Ex = 0;
Int_t detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x[24];
Float_t thetaCM = 0;


// Constantes de masa y conversiones
const double u_to_kg = 1.66053906660e-27; // Unidad de masa atómica a kg
const double MeV_to_J = 1.60218e-13; // Conversión de MeV a J
const double hbar = 6.582119e-16 ; // Planck constante reducida en J·s
const double mass_B10_u = 10.012937; // Masa de boro-10 en unidades de masa atómica
const double mass_neutron_u = 1.008665; // Masa de neutrón en unidades de masa atómica
const double r = 4.1 * 1.0e-15;

/*
double binEfficiency10kdata51bin[] = {
        0.3927, 0.3945, 0.3916, 0.3894, 0.3889, 0.3874, 0.3846, 0.3821, 0.3815, 0.3789,
        0.3765, 0.3715, 0.3702, 0.3689, 0.3660, 0.3624, 0.3601, 0.3575, 0.3522, 0.3475, 
        0.3435, 0.3449, 0.3403, 0.3485, 0.3464, 0.3467, 0.343, 0.3449, 0.3394, 0.3372, 
        0.3331, 0.3307, 0.3275, 0.3249, 0.3248, 0.3193, 0.3192, 0.3161, 0.3157, 0.3056,
        0.3059, 0.2997, 0.3004, 0.2976, 0.3005, 0.3025, 0.2987, 0.299, 0.2946,
        0.2918, 0.2859
    };
// Habiendo hecho un ajuste polinómico a la funcion de eficiencia.
double binEfficiency10k51ajuste[] = {
    0.393167, 0.392714, 0.391696, 0.390214, 0.388359, 0.386213, 0.383848, 0.381326, 0.378703, 0.376026,
    0.373334, 0.370660, 0.368031, 0.365466, 0.362980, 0.360584, 0.358280, 0.356072, 0.353954, 0.351922,
    0.349965, 0.348072, 0.346231, 0.344425, 0.342639, 0.340857, 0.339063, 0.337239, 0.335371, 0.333445,
    0.331449, 0.329373, 0.327208, 0.324952, 0.322603, 0.320165, 0.317646, 0.315059, 0.312421, 0.309758,
    0.307099, 0.304482, 0.301952, 0.299560, 0.297366, 0.295441, 0.293862, 0.292716, 0.292102, 0.292127,
    0.292910
};
*/


/*
// Con 60 bines
double binEfficiency10k[] = {
    0.3828, 0.3859, 0.3891, 0.3910, 0.3934, 0.3937, 0.3945, 0.3905, 0.3888, 0.389,
    0.3855, 0.3835, 0.3804, 0.3803, 0.3770, 0.3745, 0.3695, 0.3637, 0.3645, 0.3659, 
    0.3652, 0.3591, 0.3518, 0.3489, 0.3463, 0.3433, 0.3427, 0.3434, 0.3481, 0.3483, 
    0.3486, 0.3416, 0.3417, 0.3371, 0.3314, 0.3336, 0.3315, 0.3287, 0.3246, 0.3205, 
    0.3210, 0.3158, 0.3154, 0.3097, 0.3091, 0.3002, 0.2977, 0.2960, 0.2984, 0.3032,
    0.3053, 0.3014, 0.2991, 0.2991, 0.2934, 0.2894, 0.2862, 0.2864, 0.2820, 0.2778,
    0.2773
};
*/
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
        static const double E_min_values[11] = {11.3, 11.3, 11.3, 11.3, 11.3, 11.3, 11.3, 11.3, 11.3, 11.3, 11.3};  // Energías mínimas
        static const double E_max_values[11] = {14.4, 14.4, 14.4, 14.4, 14.4, 14.4, 14.4, 14.4, 14.4, 14.4, 14.4};  // Energías máximas

    double E = x[0];                 // Energía actual en el ajuste (energía de los datos)
    double Amp = par[0];             // Amplitud
    double E0 = par[1];              // Energía central
    double Gamma0 = par[2];          // Anchura en E0
    double r = par[3];  
    double sigma = par[4];           // Resolución (desviación estándar de la Gaussiana)             // Radio nuclear 

    // Asegurarse de que bw_index está dentro del rango válido
    if (bw_index < 0 || bw_index >= 11) {
        return 0.0; // Retornar 0 si el índice es inválido
    }

    // Número total de bines en el histograma
    int num_bins = 10000; 
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
    double Gamma_pen_l1 = Gamma0 * pow((E_bin / E0), 3.0 / 2) * (2 * E0 / (E_bin + E0)) * 
        ((1 + (2 * (10e12/MeV_to_J) * mu * E0 * pow(r, 2) / pow(hbar, 2))) / 
         (1 + (2 * (10e12/MeV_to_J) * mu * E_bin * pow(r, 2) / pow(hbar, 2))));
    double Gamma_pen_l2 = Gamma0 * pow((E_bin / E0), 5.0 / 2) * (2 * E0 / (E_bin + E0)) * 
        ((9 + (6 * mu * (10e12/MeV_to_J) * E0 * pow(r, 2) / pow(hbar, 2)) + 
         (2 * mu * (10e12/MeV_to_J) * (10e12/MeV_to_J) * E0 * pow(r, 2) / pow(hbar, 2)) * (2)) / 
         (9 + (6 * mu * (10e12/MeV_to_J) * E0 * pow(r, 2) / pow(hbar, 2)) + 
          (2 * mu * (10e12/MeV_to_J) * (10e12/MeV_to_J) * E_bin * pow(r, 2) / pow(hbar, 2)) * (2)));

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
   double denominator = (E_bin * E_bin - E0 * E0) * (E_bin * E_bin - E0 * E0) + (Gamma_eff * Gamma_eff * E0 * E0);
    //double denominator = (E_bin * E_bin - E0 * E0) * (E_bin * E_bin - E0 * E0) + (Gamma_eff * Gamma_eff / 4);
    double numerator = 2 * sqrt(2) * E0 * Gamma_eff * sqrt(E0*E0*((E0*E0) + (Gamma_eff*Gamma_eff))) / (M_PI * sqrt((E0*E0) + (sqrt(E0*E0*((E0*E0) + (Gamma_eff*Gamma_eff))))));
    double gaussian = exp(-0.5 * pow((E - E_bin) / sigma, 2)) / (sigma * sqrt(2 * M_PI));  // Distribución Gaussiana
    //double smeared_bw = Amp * gaussian / denominator;
    double smeared_bw = Amp * numerator * gaussian / (denominator);

    return smeared_bw  ;
}


double TotalFunction(double *x, double *par) {
    double E = x[0];
    double total = 0.0;

    // Número de funciones BW
    const int numBW = 11;
    int l_values[11] = {0, 1, 0, 2, 3, 3, 0, 0, 0, 4, 0}; 

    for (int i = 0; i < numBW; i++) {
        total += BWModificada(&E, par + i * 5, l_values[i], i);
    }

    // Sumar el contenido del bin correspondiente
    int bin_index = static_cast<int>((E - 11) / (15 - 11) * 60);
    if (bin_index >= 0 && bin_index < 60) {
        total += binContents[bin_index];
        
    }

    return total;
}


void readData_BW_neutron_penetrability_and_efficiency3()
{
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

TH1F* correctedCoinTimeH = new TH1F("correctedCoinTimeH", "Corrected Coincidence Time", 400, -1000, 1000);


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
  

 Long64_t nentries = tree->GetEntries();
 std::cout<<" Number of entries : "<<nentries<<"\n";
   for (Long64_t i=0;i<nentries;i++) {
     tree->GetEntry(i);

    if(i % 1000000 == 0){
         std::cout<<" Entry number : "<<i<<"\n";
         //if(!std::isnan(e[2])) std::cout<<" Energy index 0 "<<e[2]<<"\n"; 
    }   

      if (x[detID] < -0.95 || x[detID] > 0.95 || thetaCM < 10 || thetaCM > 45 || e[detID] < 1)
        continue;

    //if (!cutProtonRecoil1->IsInside(rdt[0],rdt[1])  && !cutProtonRecoil2->IsInside(rdt[2],rdt[3]) && !cutProtonRecoil3->IsInside(rdt[4],rdt[5]) && !cutProtonRecoil4->IsInside(rdt[6],rdt[7]))
      //continue; 

     if (!cutBoronRecoil1->IsInside(rdt[0],rdt[1])  && !cutBoronRecoil2->IsInside(rdt[2],rdt[3]) && !cutBoronRecoil3->IsInside(rdt[4],rdt[5]) && !cutBoronRecoil4->IsInside(rdt[6],rdt[7]))
      continue; 
    
      // Calculate polynomial correction
      double poly_correction = pol4[detID][0]
      + pol4[detID][1] * x[detID]
      + pol4[detID][2] * x[detID] * x[detID]
      + pol4[detID][3] * x[detID] * x[detID] * x[detID]
      + pol4[detID][4] * x[detID] * x[detID] * x[detID] * x[detID];

// Apply the correction to coinTime
double correctedCoinTime = coinTime - poly_correction;

coinTimeH->Fill(coinTime);

if(correctedCoinTime<-20 || correctedCoinTime>15)
continue; 

correctedCoinTimeH->Fill(correctedCoinTime);
exTotalH->Fill(Ex);

  }//events


double binEfficiency10k[] = {
    /*0.0003, 0.0003, 0.0003, 0.0002, 0.0001, 0.0003, 0.0002,*/1,1,1,1,0.1856,0.1856,0.1856, 0.1856, 0.1865, 0.1831,
0.1805, 0.1755, 0.1697, 0.1652, 0.1622, 0.1585, 0.1575, 0.1547, 0.1513, 0.1498, 0.1475, 0.1427, 0.1423, 
0.1391, 0.1365, 0.1299, 0.1306, 0.1258, 0.1249, 0.1225, 0.1209, 0.1191, 0.1166, 0.1128, 0.1124, 0.1115,
0.1080, 0.1095, 0.1062, 0.1065, 0.1031, 0.1019, 0.0966, 0.0945, 0.0898, 0.0894, 0.0888, 0.0843, 0.0830, 
0.0841, 0.0852, 0.0863, 0.0848, 0.0833, 0.0787, 0.0776, 0.0755, 0.0757, 0.0717, 0.0693, 0.066};
/*
  // Con 60 bines
double binEfficiency10k[] = {
    0.3828, 0.3859, 0.3891, 0.3910, 0.3934, 0.3937, 0.3945, 0.3905, 0.3888, 0.389,
    0.3855, 0.3835, 0.3804, 0.3803, 0.3770, 0.3745, 0.3695, 0.3637, 0.3645, 0.3659, 
    0.3652, 0.3591, 0.3518, 0.3489, 0.3463, 0.3433, 0.3427, 0.3434, 0.3481, 0.3483, 
    0.3486, 0.3416, 0.3417, 0.3371, 0.3314, 0.3336, 0.3315, 0.3287, 0.3246, 0.3205, 
    0.3210, 0.3158, 0.3154, 0.3097, 0.3091, 0.3002, 0.2977, 0.2960, 0.2984, 0.3032,
    0.3053, 0.3014, 0.2991, 0.2991, 0.2934, 0.2894, 0.2862, 0.2864, 0.2820, 0.2778,
    0.2773
};
*/

for (int i = 1; i <= 60 ; i++) {  
    int efficiencyIndex = i - 1 ;
    double efficiency = binEfficiency10k[efficiencyIndex];
      double content = exTotalH->GetBinContent(i);
      double error = exTotalH->GetBinError(i);
      double content_real = content / efficiency;
      exTotalH->SetBinContent(i, content_real); 
}


gROOT->SetBatch(kFALSE);

TH1F *h = new TH1F("", "", 60, 11, 15);

for (int i = 1; i <= 60; i++) {
    h->SetBinContent(i, binContents[i-1]);
}

// Crear las funciones BW con los parámetros adecuados
TF1 *bw1 = new TF1("bw1", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 0); 
}, 11.3, 14, 5);  // E_min = 11.2, E_max = 12.0, número de parámetros = 4

TF1 *bw2 = new TF1("bw2", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 1); 
}, 11.3, 14, 5);  // E_min = 11.5, E_max = 12.5, número de parámetros = 4

TF1 *bw3 = new TF1("bw3", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 2); 
}, 11.3, 14, 5);  // E_min = 11.5, E_max = 12.5, número de parámetros = 4

TF1 *bw4 = new TF1("bw4", [=](double *x, double *par) { 
    return BWModificada(x, par, 2, 3); 
}, 11.3, 14, 5);  // E_min = 11.5, E_max = 12.5, número de parámetros = 4

TF1 *bw5 = new TF1("bw3", [=](double *x, double *par) { 
    return BWModificada(x, par, 3, 4); 
}, 11.3, 14, 5);  // E_min = 11.5, E_max = 12.5, número de parámetros = 4

TF1 *bw6 = new TF1("bw5", [=](double *x, double *par) { 
    return BWModificada(x, par, 3, 5); 
}, 11.5, 14.4, 5);  // E_min = 12.0, E_max = 13.0, número de parámetros = 4

TF1 *bw7 = new TF1("bw5", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 5); 
}, 11.5, 14.4, 5);  // E_min = 12.0, E_max = 13.0, número de parámetros = 4

TF1 *bw8 = new TF1("bw6", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 6); 
}, 11.5, 14.4, 5);  // E_min = 12.5, E_max = 14.0, número de parámetros = 4

TF1 *bw9 = new TF1("bw7", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 7); 
}, 11.5, 14.4, 5);  // E_min = 12.8, E_max = 14.4, número de parámetros = 4

TF1 *bw10 = new TF1("bw8", [=](double *x, double *par) { 
    return BWModificada(x, par, 4, 8); 
}, 11.5, 14.4, 5);  // E_min = 13.7, E_max = 14.5, número de parámetros = 4

TF1 *bw11 = new TF1("bw9", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 9); 
}, 11.5, 14.4, 5);  // E_min = 13.7, E_max = 14.5, número de parámetros = 4

// Configurar los parámetros iniciales para cada función BW
bw1->SetParameters(500, 11.6, 0.180, r, 0.0657);
bw2->SetParameters(1500, 11.9, 0.294, r, 0.0657);
bw3->SetParameters(2800, 12, 1, r, 0.0657);
bw4->SetParameters(2800, 12.3, 0.50, r, 0.0657);
bw5->SetParameters(2800, 12.6, 0.50, r, 0.0657);
bw6->SetParameters(2800, 12.9, 0.20, r, 0.0657);
bw7->SetParameters(2800, 13.1, 0.20, r, 0.0657);
bw8->SetParameters(2700, 13.3, 0.40, r, 0.0657);
bw9->SetParameters(2000, 13.6, 0.3, r, 0.0657);
bw10->SetParameters(1000, 13.8, 0.3, r, 0.0657);
bw11->SetParameters(1000, 14, 0.3, r, 0.0657);

// Crear la función total (sumatoria de las funciones Breit-Wigner)
TF1 *total = new TF1("mstotal", TotalFunction, 11.3, 14.4, 55);

// Ajustar las funciones BW a los datos experimentales
exTotalH->Fit(bw1, "R");
exTotalH->Fit(bw2, "R+");
exTotalH->Fit(bw3, "R+");
exTotalH->Fit(bw4, "R+");
exTotalH->Fit(bw5, "R+");
exTotalH->Fit(bw6, "R+");
exTotalH->Fit(bw7, "R+");
exTotalH->Fit(bw8, "R+");
exTotalH->Fit(bw9, "R+");
exTotalH->Fit(bw10, "R+");
exTotalH->Fit(bw11, "R+");

// Obtener los parámetros ajustados
Double_t par[55];
bw1->GetParameters(&par[0]); 
bw2->GetParameters(&par[5]); 
bw3->GetParameters(&par[10]); 
bw4->GetParameters(&par[15]); 
bw5->GetParameters(&par[20]);
bw6->GetParameters(&par[25]); 
bw7->GetParameters(&par[30]); 
bw8->GetParameters(&par[35]); 
bw9->GetParameters(&par[40]); 
bw10->GetParameters(&par[45]); 
bw11->GetParameters(&par[50]); 

// Configurar los parámetros para la función total
total->SetParameters(par);

// Configurar nombres de parámetros
for (int i = 0; i < 11; i++) {
    total->SetParName(i * 5 + 0, Form("Amp%d", i + 1));
    total->SetParName(i * 5 + 1, Form("Mean%d", i + 1));
    total->SetParName(i * 5 + 2, Form("Width%d", i + 1));
    total->SetParName(i * 5 + 3, Form("R%d", i + 1));
    total->SetParName(i * 5 + 4, Form("Sigma%d", i + 1));
}

/*
// Establecer límites para cada parámetro de las funciones BW
total->SetParLimits(0, 1, 10000); 
total->SetParLimits(1, 11.69, 11.71); 
total->SetParLimits(2, 0.2, 0.25);
total->FixParameter(3, 4.1e-15); 
total->FixParameter(4, 0.0657);

total->SetParLimits(5, 1, 10000); 
total->SetParLimits(6, 11.82, 11.85); 
total->SetParLimits(7, 0.24, 0.5);    
total->FixParameter(8, 4.1e-15);
total->FixParameter(9, 0.0657);

total->SetParLimits(10, 1, 10000); 
total->SetParLimits(11, 12.0, 12.05); 
total->FixParameter(12, 0.330);   
total->FixParameter(13, 4.1e-15);
total->FixParameter(14, 0.0657);

total->SetParLimits(15, 1, 5000); 
total->FixParameter(16, 12.475); 
total->FixParameter(17, 0.480);  
total->FixParameter(18, 4.1e-15);
total->FixParameter(19, 0.0657);

total->SetParLimits(20, 1, 10000); 
total->FixParameter(21, 12.876); 
total->FixParameter(22, 0.410); 
total->FixParameter(23, 4.1e-15);
total->FixParameter(24, 0.0657);

total->SetParLimits(25, 1, 10000); 
total->SetParLimits(26, 13.1, 13.25); 
total->SetParLimits(27, 0.190, 1); 
total->FixParameter(28, 4.1e-15);
total->FixParameter(29, 0.0657);

total->SetParLimits(30, 1, 10000); 
total->FixParameter(31, 13.368); 
total->SetParLimits(32, 0.17, 0.5); 
total->FixParameter(33, 4.1e-15);
total->FixParameter(34, 0.0657);

total->SetParLimits(35, 1, 10000); 
total->FixParameter(36, 13.589); 
total->FixParameter(37, 0.237);  
total->FixParameter(38, 4.1e-15);
total->FixParameter(39, 0.0657);

total->SetParLimits(40, 1, 10000); 
total->FixParameter(41, 13.8); 
total->FixParameter(42, 0.2);  
total->FixParameter(43, 4.1e-15);
total->FixParameter(44, 0.0657);

total->SetParLimits(45, 1, 10000); 
total->SetParLimits(46, 14, 14.05); 
total->SetParLimits(47, 0.23, 0.5);  
total->FixParameter(48, 4.1e-15);
total->FixParameter(49, 0.0657);

total->SetParLimits(50, 1, 1000); 
total->FixParameter(51, 14.198); 
total->SetParLimits(52, 0.18, 0.3); 
total->FixParameter(53, 4.1e-15);
total->FixParameter(54, 0.0657);
*/

total->SetParLimits(0, 1, 10000); 
total->SetParLimits(1, 11.67, 11.7); 
total->SetParLimits(2, 0.2, 0.25);
total->FixParameter(3, 4.1e-15); 
total->FixParameter(4, 0.0657);

total->SetParLimits(5, 1, 10000); 
total->SetParLimits(6, 11.8, 11.9); 
total->SetParLimits(7, 0.1, 0.5);    
total->FixParameter(8, 4.1e-15);
total->FixParameter(9, 0.0657);

total->SetParLimits(10, 1, 10000); 
total->SetParLimits(11, 12.0, 12.1); 
total->SetParLimits(12, 0.1, 0.5);   
total->FixParameter(13, 4.1e-15);
total->FixParameter(14, 0.0657);

total->SetParLimits(15, 1, 5000); 
total->SetParLimits(16, 12.45, 12.5); 
total->SetParLimits(17, 0.3, 0.6);  
total->FixParameter(18, 4.1e-15);
total->FixParameter(19, 0.0657);

total->SetParLimits(20, 1, 10000); 
total->SetParLimits(21, 12.85, 12.9); 
total->SetParLimits(22, 0.2, 0.5); 
total->FixParameter(23, 4.1e-15);
total->FixParameter(24, 0.0657);

total->SetParLimits(25, 1, 10000); 
total->SetParLimits(26, 13.1, 13.25); 
total->SetParLimits(27, 0.190, 1); 
total->FixParameter(28, 4.1e-15);
total->FixParameter(29, 0.0657);

total->SetParLimits(30, 1, 10000); 
total->SetParLimits(31, 13.35, 13.38); 
total->SetParLimits(32, 0.17, 0.5); 
total->FixParameter(33, 4.1e-15);
total->FixParameter(34, 0.0657);

total->SetParLimits(35, 1, 10000); 
total->SetParLimits(36, 13.55, 13.6); 
total->SetParLimits(37, 0.18, 0.3);  
total->FixParameter(38, 4.1e-15);
total->FixParameter(39, 0.0657);

total->SetParLimits(40, 1, 10000); 
total->SetParLimits(41, 13.75, 13.8); 
total->SetParLimits(42, 0.15, 0.4);  
total->FixParameter(43, 4.1e-15);
total->FixParameter(44, 0.0657);

total->SetParLimits(45, 19, 10000); 
total->SetParLimits(46, 13.85, 13.95); 
total->SetParLimits(47, 0.17, 0.5);  
total->FixParameter(48, 4.1e-15);
total->FixParameter(49, 0.0657);

total->SetParLimits(50, 20, 1000); 
total->SetParLimits(51, 14., 14.2); 
total->SetParLimits(52, 0.15, 0.3); 
total->FixParameter(53, 4.1e-15);
total->FixParameter(54, 0.0657);
// Ajustar la función total y visualizar el ajuste
exTotalH->Fit(total, "R+");

// Configurar opciones de visualización
total->SetLineColor(kRed);
total->SetLineWidth(2);
total->SetNpx(10000);
exTotalH->SetTitle("Energy distribution ^{10}B + n"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");
exTotalH->SetLineColor(kBlue);
exTotalH->Draw("HIST");
exTotalH->SetStats(0);

h->SetFillColor(kGreen);
h->Draw("HIST SAME");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[55];
total->GetParameters(par_total);

// Crear nuevas funciones BW con los parámetros ajustados (newBW)
TF1 *new_bw1 = new TF1("new_bw1", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 0); 
}, 11.3, 14, 5);

TF1 *new_bw2 = new TF1("new_bw2", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 1); 
}, 11.3, 14, 5);

TF1 *new_bw3 = new TF1("new_bw3", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 2); 
}, 11.3, 14, 5);

TF1 *new_bw4 = new TF1("new_bw4", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 2, 3); 
}, 11.3, 14, 5);

TF1 *new_bw5 = new TF1("new_bw3", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 3, 4); 
}, 11.3, 14, 5);

TF1 *new_bw6 = new TF1("new_bw5", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 3, 5); 
}, 11.5, 14.4, 5);

TF1 *new_bw7 = new TF1("new_bw5", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 6); 
}, 11.5, 14.4, 5);

TF1 *new_bw8 = new TF1("new_bw6", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 7); 
}, 11.5, 14.4, 5);

TF1 *new_bw9 = new TF1("new_bw7", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 8); 
}, 11.5, 14.4, 5);

TF1 *new_bw10 = new TF1("new_bw8", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 4, 9); 
}, 11.5, 14.4, 5);

TF1 *new_bw11 = new TF1("new_bw9", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 10); 
}, 11.5, 14.4, 5);


new_bw1->SetNpx(10000);
new_bw2->SetNpx(10000);
new_bw3->SetNpx(10000);
new_bw4->SetNpx(10000);
new_bw5->SetNpx(10000);
new_bw6->SetNpx(10000);
new_bw7->SetNpx(10000);
new_bw8->SetNpx(10000);
new_bw9->SetNpx(10000);
new_bw10->SetNpx(10000);
new_bw11->SetNpx(10000);

new_bw1->SetParameters(par_total);
new_bw2->SetParameters(par_total + 5);
new_bw3->SetParameters(par_total + 10);
new_bw4->SetParameters(par_total + 15);
new_bw5->SetParameters(par_total + 20);
new_bw6->SetParameters(par_total + 25);
new_bw7->SetParameters(par_total + 30);
new_bw8->SetParameters(par_total + 35);
new_bw9->SetParameters(par_total + 40);
new_bw10->SetParameters(par_total + 45);
new_bw11->SetParameters(par_total + 50);


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

new_bw7->SetLineColor(kBlack);
new_bw7->SetLineWidth(2);
new_bw7->Draw("SAME");

new_bw8->SetLineColor(kBlack);
new_bw8->SetLineWidth(2);
new_bw8->Draw("SAME");

new_bw9->SetLineColor(kBlack);
new_bw9->SetLineWidth(2);
new_bw9->Draw("SAME");

new_bw10->SetLineColor(kBlack);
new_bw10->SetLineWidth(2);
new_bw10->Draw("SAME");

new_bw11->SetLineColor(kBlack);
new_bw11->SetLineWidth(2);
new_bw11->Draw("SAME");

// Configurar la leyenda
TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.85); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
legend->AddEntry(new_bw1, "BW individual fits"); // Agregar entrada para la primera función BW ajustada
legend->AddEntry(h, "Phase Space"); // Agregar entrada para la primera función BW ajustada
legend->AddEntry(total, "Total Fit"); // Agregar entrada para el ajuste total
legend->SetBorderSize(0); // Sin borde
legend->Draw(); // Dibujar la leyenda

// Mostrar el Canvas
gPad->Update();
/*
// Calcular la integral de cada función Breit-Wigner ajustada
Double_t integral_bw1 = new_bw1->Integral(11.3, 14.45) * 15;
Double_t integral_bw2 = new_bw2->Integral(11.3, 14.45) * 15;
Double_t integral_bw3 = new_bw3->Integral(11.3, 14.45) * 15;
Double_t integral_bw4 = new_bw4->Integral(11.3, 14.45) * 15;
Double_t integral_bw5 = new_bw5->Integral(11.3, 14.45) * 15;
Double_t integral_bw6 = new_bw6->Integral(11.3, 14.45) * 15;
Double_t integral_bw7 = new_bw7->Integral(11.3, 14.45) * 15;
Double_t integral_bw8 = new_bw8->Integral(11.3, 14.45) * 15;
Double_t integral_bw9 = new_bw9->Integral(11.3, 14.45) * 15;
Double_t integral_bw10 = new_bw10->Integral(11.3, 14.45) * 15;

// Imprimir los resultados
cout << "Integral of new_bw1: " << integral_bw1 << endl;
cout << "Integral of new_bw2: " << integral_bw2 << endl;
cout << "Integral of new_bw3: " << integral_bw3 << endl;
cout << "Integral of new_bw4: " << integral_bw4 << endl;
cout << "Integral of new_bw5: " << integral_bw5 << endl;
cout << "Integral of new_bw6: " << integral_bw6 << endl;
cout << "Integral of new_bw7: " << integral_bw7 << endl;
cout << "Integral of new_bw8: " << integral_bw8 << endl;
cout << "Integral of new_bw9: " << integral_bw9 << endl;
cout << "Integral of new_bw10: " << integral_bw10 << endl;
*/
}

