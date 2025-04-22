#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCutG.h"
#include "binContents_proton.C"
#include "penetrabilities_L_0.C"
#include "penetrabilities_L_1.C"
#include "penetrabilities_L_2.C"
#include "penetrabilities_L_3.C"
#include "penetrabilities_L_4.C"
#include "penetrabilities_L_5.C"


Float_t Ex = 0;
Int_t detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x[24];
Float_t thetaCM = 0;

const double r = 4.1 * 1.0e-15;

double BWModificada(double *x, double *par, int l, int bw_index) {
  static const double E_min_values[10] = {11.23, 11.23, 11.23, 11.23, 11.23, 11.23, 11.23, 11.23, 11.23, 11.23};
  static const double E_max_values[10] = {14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5};

  if (bw_index < 0 || bw_index > 10) {
    return 0.0;
  }

  double E = x[0];
  double Amp = par[0];
  double E0 = par[1];
  double Gamma0 = par[2];
  double sigma = par[3];

  int num_bins_bw = (E_max_values[bw_index] - E_min_values[bw_index]) / 0.000327;
  int num_bins_pen = 10000;
  double bin_width_pen = (E_max_values[bw_index] - E_min_values[bw_index]) / num_bins_pen;

  int bin_index_bw = (E - E_min_values[bw_index]) / 0.000327;
  if (bin_index_bw < 0 || bin_index_bw >= num_bins_bw) {
    return 0.0;
  }

  int bin_index_pen = (E - E_min_values[bw_index]) / bin_width_pen;
  bin_index_pen = std::min(bin_index_pen, num_bins_pen - 1);

  double E_bin = E_min_values[bw_index] + (bin_index_bw + 0.5) * 0.000327;

  double Gamma_eff;

  if (E_bin >= 11.23) {
    double Gamma_pen_l0 = Gamma0 * T0_values[bin_index_pen];
    double Gamma_pen_l1 = Gamma0 * T1_values[bin_index_pen];
    double Gamma_pen_l2 = Gamma0 * T2_values[bin_index_pen];
    double Gamma_pen_l3 = Gamma0 * T3_values[bin_index_pen];
    double Gamma_pen_l4 = Gamma0 * T4_values[bin_index_pen];
    double Gamma_pen_l5 = Gamma0 * T5_values[bin_index_pen];

    if (l == 0) {
      Gamma_eff = Gamma_pen_l0;
    } else if (l == 1) {
      Gamma_eff = Gamma_pen_l1;
    } else if (l == 2) {
      Gamma_eff = Gamma_pen_l2;
    } else if (l == 3) {
      Gamma_eff = Gamma_pen_l3;
    } else if (l == 4) {
      Gamma_eff = Gamma_pen_l4;
    } else if (l == 5) {
      Gamma_eff = Gamma_pen_l5;
    }
    else {
      Gamma_eff = Gamma0;  
  }
  }

  double denominator = (E_bin * E_bin - E0 * E0) * (E_bin * E_bin - E0 * E0) + (Gamma_eff * Gamma_eff * E0 * E0);
  double gaussian = exp(-0.5 * pow((E - E_bin) / sigma, 2)) / (sigma * sqrt(2 * M_PI));
  double smeared_bw = Amp * gaussian / denominator;
  return smeared_bw;
}



double TotalFunction(double *x, double *par) {
  double E = x[0];
  double total = 0.0;
  // Número de funciones BW
  const int numBW = 10;
  int l_values[10] = {4, 0, 2, 3, 4, 1, 2, 2, 4, 2};
  for (int i = 0; i < numBW; i++) {
    total += BWModificada(
        &E, par + i * 4, l_values[i], i); 
  }

  int bin_index = static_cast<int>((E - 11) / (14.5 - 11) * 15);
  if (bin_index >= 0 && bin_index < 15) {
    total += binContents[bin_index];
  }
  return total;
}

void readData_BW_proton_penetrability() {
  TFile *f = new TFile("h082_10BDP_trace_run013_015-019_025-041.root");
  TTree *tree = (TTree *)f->Get("tree");
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
  auto cutProtonRecoil1 = new TCutG("CUTPROTONRECOIL1", 6);
  cutProtonRecoil1->SetVarX("rdtH[0]");
  cutProtonRecoil1->SetVarY("");
  cutProtonRecoil1->SetTitle("Graph");
  cutProtonRecoil1->SetFillStyle(1000);
  cutProtonRecoil1->SetPoint(0, 3220.905, 1444.03);
  cutProtonRecoil1->SetPoint(1, 3311.691, 1170.896);
  cutProtonRecoil1->SetPoint(2, 3262.632, 1016.418);
  cutProtonRecoil1->SetPoint(3, 3164.516, 1311.94);
  cutProtonRecoil1->SetPoint(4, 3219.777, 1448.507);
  cutProtonRecoil1->SetPoint(5, 3220.905, 1444.03);

  auto cutProtonRecoil2 = new TCutG("CUTPROTONRECOIL2", 6);
  cutProtonRecoil2->SetVarX("rdtH[1]");
  cutProtonRecoil2->SetVarY("");
  cutProtonRecoil2->SetTitle("Graph");
  cutProtonRecoil2->SetFillStyle(1000);
  cutProtonRecoil2->SetPoint(0, 3119.599, 1380.306);
  cutProtonRecoil2->SetPoint(1, 3046.031, 1344.799);
  cutProtonRecoil2->SetPoint(2, 3112.911, 1077.659);
  cutProtonRecoil2->SetPoint(3, 3197.625, 1119.928);
  cutProtonRecoil2->SetPoint(4, 3120.713, 1381.996);
  cutProtonRecoil2->SetPoint(5, 3119.599, 1380.306);

  auto cutProtonRecoil3 = new TCutG("CUTPROTONRECOIL3", 6);
  cutProtonRecoil3->SetVarX("rdtH[2]");
  cutProtonRecoil3->SetVarY("");
  cutProtonRecoil3->SetTitle("Graph");
  cutProtonRecoil3->SetFillStyle(1000);
  cutProtonRecoil3->SetPoint(0, 3261.741, 1338.105);
  cutProtonRecoil3->SetPoint(1, 3168.548, 1306.482);
  cutProtonRecoil3->SetPoint(2, 3249.519, 1100.931);
  cutProtonRecoil3->SetPoint(3, 3348.823, 1144.413);
  cutProtonRecoil3->SetPoint(4, 3264.797, 1342.058);
  cutProtonRecoil3->SetPoint(5, 3261.741, 1338.105);

  auto cutProtonRecoil4 = new TCutG("CUTPROTONRECOIL4", 6);
  cutProtonRecoil4->SetVarX("rdtH[3]");
  cutProtonRecoil4->SetVarY("");
  cutProtonRecoil4->SetTitle("Graph");
  cutProtonRecoil4->SetFillStyle(1000);
  cutProtonRecoil4->SetPoint(0, 3196.594, 1393.424);
  cutProtonRecoil4->SetPoint(1, 3121.597, 1296.035);
  cutProtonRecoil4->SetPoint(2, 3229.732, 1090.438);
  cutProtonRecoil4->SetPoint(3, 3311.706, 1178.358);
  cutProtonRecoil4->SetPoint(4, 3194.85, 1389.366);
  cutProtonRecoil4->SetPoint(5, 3196.594, 1393.424);

  auto cutBoronRecoil1 = new TCutG("cutBoronRecoil1", 8);
  cutBoronRecoil1->SetVarX("rdtH[0]");
  cutBoronRecoil1->SetVarY("");
  cutBoronRecoil1->SetTitle("Graph");
  cutBoronRecoil1->SetFillStyle(1000);
  cutBoronRecoil1->SetPoint(0, 65.059, 4708.98);
  cutBoronRecoil1->SetPoint(1, 410.621, 3944.01);
  cutBoronRecoil1->SetPoint(2, 2157.11, 2218.75);
  cutBoronRecoil1->SetPoint(3, 3333.89, 1795.57);
  cutBoronRecoil1->SetPoint(4, 2913.61, 2414.06);
  cutBoronRecoil1->SetPoint(5, 46.38, 5083.33);
  cutBoronRecoil1->SetPoint(6, 93.0775, 4660.16);
  cutBoronRecoil1->SetPoint(7, 65.059, 4708.98);

  auto cutBoronRecoil2 = new TCutG("CUTBORONRECOIL2", 8);
  cutBoronRecoil2->SetVarX("rdtH[1]");
  cutBoronRecoil2->SetVarY("");
  cutBoronRecoil2->SetTitle("Graph");
  cutBoronRecoil2->SetFillStyle(1000);
  cutBoronRecoil2->SetPoint(0, 32.3707, 4660.16);
  cutBoronRecoil2->SetPoint(1, 1302.54, 2739.58);
  cutBoronRecoil2->SetPoint(2, 3067.71, 1567.71);
  cutBoronRecoil2->SetPoint(3, 3086.39, 1990.89);
  cutBoronRecoil2->SetPoint(4, 2768.85, 2479.17);
  cutBoronRecoil2->SetPoint(5, 23.0312, 5115.89);
  cutBoronRecoil2->SetPoint(6, 41.7102, 4578.78);
  cutBoronRecoil2->SetPoint(7, 32.3707, 4660.16);

  auto cutBoronRecoil3 = new TCutG("CUTBORONRECOIL3", 9);
  cutBoronRecoil3->SetVarX("rdtH[2]");
  cutBoronRecoil3->SetVarY("");
  cutBoronRecoil3->SetTitle("Graph");
  cutBoronRecoil3->SetFillStyle(1000);
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
  cutBoronRecoil4->SetVarX("rdtH[3]");
  cutBoronRecoil4->SetVarY("");
  cutBoronRecoil4->SetTitle("Graph");
  cutBoronRecoil4->SetFillStyle(1000);
  cutBoronRecoil4->SetPoint(0, 41.7102, 4692.71);
  cutBoronRecoil4->SetPoint(1, 1106.41, 3097.66);
  cutBoronRecoil4->SetPoint(2, 2563.38, 1828.12);
  cutBoronRecoil4->SetPoint(3, 3329.22, 1665.36);
  cutBoronRecoil4->SetPoint(4, 3198.46, 1974.61);
  cutBoronRecoil4->SetPoint(5, 2815.54, 2511.72);
  cutBoronRecoil4->SetPoint(6, 41.7102, 5115.89);
  cutBoronRecoil4->SetPoint(7, 41.7102, 4660.16);
  cutBoronRecoil4->SetPoint(8, 41.7102, 4692.71);

  // Histograms
  TH1F *coinTimeH = new TH1F("coinTimeH", "coinTimeH", 1000, -1000, 1000);
  TH1F *eH[24];
  for (auto i = 0; i < 24; ++i)
    eH[i] = new TH1F(Form("eH[%i]", i), Form("eH[%i]", i), 1000, -2, 18);
  TH2F *rdtH[4];
  for (auto i = 0; i < 4; ++i)
    rdtH[i] = new TH2F(Form("rdtH[%i]", i), Form("rdtH[%i]", i), 1000, 0, 6000,
                       1000, 0, 6000);
  TH1F *exH[24];
  for (auto i = 0; i < 24; i++)
    exH[i] = new TH1F(Form("exH[%i]", i), Form("exH[%i]", i), 1000, -2, 18);
  TH1F *exTotalH = new TH1F("exTotalH", "exTotalH", 15, 11, 14.5);

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
  std::cout << " Number of entries : " << nentries << std::endl;

  for (Long64_t i = 0; i < nentries; i++) {
    tree->GetEntry(i);

    if (i % 1000000 == 0) {
      std::cout << " Entry number : " << i << std::endl;
    }



    if (x[detID] < -0.95 || x[detID] > 0.95 || thetaCM < 10 || thetaCM > 45 || e[detID] < 1)
    continue;

    if (!cutProtonRecoil1->IsInside(rdt[0], rdt[1]) &&
        !cutProtonRecoil2->IsInside(rdt[2], rdt[3]) &&
        !cutProtonRecoil3->IsInside(rdt[4], rdt[5]) &&
        !cutProtonRecoil4->IsInside(rdt[6], rdt[7]))
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

exTotalH->Fill(Ex);

  } // events

 /* 
  // 35 bines
  double binEfficiency10k[] = {
            1,      1,      0.1297,      0.1297, 0.1546, 0.1623, 0.1641,
      0.1622, 0.1539, 0.1495, 0.1461, 0.1452, 0.1410, 0.1390, 0.1317,
      0.1256, 0.1258, 0.1216, 0.1171, 0.1167, 0.1117, 0.1106, 0.1074,
      0.1062, 0.1009, 0.1007, 0.0951, 0.0931, 0.0900, 0.0861, 0.0821,
      0.0792, 0.0765, 0.0763, 0.0801};
// 60 bines
      double binEfficiency10k[] = {
      1,1,1,0.1333, 0.1333, 0.1275, 0.1406, 0.1525, 0.1594, 0.1651, 0.1644, 0.1616, 0.1601, 0.1570, 0.1507, 0.1522, 0.1500, 0.1477, 0.1453, 0.1432, 0.1431, 0.1415, 0.1368, 0.1330, 0.1327, 0.1300, 
    0.1246, 0.1246, 0.1241, 0.1216, 0.1184, 0.1167, 0.1181, 0.1135, 0.1120, 0.1189, 0.1108, 0.1094, 0.1067, 0.1071, 0.1043, 0.1022, 0.0992, 0.0998, 0.0944, 0.0948, 0.0914, 0.0917, 0.0910, 0.0883, 
  0.0849, 0.0844, 0.0840, 0.0802, 0.0785, 0.0756, 0.0780, 0.0781, 0.0770, 0.0798, 0.0757};
*/

// 15 bines

double binEfficiency10k[] = {
  1, 0.1337, 0.1634, 0.1622, 0.1510, 0.1431, 0.1317, 0.1231, 0.1164, 0.1106, 0.1036, 0.0938, 0.0900, 0.0821, 0.0756};

  // Después de llenar el histograma exTotalH
  for (int i = 1; i <= exTotalH->GetNbinsX(); i++) {
    double binCenter = exTotalH->GetBinCenter(i);
    if (binCenter >= 11 && binCenter <= 14.5) {
      int efficiencyIndex =
          static_cast<int>((binCenter - 11) / (14.5 - 11) * 15);
      if (efficiencyIndex >= 0 && efficiencyIndex < 15) {
        double efficiency = binEfficiency10k[efficiencyIndex];
        if (efficiency > 0) {
          double content = exTotalH->GetBinContent(i);
          double error = exTotalH->GetBinError(i);
          exTotalH->SetBinContent(i, content / efficiency);
        }
      }
    }
  }

gROOT->SetBatch(kFALSE);

TH1F *h = new TH1F("", "", 15, 11, 14.5);

for (int i = 1; i <= 15; i++) {
    h->SetBinContent(i, binContents[i-1]);
}

TF1 *bw1 = new TF1("bw1", [=](double *x, double *par) { 
    return BWModificada(x, par, 4, 0); 
}, 11.23, 14.5, 4);

TF1 *bw2 = new TF1("bw2", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 1); 
}, 11.23, 14.5, 4);

TF1 *bw3 = new TF1("bw3", [=](double *x, double *par) { 
    return BWModificada(x, par, 2, 2); 
}, 11.23, 14.5, 4);

TF1 *bw4 = new TF1("bw4", [=](double *x, double *par) { 
    return BWModificada(x, par, 3, 3); 
}, 11.23, 14.5, 4);

TF1 *bw5 = new TF1("bw5", [=](double *x, double *par) { 
    return BWModificada(x, par, 4, 4); 
}, 11.23, 14.5, 4);

TF1 *bw6 = new TF1("bw6", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 5); 
}, 11.23, 14.5, 4);

TF1 *bw7 = new TF1("bw7", [=](double *x, double *par) { 
    return BWModificada(x, par, 2, 6); 
}, 11.23, 14.5, 4);

TF1 *bw8 = new TF1("bw8", [=](double *x, double *par) { 
  return BWModificada(x, par, 2, 7); 
}, 11.23, 14.5, 4);

TF1 *bw9 = new TF1("bw9", [=](double *x, double *par) { 
  return BWModificada(x, par, 4, 8); 
}, 11.23, 14.5, 4);

TF1 *bw10 = new TF1("bw10", [=](double *x, double *par) { 
  return BWModificada(x, par, 2, 9); 
}, 11.23, 14.5, 4);

bw1->SetParameters(300, 11.2, 0.2, 0.0657);  
bw2->SetParameters(300, 11.5, 0.2, 0.0657); 
bw3->SetParameters(300, 11.9, 0.2, 0.0657); 
bw4->SetParameters(300, 12.1, 0.2, 0.0657); 
bw5->SetParameters(300, 12.5, 0.2, 0.0657); 
bw6->SetParameters(300, 13.2, 0.2, 0.0657); 
bw7->SetParameters(300, 13.3, 0.2, 0.0657);
bw8->SetParameters(300, 13.6, 0.2, 0.0657);
bw9->SetParameters(300, 14, 0.2, 0.0657);  
bw10->SetParameters(300, 14.2, 0.2, 0.0657); 

// Definir total dada por la suma de los 16 picos
TF1 *total = new TF1("mstotal", TotalFunction, 11.23, 14.5, 40);

// Ajustar cada función a los datos teniendo en cuenta la anterior
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

Double_t par[40]; // 

bw1->GetParameters(&par[0]);
bw2->GetParameters(&par[4]);
bw3->GetParameters(&par[8]);
bw4->GetParameters(&par[12]);
bw5->GetParameters(&par[16]);
bw6->GetParameters(&par[20]);
bw7->GetParameters(&par[24]);
bw8->GetParameters(&par[28]);
bw9->GetParameters(&par[32]);
bw10->GetParameters(&par[36]);

// Configurar nombres de parámetros
for (int i = 0; i < 10; i++) {
    total->SetParName(i * 4 + 0, Form("Amp%d", i + 1));
    total->SetParName(i * 4 + 1, Form("Mean%d", i + 1));
    total->SetParName(i * 4 + 2, Form("Width%d", i + 1));
    total->SetParName(i * 4 + 3, Form("Sigma%d", i + 1));
}

// Establecer límites de parámetros (si es necesario)
total->SetParLimits(0, 1, 100);
total->SetParLimits(1, 11.25, 11.3);
total->SetParLimits(2, 0.1, 0.3);
total->FixParameter(3, 0.0657);
total->SetParLimits(4, 1, 100);
total->SetParLimits(5, 11.5, 11.55);
total->SetParLimits(6, 0., 0.3);
total->FixParameter(7, 0.0657);
total->SetParLimits(9, 11.82, 11.88);
total->SetParLimits(10, 0.190, 0.290);
total->FixParameter(11, 0.0657);
total->SetParLimits(13, 12, 12.1);
total->SetParLimits(14, 0.21, 0.450);
total->FixParameter(15, 0.0657);
total->SetParLimits(17, 12.45, 12.6);
total->SetParLimits(18, 0.375, 0.590);
total->FixParameter(19, 0.0657);
total->SetParLimits(21, 13.19, 13.21);
total->SetParLimits(22, 0.150, 0.230);
total->FixParameter(23, 0.0657);
total->SetParLimits(25, 13.357, 13.380);
total->SetParLimits(26, 0.144, 0.196);
total->FixParameter(27, 0.0657);
total->SetParLimits(29, 13.57, 13.6);
total->SetParLimits(30, 0.15, 0.19);
total->FixParameter(31, 0.0657);
total->SetParLimits(33, 13.98, 14.02);
total->SetParLimits(34, 0.2, 0.26);
total->FixParameter(35, 0.0657);
total->SetParLimits(37, 14.186, 14.210);
total->SetParLimits(38, 0.160, 0.2);
total->FixParameter(39, 0.0657);

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
Double_t par_total[40];
total->GetParameters(par_total);

TF1 *new_bw1 = new TF1("new_bw1", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 4, 0); 
}, 11.23, 14.5, 4);

TF1 *new_bw2 = new TF1("new_bw2", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 1); 
}, 11.23, 14.5, 4);

TF1 *new_bw3 = new TF1("new_bw3", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 2, 2); 
}, 11.23, 14.5, 4);

TF1 *new_bw4 = new TF1("new_bw4", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 3, 3); 
}, 11.23, 14.5, 4);

TF1 *new_bw5 = new TF1("new_bw5", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 4, 4); 
}, 11.23, 14.5, 4);

TF1 *new_bw6 = new TF1("new_bw6", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 5); 
}, 11.23, 14.5, 4);

TF1 *new_bw7 = new TF1("new_bw7", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 2, 6); 
}, 11.23, 14.5, 4);

TF1 *new_bw8 = new TF1("new_bw8", [=](double *x, double *par_total) { 
  return BWModificada(x, par_total, 2, 7); 
}, 11.23, 14.5, 4);

TF1 *new_bw9 = new TF1("new_bw9", [=](double *x, double *par_total) { 
  return BWModificada(x, par_total, 4, 8); 
}, 11.23, 14.5, 4);

TF1 *new_bw10 = new TF1("new_bw10", [=](double *x, double *par_total) { 
  return BWModificada(x, par_total, 2, 9); 
}, 11.23, 14.5, 4);

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

new_bw1->SetParameters(par_total);
new_bw2->SetParameters(par_total+4);
new_bw3->SetParameters(par_total+8);
new_bw4->SetParameters(par_total+12);
new_bw5->SetParameters(par_total+16);
new_bw6->SetParameters(par_total+20);
new_bw7->SetParameters(par_total+24);
new_bw8->SetParameters(par_total+28);
new_bw9->SetParameters(par_total+32);
new_bw10->SetParameters(par_total+36);

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

// Dibujar solo la línea de ajuste resultante

TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
//legend->AddEntry(new_bw1, "BW individual fits"); // Agregar entrada para la primera función BW ajustada
//legend->AddEntry(total, "Total Fit"); 
legend->AddEntry(h, "Phase Space"); 
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

