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
#include "binContents_alpha.C"
#include "penetrabilities_alpha_L_0.C"
#include "penetrabilities_alpha_L_1.C"
#include "penetrabilities_alpha_L_2.C"
#include "penetrabilities_alpha_L_3.C"
#include "penetrabilities_alpha_L_4.C"
#include "penetrabilities_alpha_L_5.C"

Float_t Ex = 0;
Int_t detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x[24];
Float_t thetaCM = 0;


    double BWModificada(double *x, double *par, int l, int bw_index) {

      // Definir los rangos de energía para cada Breit-Wigner (E_min y E_max)
        static const double E_min_values[16] = {8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5,8.5, 8.5};  // Energías mínimas
        static const double E_max_values[16] = {14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5};  // Energías máximas
      
        if (bw_index < 0 || bw_index >=16){
            return 0.0;
        }

        double E = x[0];
        double Amp = par[0];
        double E0 = par[1];
        double Gamma0 = par[2];
        double sigma = par[3];

        int num_bins_bw = (E_max_values[bw_index] - E_min_values[bw_index]) / 0.000633;
        int num_bins_pen = 10000;
        double bin_width_pen = (E_max_values[bw_index] - E_min_values[bw_index]) / num_bins_pen;
      
        int bin_index_bw = (E - E_min_values[bw_index]) / 0.000633;
        if (bin_index_bw < 0 || bin_index_bw >= num_bins_bw) {
          return 0.0;
        }

          // Cálculo de bin para E0
        int bin_index_E0 = static_cast<int>((E0 - E_min_values[bw_index])/bin_width_pen);
        bin_index_E0 = std::max(0, std::min(bin_index_E0, num_bins_pen - 1));
      
        int bin_index_pen = (E - E_min_values[bw_index]) / bin_width_pen;
        bin_index_pen = std::min(bin_index_pen, num_bins_pen - 1);
      
        double E_bin = E_min_values[bw_index] + (bin_index_bw + 0.5) * 0.000633;
      
        double Gamma_eff;
      
        if (E_bin >= 8.66) {
            double Gamma_pen_l0 = Gamma0 * T0_values_alpha[bin_index_pen] / T0_values_alpha[bin_index_E0];
            double Gamma_pen_l1 = Gamma0 * T1_values_alpha[bin_index_pen] / T1_values_alpha[bin_index_E0];
            double Gamma_pen_l2 = Gamma0 * T2_values_alpha[bin_index_pen] / T2_values_alpha[bin_index_E0];
            double Gamma_pen_l3 = Gamma0 * T3_values_alpha[bin_index_pen] / T3_values_alpha[bin_index_E0];
            double Gamma_pen_l4 = Gamma0 * T4_values_alpha[bin_index_pen] / T4_values_alpha[bin_index_E0];
            double Gamma_pen_l5 = Gamma0 * T5_values_alpha[bin_index_pen] / T5_values_alpha[bin_index_E0];

          if (l == 0) {
            Gamma_eff = Gamma_pen_l0;
          } else if (l == 1) {
            Gamma_eff = Gamma_pen_l1;
          } else if (l == 2) {
            Gamma_eff = Gamma_pen_l2;
          } else if (l == 3) {
            Gamma_eff = Gamma_pen_l3;
          }
          else if (l == 4) {
            Gamma_eff = Gamma_pen_l4;
          } else if (l == 5) {
            Gamma_eff = Gamma_pen_l5;
          }
          else {
            Gamma_eff = Gamma0;  
        }
          
        }
      
        double denominator = (E_bin * E_bin - E0 * E0) * (E_bin * E_bin - E0 * E0) + (Gamma_eff * Gamma_eff * E0 * E0);
        double numerator = 2 * sqrt(2) * E0 * Gamma_eff * sqrt(E0*E0*((E0*E0) + (Gamma_eff*Gamma_eff))) / (M_PI * sqrt((E0*E0) + (sqrt(E0*E0*((E0*E0) + (Gamma_eff*Gamma_eff))))));
        double gaussian = exp(-0.5 * pow((E - E_bin) / sigma, 2)) / (sigma * sqrt(2 * M_PI));
        double smeared_bw = Amp*numerator * gaussian / denominator;
        return smeared_bw;
}

double TotalFunction(double *x, double *par) {
    double E = x[0];
    double total = 0.0;

    // Número de funciones BW
    const int numBW = 16;

    int l_values[16] = {3, 1, 1, 0, 3, 3, 1, 1, 2, 3, 1, 2, 1, 1, 1, 1};
  for (int i = 0; i < numBW; i++) {
    total += BWModificada(
        &E, par + i * 4, l_values[i],
        i); // Pass the index 'i' to BWModificada to access E_min/max
  }
// Sumar el contenido del bin correspondiente
    int bin_index = static_cast<int>((E - 8.5) / (14.5 - 8.5) * 90);
    if (bin_index >= 0 && bin_index < 90) {
        total += binContents[bin_index];
    }


    return total;
}

void readData_BW_alpha_def_penetrability()
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
   cutLitium1->Draw("");

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
   cutLitium3->Draw("");

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
   cutLitium2->Draw("");

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
   cutLitium4->Draw("");

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
   cutalpha1->Draw("");

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
   cutalpha3->Draw("");

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
   cutalpha2->Draw("");

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
   cutalpha4->Draw("");

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

      if (x[detID] < -0.95 || x[detID] > 0.95 || thetaCM < 10 || e[detID] < 1)
        continue;

      bool passLitiumCut = cutLitium1->IsInside(rdt[0], rdt[1]) || cutLitium2->IsInside(rdt[2], rdt[3]) || cutLitium3->IsInside(rdt[4], rdt[5]) || cutLitium4->IsInside(rdt[6], rdt[7]);
      bool passAlphaCut = cutalpha1->IsInside(rdt[0], rdt[1]) || cutalpha2->IsInside(rdt[2], rdt[3]) || cutalpha3->IsInside(rdt[4], rdt[5]) || cutalpha4->IsInside(rdt[6], rdt[7]);

        if (!passLitiumCut && !passAlphaCut)
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

  }//events
/*
  double binEfficiency10k[] = {0.0005, 0.0007, 0.0008, 1,1,1, 0.1951, 0.1931, 0.1837, 0.1783, 0.1725, 
0.1731, 0.1743, 0.171, 0.1721, 0.1723, 0.1697, 0.1658, 0.1615, 0.1544, 0.1518, 0.1526,
  0.1484, 0.1437, 0.1409, 0.1355, 0.1335, 0.1264, 0.1207, 0.1154, 0.1122, 0.1106, 0.1072, 
   0.105, 0.1022, 0.0987, 0.0907, 0.0874, 0.0866, 0.0832, 0.0819, 0.0791, 0.0776, 0.0761, 
   0.0755, 0.0723, 0.0683, 0.0667, 0.0631, 0.0616, 0.0599, 0.0573, 0.0576, 0.0535, 0.0517, 
   0.0501, 0.0486, 0.0479, 0.046, 0.0458, 0.0435, 0.0414, 0.0394, 0.0370, 0.0360, 0.0336, 
   0.0342, 0.0336, 0.0328, 0.0311, 0.0305, 0.0294, 0.0276, 0.0287, 0.0269, 0.0265, 0.0251, 
   0.0255, 0.0236, 0.0228, 0.0209, 0.0198, 0.0183, 0.0188, 0.0192, 0.0178, 0.0177, 0.0183,
   0.0165, 0.0164, 0.0167, 0.0164, 0.0149
   };
   */
    // Para el GS sin coincidencia ya hechas las cuentas.
   double binEfficiency10k[] = {0.0441,0.0441,0.0441, 0.0441, 0.1095, 0.1422, 0.1501, 0.1561, 0.1591, 0.1604, 0.1593, 0.1570, 0.1532, 0.1531, 0.1527, 
    0.1527, 0.1503, 0.1464, 0.1443, 0.1394, 0.1330, 0.1277, 0.1256, 0.1212, 0.1168, 0.1153, 0.1132, 0.1085, 0.1053, 0.1024, 0.0994, 
    0.0970, 0.0963, 0.0910, 0.0880, 0.0872, 0.0850, 0.0840, 0.0816, 0.0770, 0.0749, 0.0715, 0.0696, 0.0676, 0.0654, 0.0618, 0.0626, 
    0.0606, 0.0577, 0.0545, 0.0504, 0.0492, 0.0498, 0.0474, 0.0474, 0.0421, 0.0414, 0.0423, 0.0387, 0.0370, 0.0347, 0.0334, 0.0303,
    0.0314, 0.0297, 0.0291, 0.0278, 0.0261, 0.0270, 0.0264, 0.0252, 0.0241, 0.0223, 0.0226, 0.0212, 0.0195, 0.0193, 0.0190, 0.0171, 
    0.0166, 0.0166, 0.0155, 0.0144, 0.0147, 0.0155, 0.0146, 0.0145, 0.0147, 0.0152, 0.0135
       };

    // Para el GS con coincidencia únicamente.
       double binEfficiency10k_coin[] = {1,1,1, 0.1539, 0.0865, 0.0516, 0.0302, 0.0232, 0.0181, 0.0157, 0.0128, 0.0116, 0.0103, 0.0100, 0.0092, 
    0.0085, 0.0068, 0.0065, 0.0057, 0.0052, 0.0051, 0.0049, 0.0041, 0.0036, 0.0034, 0.0030, 0.0026, 0.0027, 0.0024, 0.0022, 0.0019, 
    0.0016, 0.0021, 0.0020, 0.0019, 0.0018, 0.0019, 0.0019, 0.0019, 0.0019, 0.0020, 0.0018, 0.0018, 0.0016, 0.0016, 0.0017, 0.0017, 
    0.0014, 0.0012, 0.0013, 0.0012, 0.001, 0.001, 0.0006, 0.0006, 0.0009, 0.0009, 0.0009, 0.0011, 0.0012, 0.0010, 0.0010, 0.0011, 
    0.0011, 0.0010, 0.0006, 0.0008, 0.0007, 0.0009, 0.0006, 0.0007, 0.0006, 0.0006, 0.0005, 0.0007, 0.0005, 0.0004, 0.0006, 0.0005, 
    0.0007, 0.0006, 0.0005, 0.0005, 0.0006, 0.0005, 0.0004, 0.0004, 0.0004, 0.0004, 0.0003
    };
/*
        // Para el *Li sin coincidencia ya hechas las cuentas.
   double binEfficiency10k_exc[] = {1,1,1,1,1,1,1,1,1,1, 0.0310, 0.1043, 0.1431, 0.1521, 0.1585, 0.1627, 0.1697, 
0.1678, 0.1704, 0.1654, 0.1601, 0.1579, 0.1544, 0.1491, 0.1445, 0.1379, 0.1336, 0.1310, 0.1247, 0.124, 0.118, 
0.1146, 0.1100, 0.1065, 0.1021, 0.1002, 0.0988, 0.0956, 0.0953, 0.0905, 0.0875, 0.0842, 0.0833, 0.0806, 0.0773, 
0.0735, 0.0718, 0.0705, 0.0663, 0.0651, 0.0630, 0.0609, 0.0593, 0.0574, 0.0563, 0.0560, 0.0530, 0.0492, 0.0508,
0.0468, 0.0441, 0.0428, 0.0404, 0.0388, 0.0368, 0.0355, 0.0340, 0.0327, 0.0308, 0.0301, 0.0296, 0.0282, 0.0291, 
0.0273, 0.0243, 0.0238, 0.0233, 0.0215, 0.0210, 0.0193, 0.0191, 0.0182, 0.0164, 0.0166, 0.0167, 0.0158, 0.0151, 
0.0156, 0.0156, 0.0160, 0.0154
};

// Para el *Li con coincidencia únicamente.
double binEfficiency10k_coin_exc[] = {1,1,1,1,1,1,1,1,1,1, 0.1822, 0.1043, 0.0597, 0.0369, 0.0277, 0.0212, 0.0170, 
0.0142, 0.0120, 0.0107, 0.0101, 0.0082, 0.0078, 0.0076, 0.0064, 0.0055, 0.0049, 0.0046, 0.0041, 0.0041, 0.0034, 
0.0029, 0.0028, 0.0027, 0.0026, 0.0022, 0.0022, 0.0024, 0.0023, 0.0024, 0.0022, 0.0019, 0.0018, 0.0018, 0.0017, 
0.0018, 0.0018, 0.0017, 0.0015, 0.0013, 0.0011, 0.0012, 0.0010, 0.0009, 0.0011, 0.0007, 0.0010, 0.0009, 0.0011,
0.0011, 0.0010, 0.0013, 0.0012, 0.0013, 0.0012, 0.0009, 0.0009, 0.0008, 0.0007, 0.0008, 0.0008, 0.0006, 0.0006, 
0.0007, 0.0006, 0.0006, 0.0007, 0.0005, 0.0008, 0.0006, 0.0006, 0.0006, 0.0005, 0.0005, 0.0006, 0.0006, 0.0006,
0.0006, 0.0007, 0.0003, 0.0004
};
*/
// Después de llenar el histograma exTotalH
for (int i = 1; i <= exTotalH->GetNbinsX(); i++) {
  double binCenter = exTotalH->GetBinCenter(i);
  if (binCenter >= 8.5 && binCenter <= 14.5) {
    int efficiencyIndex = static_cast<int>((binCenter - 8.5) / (14.5 - 8.5) * 90);
    if (efficiencyIndex >= 0 && efficiencyIndex < 90) {
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



gROOT->SetBatch(kFALSE);

TH1F *h = new TH1F("", "", 90, 8.5, 14.5);

for (int i = 1; i <= 90; i++) {
    h->SetBinContent(i, binContents[i-1]);
}

TF1 *bw1 = new TF1("bw1", [=](double *x, double *par) { 
    return BWModificada(x, par, 3, 0); 
}, 8.5, 14.5, 4);

TF1 *bw2 = new TF1("bw2", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 1); 
}, 8.5, 14.5, 4);

TF1 *bw3 = new TF1("bw3", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 2); 
}, 8.5, 14.5, 4);

TF1 *bw4 = new TF1("bw4", [=](double *x, double *par) { 
    return BWModificada(x, par, 0, 3); 
}, 8.5, 14.5, 4);

TF1 *bw5 = new TF1("bw5", [=](double *x, double *par) { 
    return BWModificada(x, par, 3, 4); 
}, 8.5, 14.5, 4);

TF1 *bw6 = new TF1("bw6", [=](double *x, double *par) { 
    return BWModificada(x, par, 3, 5); 
}, 8.5, 14.5, 4);

TF1 *bw7 = new TF1("bw7", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 6); 
}, 8.5, 14.5, 4);

TF1 *bw8 = new TF1("bw8", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 7); 
}, 8.5, 14.5, 4);

TF1 *bw9 = new TF1("bw9", [=](double *x, double *par) { 
    return BWModificada(x, par, 2, 8); 
}, 8.5, 14.5, 4);

TF1 *bw10 = new TF1("bw10", [=](double *x, double *par) { 
    return BWModificada(x, par, 3, 9); 
}, 8.5, 14.5, 4);

TF1 *bw11 = new TF1("bw11", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 10); 
}, 8.5, 14.5, 4);

TF1 *bw12 = new TF1("bw12", [=](double *x, double *par) { 
    return BWModificada(x, par, 2, 11); 
}, 8.5, 14.5, 4);

TF1 *bw13 = new TF1("bw13", [=](double *x, double *par) { 
    return BWModificada(x, par, 1,12); 
}, 8.5, 14.5, 4);

TF1 *bw14 = new TF1("bw14", [=](double *x, double *par) { 
    return BWModificada(x, par, 1,13); 
}, 8.5, 14.5, 4);

TF1 *bw15 = new TF1("bw15", [=](double *x, double *par) { 
    return BWModificada(x, par, 1,14); 
}, 8.5, 14.5, 4);

TF1 *bw16 = new TF1("bw16", [=](double *x, double *par) { 
    return BWModificada(x, par, 1, 15); 
}, 8.5, 14.5, 4);



// Definir parámetros iniciales para cada Breit-Wigner
bw1->SetParameters(300, 9.18, 0.004, 0.0657);  // Parámetros iniciales: amplitud, media, anchura (ancho a media altura)
bw2->SetParameters(300, 9.27, 0.004, 0.0657); 
bw3->SetParameters(300, 9.85, 0.10, 0.0657);
bw4->SetParameters(300, 10.28, 0.10, 0.0657);
bw5->SetParameters(300, 10.6, 0.2, 0.0657);
bw6->SetParameters(300, 11.2, 0.2, 0.0657);
bw7->SetParameters(300, 11.45, 0.2, 0.0657);
bw8->SetParameters(300, 11.66, 0.2, 0.0657);
bw9->SetParameters(300, 11.84, 0.2, 0.0657);
bw10->SetParameters(300, 12., 0.2, 0.0657);
bw11->SetParameters(300, 12.5, 0.2, 0.0657);
bw12->SetParameters(300, 12.95, 0.2, 0.0657);
bw13->SetParameters(300, 13.2, 0.2, 0.0657);
bw14->SetParameters(300, 13.3, 0.2, 0.0657);
bw15->SetParameters(300, 13.6, 0.2, 0.0657);
bw16->SetParameters(300, 13.8, 0.2, 0.0657);

// Definir total dada por la suma de los 16 picos
TF1 *total = new TF1("mstotal", TotalFunction, 8.5, 14.5, 68);

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
exTotalH->Fit(bw11, "R+");
exTotalH->Fit(bw12, "R+");
exTotalH->Fit(bw13, "R+");
exTotalH->Fit(bw14, "R+");
exTotalH->Fit(bw15, "R+");
exTotalH->Fit(bw16, "R+");


// Crear un array para almacenar los parámetros
Double_t par[64]; // 

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
bw11->GetParameters(&par[40]);
bw12->GetParameters(&par[44]);
bw13->GetParameters(&par[48]);
bw14->GetParameters(&par[52]);
bw15->GetParameters(&par[56]);
bw16->GetParameters(&par[60]);


total->SetParameters(par);

// Configurar nombres de parámetros
for (int i = 0; i < 16; i++) {
    total->SetParName(i * 4 + 0, Form("Amp%d", i + 1));
    total->SetParName(i * 4 + 1, Form("Mean%d", i + 1));
    total->SetParName(i * 4 + 2, Form("Width%d", i + 1));
    total->SetParName(i * 4 + 3, Form("Sigma%d", i + 1));
}


// Establecer límites de parámetros (si es necesario)
total->SetParLimits(0, 1, 10000);
total->SetParLimits(1, 9.15, 9.25);
total->SetParLimits(2, 0.05, 0.3);
total->FixParameter(3, 0.0657);
total->SetParLimits(4, 1, 10000);
total->SetParLimits(5, 9.25, 9.40);
total->SetParLimits(6, 0.05, 0.3);
total->FixParameter(7, 0.0657);
total->SetParLimits(8, 1, 10000);
total->SetParLimits(9, 9.85, 9.90);
total->SetParLimits(10, 0.05, 0.3);
total->FixParameter(11, 0.0657);
total->SetParLimits(12, 1, 10000);
total->SetParLimits(13, 10.28, 10.38);
total->SetParLimits(14, 0.05, 0.3);
total->FixParameter(15, 0.0657);
total->SetParLimits(16, 1, 10000);
total->SetParLimits(17, 10.6, 10.7);
total->SetParLimits(18, 0.01, 0.3);
total->FixParameter(19, 0.0657);
total->SetParLimits(20, 1, 10000);
total->SetParLimits(21, 11.27, 11.31);
total->SetParLimits(22, 0.01, 0.3);
total->FixParameter(23, 0.0657);
total->SetParLimits(24, 1, 10000);
total->SetParLimits(25, 11.4, 11.5);
total->SetParLimits(26, 0.1, 0.15);
total->FixParameter(27, 0.0657);
total->SetParLimits(28, 1, 10000);
total->SetParLimits(29, 11.67, 11.71);
total->SetParLimits(30, 0.150, 0.250);
total->FixParameter(31, 0.0657);
total->SetParLimits(32, 1, 10000);
total->SetParLimits(33, 11.82, 11.88);
total->SetParLimits(34, 0.190, 0.290);
total->FixParameter(35, 0.0657);
total->SetParLimits(36, 1, 10000);
total->SetParLimits(37, 12, 12.1);
total->SetParLimits(38, 0.208, 0.452);
total->FixParameter(39, 0.0657);
total->SetParLimits(40, 1, 10000);
total->SetParLimits(41, 12.45, 12.52);
total->SetParLimits(42, 0.375, 0.585);
total->FixParameter(43, 0.0657);
total->SetParLimits(44, 1, 10000);
total->SetParLimits(45, 12.85, 12.91);
total->SetParLimits(46, 0.264, 0.556);
total->FixParameter(47, 0.0657);
total->SetParLimits(48, 1, 10000);
total->SetParLimits(49, 13.19, 13.21);
total->SetParLimits(50, 0.150, 0.230);
total->FixParameter(51, 0.0657);
total->SetParLimits(52, 1, 10000);
total->SetParLimits(53, 13.36, 13.38);
total->SetParLimits(54, 0.144, 0.196);
total->FixParameter(55, 0.0657);
total->SetParLimits(56, 1, 10000);
total->SetParLimits(57, 13.57, 13.6);
total->SetParLimits(58, 0.175, 0.25);
total->FixParameter(59, 0.0657);
total->SetParLimits(60, 1, 10000);
total->SetParLimits(61, 13.79, 13.82);
total->SetParLimits(62, 0.160, 0.240);
total->FixParameter(63, 0.0657);

// Ajustar solo la función total y desactivar la visualización de la línea de ajuste resultante
exTotalH->Fit(total, "R+");

// Configurar opciones de visualización para la línea de ajuste resultante
total->SetLineColor(kRed);
total->SetLineWidth(2);
total->SetNpx(10000);

exTotalH->SetTitle("Energy distribution ^{7}Li and #alpha"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");
exTotalH->SetStats(0);

// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
h->SetFillColor(kGreen);
h->Draw("HIST SAME");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[64];
total->GetParameters(par_total);

TF1 *new_bw1 = new TF1("new_bw1", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 3, 0); 
}, 8.5, 14.5, 4);

TF1 *new_bw2 = new TF1("new_bw2", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 1); 
}, 8.5, 14.5, 4);

TF1 *new_bw3 = new TF1("new_bw3", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 2); 
}, 8.5, 14.5, 4);

TF1 *new_bw4 = new TF1("new_bw4", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 0, 3); 
}, 8.5, 14.5, 4);

TF1 *new_bw5 = new TF1("new_bw5", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 3, 4); 
}, 8.5, 14.5, 4);

TF1 *new_bw6 = new TF1("new_bw6", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 3, 5); 
}, 8.5, 14.5, 4);

TF1 *new_bw7 = new TF1("new_bw7", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 6); 
}, 8.5, 14.5, 4);

TF1 *new_bw8 = new TF1("new_bw8", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 7); 
}, 8.5, 14.5, 4);

TF1 *new_bw9 = new TF1("new_bw9", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 2, 8); 
}, 8.5, 14.5, 4);

TF1 *new_bw10 = new TF1("new_bw10", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 3, 9); 
}, 8.5, 14.5, 4);

TF1 *new_bw11 = new TF1("new_bw11", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 10); 
}, 8.5, 14.5, 4);

TF1 *new_bw12 = new TF1("new_bw12", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 2, 11); 
}, 8.5, 14.5, 4);

TF1 *new_bw13 = new TF1("new_bw13", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 12); 
}, 8.5, 14.5, 4); 

TF1 *new_bw14 = new TF1("new_bw14", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 13); 
}, 8.5, 14.5, 4);

TF1 *new_bw15 = new TF1("new_bw15", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 14); 
}, 8.5, 14.5, 4);

TF1 *new_bw16 = new TF1("new_bw16", [=](double *x, double *par_total) { 
    return BWModificada(x, par_total, 1, 15); 
}, 8.5, 14.5, 4);

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
new_bw12->SetNpx(10000);
new_bw13->SetNpx(10000);
new_bw14->SetNpx(10000);
new_bw15->SetNpx(10000);
new_bw16->SetNpx(10000);


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
new_bw11->SetParameters(par_total+40);
new_bw12->SetParameters(par_total+44);
new_bw13->SetParameters(par_total+48);
new_bw14->SetParameters(par_total+52);
new_bw15->SetParameters(par_total+56);
new_bw16->SetParameters(par_total+60);

// Dibujar las nuevas funciones en el mismo Canvas
new_bw1->SetLineColor(kBlack);
new_bw2->SetLineColor(kBlack);
new_bw3->SetLineColor(kBlack);
new_bw4->SetLineColor(kBlack);
new_bw5->SetLineColor(kBlack);
new_bw6->SetLineColor(kBlack);
new_bw7->SetLineColor(kBlack);
new_bw8->SetLineColor(kBlack);
new_bw9->SetLineColor(kBlack);
new_bw10->SetLineColor(kBlack);
new_bw11->SetLineColor(kBlack);
new_bw12->SetLineColor(kBlack);
new_bw13->SetLineColor(kBlack);
new_bw14->SetLineColor(kBlack);
new_bw15->SetLineColor(kBlack);
new_bw16->SetLineColor(kBlack);


new_bw1->Draw("SAME");
new_bw2->Draw("SAME");
new_bw3->Draw("SAME");
new_bw4->Draw("SAME");
new_bw5->Draw("SAME");
new_bw6->Draw("SAME");
new_bw7->Draw("SAME");
new_bw8->Draw("SAME");
new_bw9->Draw("SAME");
new_bw10->Draw("SAME");
new_bw11->Draw("SAME");
new_bw12->Draw("SAME");
new_bw13->Draw("SAME");
new_bw14->Draw("SAME");
new_bw15->Draw("SAME");
new_bw16->Draw("SAME");

// Configurar la leyenda
TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
legend->AddEntry(new_bw1, "BW individual fits"); // Agregar entrada para la primera función BW ajustada
legend->AddEntry(total, "Total Fit"); 
legend->AddEntry(h, "Phase Space"); 
legend->SetBorderSize(0); // Sin borde
legend->Draw(); // Dibujar la leyenda

// Mostrar el Canvas
gPad->Update();
/*
// Calcular la integral de cada función Breit-Wigner ajustada
Double_t integral_bw1 = new_bw1->Integral(8.5, 14.5) * 15;
Double_t integral_bw2 = new_bw2->Integral(8.5, 14.5) * 15;
Double_t integral_bw3 = new_bw3->Integral(8.5, 14.5) * 15;
Double_t integral_bw4 = new_bw4->Integral(8.5, 14.5) * 15;
Double_t integral_bw5 = new_bw5->Integral(8.5, 14.5) * 15;
Double_t integral_bw6 = new_bw6->Integral(8.5, 14.5) * 15;
Double_t integral_bw7 = new_bw7->Integral(8.5, 14.5) * 15;
Double_t integral_bw8 = new_bw8->Integral(8.5, 14.5) * 15;
Double_t integral_bw9 = new_bw9->Integral(8.5, 14.5) * 15;
Double_t integral_bw10 = new_bw10->Integral(8.5, 14.5) * 15;
Double_t integral_bw11 = new_bw11->Integral(8.5, 14.5) * 15;
Double_t integral_bw12 = new_bw12->Integral(8.5, 14.5) * 15;
Double_t integral_bw13 = new_bw13->Integral(8.5, 14.5) * 15;
Double_t integral_bw14 = new_bw14->Integral(8.5, 14.5) * 15;
Double_t integral_bw15 = new_bw15->Integral(8.5, 14.5) * 15;



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
cout << "Integral of new_bw11: " << integral_bw11 << endl;
cout << "Integral of new_bw12: " << integral_bw12 << endl;
cout << "Integral of new_bw13: " << integral_bw13 << endl;
cout << "Integral of new_bw14: " << integral_bw14 << endl;
cout << "Integral of new_bw15: " << integral_bw15 << endl;

*/

}

