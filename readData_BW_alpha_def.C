Float_t Ex = 0;
Int_t detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x = 0;
Float_t thetaCM = 0;

double binContents[] = {
    2.53001,    1.95318,    1.70915,    1.39168,    2.59658,    2.04515,    1.50684,    1.6893,
    1.7507,    1.37287,    1.61848,    1.50365,    2.20909,    2.12121,    2.15952,    1.41156,
    1.6826,    2.16727,    2.0052,    1.74005,    1.1083,    1.83284,    1.7132,    1.44595,
    1.54786,    1.97799,    1.92393,    1.48235,    2.02284,    2.1438,    1.53354,    1.27626,
    1.52636,    1.4259,    2.20378,    2.11581,    2.41934,    2.21713,    1.67165,    1.43407,
    1.75284,    2.04898,    1.98122,    2.07589,    1.85864,    1.66806,    1.46558,    1.76489,
    1.84274,    2.13718,    1.91634,    1.38483,    1.29042,    1.10369,    1.84644,    1.37262,
    1.74341,    1.9234,    1.70242,    1.69708,    1.4763,    1.71216,    1.35865,    1.67331,
    1.48487,    1.55526,    1.24905,    1.306,    1.44196,    1.22488,    1.77157,    1.47247,
    1.42804,    1.62576,    1.52399,    1.14709,    1.61432,    1.25628,    1.26631,    1.18114,
    1.52448,    2.0716,    1.38945,    1.71122,    1.43502,    1.00044,    0.916956,    1.732,
    1.49341,    1.19233
    };

double binEfficiency10k[] = { 0.3353, 0.339, 0.3404, 0.3423, 0.3407, 0.3407, 0.3414, 0.3415, 0.3448, 0.3448, 0.3484, 0.3498, 0.3517, 0.3527, 0.3591, 0.3584, 0.3604, 0.3598, 0.356, 0.3557,
   0.3572, 0.3584, 0.3621, 0.3589, 0.3642, 0.3661, 0.3671, 0.369, 0.3689, 0.375,
   0.3723, 0.3762, 0.3765, 0.3793, 0.38, 0.3824, 0.3820, 0.3778, 0.3833, 0.3883,
   0.3903, 0.3895, 0.3912, 0.3894, 0.3919, 0.3923, 0.3868, 0.3849, 0.3812,
   0.3824, 0.3765, 0.3752, 0.3671, 0.3671, 0.3671, 0.3673, 0.3665, 0.3668,
   0.3609, 0.3592, 0.3520, 0.3420, 0.3441, 0.3412, 0.3434, 0.3478, 0.349, 0.349,
   0.3444, 0.3397, 0.3373, 0.3358, 0.3319, 0.3297, 0.3265, 0.3281, 0.3234, 0.3212,
   0.3192, 0.3124, 0.3144, 0.3078, 0.3034, 0.3023, 0.2956, 0.2987, 0.3038, 0.3025, 0.3023, 0.3014, 0.2993
                                
    };


    double BW(double *x, double *par, int bw_index) {

      // Definir los rangos de energía para cada Breit-Wigner (E_min y E_max)
        static const double E_min_values[17] = {8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5};  // Energías mínimas
        static const double E_max_values[17] = {14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5, 14.5};  // Energías máximas
      
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
    double denominator = (E * E - E0 * E0) * (E * E - E0 * E0) + Gamma0 * Gamma0 * E0 * E0;

    double gaussian = exp(-0.5 * pow((E - E_bin) / sigma, 2)) / (sigma * sqrt(2 * M_PI));  // Distribución Gaussiana

    return Amp *gaussian / denominator  ;
}

double TotalFunction(double *x, double *par) {
    double E = x[0];
    double total = 0.0;

    // Número de funciones BW
    const int numBW = 17;

    for (int i = 0; i < numBW; i++) {
        total += BW(&E, par + i * 4, i);
        total += binContents[i];
    }

    return total;
}

void readData_BW_alpha_def()
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

      bool passLitiumCut = cutLitium1->IsInside(rdt[0], rdt[1]) || cutLitium2->IsInside(rdt[2], rdt[3]) || cutLitium3->IsInside(rdt[4], rdt[5]) || cutLitium4->IsInside(rdt[6], rdt[7]);
      bool passAlphaCut = cutalpha1->IsInside(rdt[0], rdt[1]) || cutalpha2->IsInside(rdt[2], rdt[3]) || cutalpha3->IsInside(rdt[4], rdt[5]) || cutalpha4->IsInside(rdt[6], rdt[7]);

        if (!passLitiumCut && !passAlphaCut)
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

  TCanvas *c4 = new TCanvas();
  coinTimeH->Draw();

*/
/*
TH1F *h = new TH1F("", "", 90, 8.5, 14.5);

for (int i = 1; i <= 90; i++) {
    h->SetBinContent(i, binContents[i-1]);
}
*/

const int numFunctions = 17;

    // Crear el array de funciones Breit-Wigner
    TF1 *bw[numFunctions];

    // Valores iniciales de los parámetros para cada Breit-Wigner
    double amplitudes[numFunctions] = {300, 100, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300};
    double E0_values[numFunctions] = {9.18, 9.27, 9.85, 10.28, 10.6, 11.2, 11.45, 11.64, 11.85, 12.0, 12.55, 12.87, 13.3, 13.6, 13.85, 14, 14.2};
    double widths[numFunctions] = {0.004, 0.004, 0.10, 0.10, 0.2, 0.2, 0.2, 0.35, 0.4, 0.4, 1.0, 0.4, 0.6, 0.3, 0.5, 0.4, 0.2};
    double sigma = 0.0657;  // Resolución común para todas

    // Crear las funciones en el bucle
    for (int i = 0; i < numFunctions; i++) {
        // Usar lambdas para pasar el índice `i` y crear cada TF1 dinámicamente
        bw[i] = new TF1(Form("bw%d", i + 1), [i](double *x, double *par) {
            return BW(x, par, i);
        }, 8.5, 14.5, 4); // Intervalo de energía y 4 parámetros

        // Configurar parámetros iniciales para cada Breit-Wigner
        bw[i]->SetParameters(amplitudes[i], E0_values[i], widths[i], sigma);
    }


// Definir total dada por la suma de los 8 picos
TF1 *total = new TF1("mstotal", TotalFunction, 8.5, 14.5, 68);

// Ajustar todas las funciones Breit-Wigner de manera dinámica
for (int i = 0; i < 17; i++) {
    if (i == 0) {
        // Para la primera función, sin el "+" (sin conservar ajustes previos)
        exTotalH->Fit(bw[i], "R");
    } else {
        // Para las siguientes funciones, incluir el "+" para sumar los ajustes
        exTotalH->Fit(bw[i], "R+");
    }
}


// Crear un array para almacenar los parámetros
Double_t par[68]; // 

// Obtener los parámetros de todas las funciones Breit-Wigner de manera dinámica
for (int i = 0; i < 17; i++) {
    bw[i]->GetParameters(&par[i * 4]); // Guardar los 4 parámetros de cada función en el lugar correspondiente
}


total->SetParameters(par);

// Configurar nombres de parámetros
for (int i = 0; i < 17; i++) {
    total->SetParName(i * 4 + 0, Form("Amp%d", i + 1));
    total->SetParName(i * 4 + 1, Form("Mean%d", i + 1));
    total->SetParName(i * 4 + 2, Form("Width%d", i + 1));
    total->SetParName(i * 4 + 3, Form("Sigma%d", i + 1));
}


// Establecer límites de parámetros (si es necesario)
total->SetParLimits(1, 9.17, 9.19);
total->SetParLimits(2, 0., 0.3);
total->FixParameter(3, 0.0657);
total->SetParLimits(5, 9.25, 9.30);
total->SetParLimits(6, 0., 0.3);
total->FixParameter(7, 0.0657);
total->SetParLimits(9, 9.85, 9.90);
total->SetParLimits(10, 0.05, 0.3);
total->FixParameter(11, 0.0657);
total->SetParLimits(13, 10.28, 10.38);
total->SetParLimits(14, 0.05, 0.3);
total->FixParameter(15, 0.0657);
total->SetParLimits(17, 10.6, 10.7);
total->SetParLimits(18, 0.01, 0.3);
total->FixParameter(19, 0.0657);
total->SetParLimits(21, 11.2, 11.3);
total->SetParLimits(22, 0.01, 0.18);
total->FixParameter(23, 0.0657);
total->SetParLimits(25, 11.45, 11.50);
total->SetParLimits(26, 0.05, 0.15);
total->FixParameter(27, 0.0657);
total->SetParLimits(29, 11.668, 11.672);
total->SetParLimits(30, 0.178, 0.182);
total->FixParameter(31, 0.0657);
total->SetParLimits(33, 11.838, 11.842);
total->SetParLimits(34, 0.237, 0.243);
total->FixParameter(35, 0.0657);
total->SetParLimits(37, 12.046, 12.056);
total->SetParLimits(38, 0.247, 0.255);
total->FixParameter(39, 0.0657);
total->SetParLimits(41, 12.48, 12.56);
total->SetParLimits(42, 0.590, 0.616);
total->FixParameter(43, 0.0657);
total->SetParLimits(45, 12.94, 12.96);
total->SetParLimits(46, 0.287, 0.311);
total->FixParameter(47, 0.0657);
total->SetParLimits(49, 13.298, 13.302);
total->SetParLimits(50, 0.293, 0.331);
total->FixParameter(51, 0.0657);
total->SetParLimits(53, 13.585, 13.601);
total->SetParLimits(54, 0.222, 0.228);
total->FixParameter(55, 0.0657);
total->SetParLimits(57, 13.793, 13.807);
total->SetParLimits(58, 0.198, 0.212);
total->FixParameter(59, 0.0657);
total->SetParLimits(61, 13.997, 14.003);
total->SetParLimits(62, 0.256, 0.274);
total->FixParameter(63, 0.0657);
total->SetParLimits(65, 14.1, 14.3);
total->SetParLimits(66, 0.177, 0.183);
total->FixParameter(67, 0.0657);



// Ajustar solo la función total y desactivar la visualización de la línea de ajuste resultante
exTotalH->Fit(total, "R+");

// Configurar opciones de visualización para la línea de ajuste resultante
total->SetLineColor(kRed);
total->SetLineWidth(2);
total->SetNpx(10000);

exTotalH->SetTitle("Energy distribution ^{7}Li and #alpha"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");

// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
//h->SetFillColor(kGreen);
//h->Draw("HIST SAME");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[68];
total->GetParameters(par_total);

    // Crear el array de funciones Breit-Wigner
    TF1 *new_bw[numFunctions];

// Crear las funciones en el bucle
    for (int i = 0; i < numFunctions; i++) {
        // Usar lambdas para pasar el índice `i` y crear cada TF1 dinámicamente
        new_bw[i] = new TF1(Form("bw%d", i + 1), [i](double *x, double *par) {
            return BW(x, par, i);
        }, 8.5, 14.5, 4); // Intervalo de energía y 4 parámetros

        // Configurar parámetros iniciales para cada Breit-Wigner
      new_bw[i]->SetParameters(par_total[i * 4], par_total[i * 4 + 1], par_total[i * 4 + 2], par_total[i * 4 + 3]);    
      new_bw[i]->SetLineColor(kBlack);
      new_bw[i]->Draw("SAME");
      new_bw[i]->SetNpx(10000);
      }

// Mostrar el Canvas
gPad->Update();

// Configurar la leyenda
TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
legend->AddEntry(new_bw[0], "BW individual fits"); // Agregar entrada para la primera función BW ajustada
legend->AddEntry(total, "Total Fit"); 
//legend->AddEntry(h, "Phase Space"); 
legend->SetBorderSize(0); // Sin borde
legend->Draw(); // Dibujar la leyenda


/*
// Calcular la integral de cada función Breit-Wigner ajustada
Double_t integral_bw1 = new_bw1->Integral(8.5, 12);
Double_t integral_bw2 = new_bw2->Integral(8.5, 12);
Double_t integral_bw3 = new_bw3->Integral(9, 14.0);
Double_t integral_bw4 = new_bw4->Integral(9, 14);
Double_t integral_bw5 = new_bw5->Integral(9, 14);
Double_t integral_bw6 = new_bw6->Integral(9, 14);
Double_t integral_bw7 = new_bw7->Integral(9, 14);
Double_t integral_bw8 = new_bw8->Integral(9, 14);
Double_t integral_bw9 = new_bw9->Integral(9, 14);
Double_t integral_bw10 = new_bw10->Integral(10, 14.5);
Double_t integral_bw11 = new_bw11->Integral(10, 14.5);
Double_t integral_bw12 = new_bw12->Integral(10, 14.5);


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

*/

}

