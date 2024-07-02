Float_t Ex = 0;
Int_t detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x = 0;
Float_t thetaCM = 0;

void readData_cut_alpha()
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

TF1 *bw1 = new TF1("bw1", "[0] / ((x * x - [1] * [1]) * (x * x - [1] * [1]) + [1] * [1] * [2] * [2])", 8.5, 9.5);
TF1 *bw2 = new TF1("bw2", "[3] / ((x * x - [4] * [4]) * (x * x - [4] * [4]) + [4] * [4] * [5] * [5])", 8.5, 11);
TF1 *bw3 = new TF1("bw3", "[6] / ((x * x - [7] * [7]) * (x * x - [7] * [7]) + [7] * [7] * [8] * [8])", 9, 12);
TF1 *bw4 = new TF1("bw4", "[9] / ((x * x - [10] * [10]) * (x * x - [10] * [10]) + [10] * [10] * [11] * [11])", 9, 13);
TF1 *bw5 = new TF1("bw5", "[12] / ((x * x - [13] * [13]) * (x * x - [13] * [13]) + [13] * [13] * [14] * [14])", 9, 13);
TF1 *bw6 = new TF1("bw6", "[15] / ((x * x - [16] * [16]) * (x * x - [16] * [16]) + [16] * [16] * [17] * [17])", 9, 13);
TF1 *bw7 = new TF1("bw7", "[18]/ ((x * x - [19] * [19]) * (x * x - [19] * [19]) + [19] * [19] * [20] * [20])", 9, 13);
TF1 *bw8 = new TF1("bw8", "[21]/ ((x * x - [22] * [22]) * (x * x - [22] * [22]) + [22] * [22] * [23] * [23])", 9, 13);
TF1 *bw9 = new TF1("bw9", "[24]/ ((x * x - [25] * [25]) * (x * x - [25] * [25]) + [25] * [25] * [26] * [26])", 9, 14);
TF1 *bw10 = new TF1("bw10", "[27]/ ((x * x - [28] * [28]) * (x * x - [28] * [28]) + [28] * [28] * [29] * [29])", 10, 14.5);
TF1 *bw11 = new TF1("bw11", "[30]/ ((x * x - [31] * [31]) * (x * x - [31] * [31]) + [31] * [31] * [32] * [32])", 10, 14.5);
TF1 *bw12 = new TF1("bw12", "[33]/ ((x * x - [34] * [34]) * (x * x - [34] * [34]) + [34] * [34] * [35] * [35])", 10, 14.5);

// Definir parámetros iniciales para cada Breit-Wigner
bw1->SetParameters(300, 9.18, 0.004);  // Parámetros iniciales: amplitud, media, anchura (ancho a media altura)
bw2->SetParameters(100, 9.27, 0.004);
bw3->SetParameters(300, 9.85, 0.10);
bw4->SetParameters(300, 10.28, 0.10);
bw5->SetParameters(300, 10.6, 0.2);
bw6->SetParameters(300, 11.2, 0.2);
bw7->SetParameters(300, 11.45, 0.2);
bw8->SetParameters(300, 11.66, 0.2);
bw9->SetParameters(300, 11.9, 0.2);
bw10->SetParameters(300, 12.5, 0.2);
bw11->SetParameters(300, 13.25, 0.2);
bw12->SetParameters(300, 13.61, 0.2);
// Definir total dada por la suma de los 8 picos
TF1 *total = new TF1("mstotal", "[0] / ((x * x - [1] * [1]) * (x * x - [1] * [1]) + [1] * [1] * [2] * [2]) + [3] / ((x * x - [4] * [4]) * (x * x - [4] * [4]) + [4] * [4] * [5] * [5]) + [6] / ((x * x - [7] * [7]) * (x * x - [7] * [7]) + [7] * [7] * [8] * [8]) + [9] / ((x * x - [10] * [10]) * (x * x - [10] * [10]) + [10] * [10] * [11] * [11]) + [12] / ((x * x - [13] * [13]) * (x * x - [13] * [13]) + [13] * [13] * [14] * [14]) + [15] / ((x * x - [16] * [16]) * (x * x - [16] * [16]) + [16] * [16] * [17] * [17])+[18]/ ((x * x - [19] * [19]) * (x * x - [19] * [19]) + [19] * [19] * [20] * [20])+[21]/ ((x * x - [22] * [22]) * (x * x - [22] * [22]) + [22] * [22] * [23] * [23])+[24]/ ((x * x - [25] * [25]) * (x * x - [25] * [25]) + [25] * [25] * [26] * [26])+[27]/ ((x * x - [28] * [28]) * (x * x - [28] * [28]) + [28] * [28] * [29] * [29])+[30]/ ((x * x - [31] * [31]) * (x * x - [31] * [31]) + [31] * [31] * [32] * [32])+[33]/ ((x * x - [34] * [34]) * (x * x - [34] * [34]) + [34] * [34] * [35] * [35])", 8.5, 14.5);

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
// Obtener los parámetros del fit
Double_t par[36];
bw1->GetParameters(&par[0]);
bw2->GetParameters(&par[3]);
bw3->GetParameters(&par[6]);
bw4->GetParameters(&par[9]);
bw5->GetParameters(&par[12]);
bw6->GetParameters(&par[15]);
bw7->GetParameters(&par[18]);
bw8->GetParameters(&par[21]);
bw9->GetParameters(&par[24]);
bw10->GetParameters(&par[27]);
bw11->GetParameters(&par[30]);
bw12->GetParameters(&par[33]);
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
total->SetParName(18, "Amp7");
total->SetParName(19, "Mean7");
total->SetParName(20, "Width7");
total->SetParName(21, "Amp8");
total->SetParName(22, "Mean8");
total->SetParName(23, "Width8");
total->SetParName(24, "Amp9");
total->SetParName(25, "Mean9");
total->SetParName(26, "Width9");
total->SetParName(27, "Amp10");
total->SetParName(28, "Mean10");
total->SetParName(29, "Width10");
total->SetParName(30, "Amp11");
total->SetParName(31, "Mean11");
total->SetParName(32, "Width11");
total->SetParName(33, "Amp12");
total->SetParName(34, "Mean12");
total->SetParName(35, "Width12");
// Establecer límites de parámetros (si es necesario)
total->SetParLimits(1, 9.17, 9.19);
total->SetParLimits(2, 0.0, 0.3);
total->SetParLimits(4, 9.25, 9.30);
total->SetParLimits(5, 0, 0.3);
total->SetParLimits(7, 9.85, 9.90);
total->SetParLimits(8, 0.05, 0.3);
total->SetParLimits(10, 10.28, 10.38);
total->SetParLimits(11, 0.05, 0.3);
total->SetParLimits(13, 10.6, 10.7);
total->SetParLimits(14, 0.01, 0.3);
total->SetParLimits(16, 11.2, 11.3);
total->SetParLimits(17, 0.01, 0.18);
total->SetParLimits(19, 11.45, 11.50);
total->SetParLimits(20, 0.05, 0.15);
total->SetParLimits(22, 11.66, 11.67);
total->SetParLimits(23, 0.153, 0.155);
total->SetParLimits(25, 11.90, 11.91);
total->SetParLimits(26, 0.282, 0.284);
total->SetParLimits(28, 12.5, 12.51);
total->SetParLimits(29, 0.576, 0.577);
total->SetParLimits(31, 13.25, 13.26);
total->SetParLimits(32, 0.39, 0.4);
total->SetParLimits(34, 13.61, 13.62);
total->SetParLimits(35, 0.34, 0.35);

// Ajustar solo la función total y desactivar la visualización de la línea de ajuste resultante
exTotalH->Fit(total, "R+");

// Configurar opciones de visualización para la línea de ajuste resultante
total->SetLineColor(kRed);
total->SetLineWidth(2);
total->SetNpx(5000);
exTotalH->SetTitle("Energy distribution ^{7}Li and #alpha"); 
exTotalH->SetXTitle("E (MeV)");
exTotalH->SetYTitle("Counts");

// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[36];
total->GetParameters(par_total);

// Crear nuevas funciones Breit-Wigner con los parámetros ajustados
TF1 *new_bw1 = new TF1("new_bw1", "[0] / ((x * x - [1] * [1]) * (x * x - [1] * [1]) + [1] * [1] * [2] * [2])", 8.5, 12);
TF1 *new_bw2 = new TF1("new_bw2", "[3] / ((x * x - [4] * [4]) * (x * x - [4] * [4]) + [4] * [4] * [5] * [5])", 8.5, 12);
TF1 *new_bw3 = new TF1("new_bw3", "[6] / ((x * x - [7] * [7]) * (x * x - [7] * [7]) + [7] * [7] * [8] * [8])", 9, 13.0);
TF1 *new_bw4 = new TF1("new_bw4", "[9] / ((x * x - [10] * [10]) * (x * x - [10] * [10]) + [10] * [10] * [11] * [11])", 9, 13);
TF1 *new_bw5 = new TF1("new_bw5", "[12] / ((x * x - [13] * [13]) * (x * x - [13] * [13]) + [13] * [13] * [14] * [14])", 9, 13);
TF1 *new_bw6 = new TF1("new_bw6", "[15] / ((x * x - [16] * [16]) * (x * x - [16] * [16]) + [16] * [16] * [17] * [17])", 9, 13);
TF1 *new_bw7 = new TF1("new_bw7", "[18] / ((x * x - [19] * [19]) * (x * x - [19] * [19]) + [19] * [19] * [20] * [20])", 9, 13);
TF1 *new_bw8 = new TF1("new_bw8", "[21] / ((x * x - [22] * [22]) * (x * x - [22] * [22]) + [22] * [22] * [23] * [23])", 9, 13);
TF1 *new_bw9 = new TF1("new_bw9", "[24] / ((x * x - [25] * [25]) * (x * x - [25] * [25]) + [25] * [25] * [26] * [26])", 9, 13);
TF1 *new_bw10 = new TF1("new_bw10", "[27] / ((x * x - [28] * [28]) * (x * x - [28] * [28]) + [28] * [28] * [29] * [29])", 10, 14.5);
TF1 *new_bw11 = new TF1("new_bw11", "[30] / ((x * x - [31] * [31]) * (x * x - [31] * [31]) + [31] * [31] * [32] * [32])", 10, 14.5);
TF1 *new_bw12 = new TF1("new_bw12", "[33] / ((x * x - [34] * [34]) * (x * x - [34] * [34]) + [34] * [34] * [35] * [35])", 10, 14.5);

new_bw1->SetNpx(5000);
new_bw2->SetNpx(5000);
new_bw3->SetNpx(5000);
new_bw4->SetNpx(5000);
new_bw5->SetNpx(5000);
new_bw6->SetNpx(5000);
new_bw7->SetNpx(5000);
new_bw8->SetNpx(5000);
new_bw9->SetNpx(5000);
new_bw10->SetNpx(5000);
new_bw11->SetNpx(5000);
new_bw12->SetNpx(5000);

new_bw1->SetParameters(par_total);
new_bw2->SetParameters(par_total);
new_bw3->SetParameters(par_total);
new_bw4->SetParameters(par_total);
new_bw5->SetParameters(par_total);
new_bw6->SetParameters(par_total);
new_bw7->SetParameters(par_total);
new_bw8->SetParameters(par_total);
new_bw9->SetParameters(par_total);
new_bw10->SetParameters(par_total);
new_bw11->SetParameters(par_total);
new_bw12->SetParameters(par_total);
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

new_bw7->SetParName(18, "Amp7");
new_bw7->SetParName(19, "Mean7");
new_bw7->SetParName(20, "Width7");

new_bw8->SetParName(21, "Amp8");
new_bw8->SetParName(22, "Mean8");
new_bw8->SetParName(23, "Width8");

new_bw9->SetParName(24, "Amp9");
new_bw9->SetParName(25, "Mean9");
new_bw9->SetParName(26, "Width9");

new_bw10->SetParName(27, "Amp10");
new_bw10->SetParName(28, "Mean10");
new_bw10->SetParName(29, "Width10");

new_bw11->SetParName(30, "Amp11");
new_bw11->SetParName(31, "Mean11");
new_bw11->SetParName(32, "Width11");

new_bw12->SetParName(33, "Amp12");
new_bw12->SetParName(34, "Mean12");
new_bw12->SetParName(35, "Width12");

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

// Mostrar el Canvas
gPad->Update();

// Configurar la leyenda
TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) donde x1,y1 son esquina inferior izquierda y x2,y2 son esquina superior derecha en fracción del canvas
legend->AddEntry(exTotalH, "Experimental Data"); // Agregar entrada para los datos experimentales
legend->AddEntry(new_bw1, "BW individual fits"); // Agregar entrada para la primera función BW ajustada
legend->AddEntry(total, "Total Fit"); // Agregar entrada para el ajuste total
legend->SetBorderSize(0); // Sin borde
legend->Draw(); // Dibujar la leyenda

// Calcular la integral y la incertidumbre de cada función Breit-Wigner ajustada
Double_t integral_bw1, error_bw1;
Double_t integral_bw2, error_bw2;
Double_t integral_bw3, error_bw3;
Double_t integral_bw4, error_bw4;
Double_t integral_bw5, error_bw5;
Double_t integral_bw6, error_bw6;
Double_t integral_bw7, error_bw7;
Double_t integral_bw8, error_bw8;
Double_t integral_bw9, error_bw9;
Double_t integral_bw10, error_bw10;
Double_t integral_bw11, error_bw11;
Double_t integral_bw12, error_bw12;

integral_bw1 = new_bw1->IntegralAndError(8.5, 12, error_bw1);
integral_bw2 = new_bw2->IntegralAndError(8.5, 12, error_bw2);
integral_bw3 = new_bw3->IntegralAndError(9, 13.0, error_bw3);
integral_bw4 = new_bw4->IntegralAndError(9, 13, error_bw4);
integral_bw5 = new_bw5->IntegralAndError(9, 13, error_bw5);
integral_bw6 = new_bw6->IntegralAndError(9, 13, error_bw6);
integral_bw7 = new_bw7->IntegralAndError(9, 13, error_bw7);
integral_bw8 = new_bw8->IntegralAndError(9, 13, error_bw8);
integral_bw9 = new_bw9->IntegralAndError(9, 13, error_bw9);
integral_bw10 = new_bw10->IntegralAndError(10, 14.5, error_bw10);
integral_bw11 = new_bw11->IntegralAndError(10, 14.5, error_bw11);
integral_bw12 = new_bw12->IntegralAndError(10, 14.5, error_bw12);

// Imprimir los resultados con incertidumbres
std::cout << "Integral of new_bw1: " << integral_bw1 << " ± " << error_bw1 << std::endl;
std::cout << "Integral of new_bw2: " << integral_bw2 << " ± " << error_bw2 << std::endl;
std::cout << "Integral of new_bw3: " << integral_bw3 << " ± " << error_bw3 << std::endl;
std::cout << "Integral of new_bw4: " << integral_bw4 << " ± " << error_bw4 << std::endl;
std::cout << "Integral of new_bw5: " << integral_bw5 << " ± " << error_bw5 << std::endl;
std::cout << "Integral of new_bw6: " << integral_bw6 << " ± " << error_bw6 << std::endl;
std::cout << "Integral of new_bw7: " << integral_bw7 << " ± " << error_bw7 << std::endl;
std::cout << "Integral of new_bw8: " << integral_bw8 << " ± " << error_bw8 << std::endl;
std::cout << "Integral of new_bw9: " << integral_bw9 << " ± " << error_bw9 << std::endl;
std::cout << "Integral of new_bw10: " << integral_bw10 << " ± " << error_bw10 << std::endl;
std::cout << "Integral of new_bw11: " << integral_bw11 << " ± " << error_bw11 << std::endl;
std::cout << "Integral of new_bw12: " << integral_bw12 << " ± " << error_bw12 << std::endl;
}

