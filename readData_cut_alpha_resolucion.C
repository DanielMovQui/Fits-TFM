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

     //if (!cutBoronRecoil1->IsInside(rdt[0],rdt[1])  && !cutBoronRecoil2->IsInside(rdt[2],rdt[3]) && !cutBoronRecoil3->IsInside(rdt[4],rdt[5]) && !cutBoronRecoil4->IsInside(rdt[6],rdt[7]))
      //continue; 
    
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

TF1 *g1 = new TF1("m1", "gaus", 5, 9.5);

// Definir parámetros iniciales para cada gaussiana
g1->SetParameters(100, 7.28, 0.1);  // Parámetros iniciales: amplitud, media, desviación estándar

// The total is the sum of the five, each has 3 parameters
TF1 *total = new TF1("mstotal", "gaus(0)", 5, 9.5);

// Fit each function and add it to the list of functions
exTotalH->Fit(g1, "R");

// Get the parameters from the fit
Double_t par[3];
g1->GetParameters(&par[0]);

// Use the parameters on the sum
total->SetParameters(par);

// Configurar nombres de parámetros uno por uno
total->SetParName(0, "A1");
total->SetParName(1, "Mean1");
total->SetParName(2, "Sigma1");

//total->SetParLimits(0, )
total->SetParLimits(1, 7.2, 7.3); 
total->SetParLimits(2, 0.05, 0.200);

// Ajustar solo la función total y desactivar la visualización de la línea de ajuste resultante
exTotalH->Fit(total, "R+");
// Configurar opciones de visualización para la línea de ajuste resultante
total->SetLineColor(kRed);  // Puedes ajustar el color según tus preferencias
total->SetLineWidth(2);     // Puedes ajustar el grosor de la línea según tus preferencias
// Dibujar solo la línea de ajuste resultante
exTotalH->Draw("HIST");
total->Draw("SAME");

// Obtener los parámetros ajustados de la función total
Double_t par_total[3];
total->GetParameters(par_total);

// Crear nuevas funciones Gaussianas con los parámetros ajustados
TF1 *new_g1 = new TF1("new_m1", "gaus", 5, 9.5);

// Establecer los parámetros ajustados en las nuevas funciones
new_g1->SetParameters(par_total);

// Configurar nombres de parámetros
new_g1->SetParName(0, "A1");
new_g1->SetParName(1, "Mean1");
new_g1->SetParName(2, "Sigma1");

// Dibujar las nuevas funciones en el mismo Canvas
new_g1->SetLineColor(kBlack);

new_g1->Draw("SAME");

// Mostrar el Canvas
gPad->Update();


