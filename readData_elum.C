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
Float_t elum = 0;


void readData_elum()
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
    tree->SetBranchAddress("elum", &elum);


    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("Ex", 1);
    tree->SetBranchStatus("x", 1);
    tree->SetBranchStatus("rdt", 1);
    tree->SetBranchStatus("thetaCM", 1);
    tree->SetBranchStatus("coinTime", 1);
    tree->SetBranchStatus("detID", 1);
    tree->SetBranchStatus("e", 1);
    tree->SetBranchStatus("elum", 1);


    // Definición de cortes de protones

   auto cutProtonRecoil1 = new TCutG("CUTPROTONRECOIL1",7);
   cutProtonRecoil1->SetVarX("rdtH[0]");
   cutProtonRecoil1->SetVarY("");
   cutProtonRecoil1->SetTitle("Graph");
   cutProtonRecoil1->SetFillStyle(1000);
   cutProtonRecoil1->SetPoint(0,3483.32,909.103);
   cutProtonRecoil1->SetPoint(1,3810.2,649.767);
   cutProtonRecoil1->SetPoint(2,3782.18,1281.9);
   cutProtonRecoil1->SetPoint(3,2549.37,1703.32);
   cutProtonRecoil1->SetPoint(4,2624.08,1298.11);
   cutProtonRecoil1->SetPoint(5,3530.02,876.686);
   cutProtonRecoil1->SetPoint(6,3483.32,909.103);
    cutProtonRecoil1->Draw("");


   auto cutProtonRecoil2 = new TCutG("CUTPROTONRECOIL2",7);
   cutProtonRecoil2->SetVarX("rdtH[2]");
   cutProtonRecoil2->SetVarY("");
   cutProtonRecoil2->SetTitle("Graph");
   cutProtonRecoil2->SetFillStyle(1000);
   cutProtonRecoil2->SetPoint(0,2512.01,1508.82);
   cutProtonRecoil2->SetPoint(1,3035.02,973.937);
   cutProtonRecoil2->SetPoint(2,3894.26,682.184);
   cutProtonRecoil2->SetPoint(3,3763.5,1427.77);
   cutProtonRecoil2->SetPoint(4,2521.35,1460.19);
   cutProtonRecoil2->SetPoint(5,2521.35,1460.19);
   cutProtonRecoil2->SetPoint(6,2512.01,1508.82);
    cutProtonRecoil2->Draw("");

    auto  cutProtonRecoil3 = new TCutG("CUTPROTONRECOIL3",7);
   cutProtonRecoil3->SetVarX("rdtH[1]");
   cutProtonRecoil3->SetVarY("");
   cutProtonRecoil3->SetTitle("Graph");
   cutProtonRecoil3->SetFillStyle(1000);
   cutProtonRecoil3->SetPoint(0,2264.51,1687.11);
   cutProtonRecoil3->SetPoint(1,2544.7,1119.81);
   cutProtonRecoil3->SetPoint(2,3628.08,747.018);
   cutProtonRecoil3->SetPoint(3,3721.48,1200.86);
   cutProtonRecoil3->SetPoint(4,3487.99,1379.15);
   cutProtonRecoil3->SetPoint(5,2264.51,1687.11);
   cutProtonRecoil3->SetPoint(6,2264.51,1687.11);
    cutProtonRecoil3->Draw("");

      
   auto cutProtonRecoil4 = new TCutG("cutProtonRecoil4",6);
   cutProtonRecoil4->SetVarX("rdtH[3]");
   cutProtonRecoil4->SetVarY("");
   cutProtonRecoil4->SetTitle("Graph");
   cutProtonRecoil4->SetFillStyle(1000);
   cutProtonRecoil4->SetPoint(0,2311.21,1654.69);
   cutProtonRecoil4->SetPoint(1,2694.13,1168.44);
   cutProtonRecoil4->SetPoint(2,3852.23,779.435);
   cutProtonRecoil4->SetPoint(3,3749.5,1346.73);
   cutProtonRecoil4->SetPoint(4,2301.87,1670.9);
   cutProtonRecoil4->SetPoint(5,2311.21,1654.69);
    cutProtonRecoil4->Draw("");

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

TH1F* elumH = new TH1F("elumH", "elumH", 3, -10, 100);

TH1F* exTotalH = new TH1F("exTotalH","exTotalH",35,10.9,14.4);      


 Long64_t nentries = tree->GetEntries();
 std::cout<<" Number of entries : "<<nentries<<"\n";
   for (Long64_t i=0;i<nentries;i++) {
     tree->GetEntry(i);
    elumH -> Fill(elum);


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

  TCanvas *c5 = new TCanvas();
  elumH->Draw();

// Mostrar el Canvas
gPad->Update();

}

