void readData()
{
TFile *f = new TFile("h082_10BDP_trace_run013_015-019_025-041.root");
TTree *tree = (TTree*)f->Get("tree");

// Asignación de variables flotantes con [] número de datos o espacio.
Float_t Ex = 0;
Float_t e[24];
Float_t rdt[8];


tree->SetBranchAddress("Ex",&Ex);
tree->SetBranchAddress("e", e);
tree->SetBranchAddress("rdt", rdt);

 tree->SetBranchStatus("*", 0); // Desactiva todas las ramas, * es para seleccionar todas las ramas
 // Activo las ramas que necesito 
 tree->SetBranchStatus("Ex", 1);
 tree->SetBranchStatus("x", 1);
 tree->SetBranchStatus("thetaCM", 1);
 tree->SetBranchStatus("coinTime", 1);
 tree->SetBranchStatus("detID", 1);
 tree->SetBranchStatus("rdt", 1);
 tree->SetBranchStatus("e", 1);

//Histograms
// Declaración de un puntero hist(nombre del hist, titulo, bins, minx, maxx)
TH1F* exH = new TH1F("exH","exH",1000,-20,20);
TH1F* eH[24]; // un arreglo de 24 punteros a objetos TH1F
// Crear un objeto tipo TH1F para cada i
for (auto i = 0; i < 24; ++i) // auto fija la variable de iteración, la condición del bucle, ++i incrementa i en 1.
      eH[i] = new TH1F(Form("eH[%i]", i), Form("eH[%i]", i), 1000, -2, 18);

TH2F* rdtt[4];
for (auto i = 0; i < 4; ++i)
    rdtt[i] = new TH2F(Form("rdt[%i]", i), Form("rdt[%i]", i), 1000, 0, 6000, 1000, 0, 6000);

 Long64_t nentries = 10000000; //tree->GetEntries(); // Obtener total de entradas y almacenar en nentries. Datos tipo Long64_t (números enteros largos)
 std::cout<<" Number of entries : "<<nentries<<"\n"; // std::cout para imprimir en C++
   for (Long64_t i=0;i<nentries;i++) { // Bucle para iterar sobre cada una de las entradas
     tree->GetEntry(i);

    if(i % 1000000 == 0){ // Printear cada millón de entradas
         std::cout<<" Entry number : "<<i<<"\n";
         //if(!std::isnan(e[2])) std::cout<<" Energy index 0 "<<e[2]<<"\n"; 
    }    
    exH->Fill(Ex); // Llenar el histograma exH con Ex (sacado de branchAddress e inicializado en 0 antes)

    for(auto j=0;j<24;++j){
        // Llenar los 24 histogramas con la variable e
        eH[j]->Fill(e[j]);
    } 

    for(Int_t j=0;j<4;++j){
        rdtt[j]->Fill(rdt[j * 2],rdt[2 * j+1]);
    }

  }

  // Primer corte 

   TCutG *cutg = new TCutG("CUTG",8);
   cutg->SetVarX("rdt[1]");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,3287.865,2088.75);
   cutg->SetPoint(1,3341.976,2002.75);
   cutg->SetPoint(2,3378.051,1877.75);
   cutg->SetPoint(3,3354.798,1827.75);
   cutg->SetPoint(4,3288.517,1955.75);
   cutg->SetPoint(5,3288.082,2089.75);
   cutg->SetPoint(6,3288.082,2089.75);
   cutg->SetPoint(7,3287.865,2088.75);
   cutg->Draw("l");


// Canvas para plotear y drawear, exH y eH
  TCanvas *c1 = new TCanvas();
  exH->Draw();

  TCanvas *c2 = new TCanvas();
  c2->Divide(5,5);
    for(auto i=0;i<24;++i){
      c2->cd(i+1);
      eH[i]->Draw();
    }

  TCanvas *c3 = new TCanvas();
  c3 -> Divide(2,2);
  for(Int_t i=0;i<4;++i){
    c3->cd(i + 1);
    rdtt[i]->Draw("colz");
  }
}

