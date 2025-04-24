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

#include "penetrabilities_neutron_L_0.C"
#include "penetrabilities_neutron_L_1.C"
#include "penetrabilities_neutron_L_2.C"
#include "penetrabilities_neutron_L_3.C"
#include "penetrabilities_neutron_L_4.C"
#include "penetrabilities_neutron_L_5.C"

/*
#include "penetrabilities_L_0.C"
#include "penetrabilities_L_1.C"
#include "penetrabilities_L_2.C"
#include "penetrabilities_L_3.C"
#include "penetrabilities_L_4.C"
#include "penetrabilities_L_5.C"
*/
// BWModificada function to implement penetrability
double BWModificada(double *x, double *par, int l, int bw_index) {
    
    static const double E_min_values[1] = {11.46};
    static const double E_max_values[1] = {14.3};
/*
    static const double E_min_values[1] = {11.23};
    static const double E_max_values[1] = {14.5};    
*/
    double E = x[0];
    double Amp = 10;
    double E0 = 11.6;
    double Gamma0 = 0.2;
    double sigma = 0.0657;

    int num_bins_bw = 10000;
    int num_bins_pen = 10000;
    double bin_width_pen = (E_max_values[bw_index] - E_min_values[bw_index]) / num_bins_pen;

    int bin_index_bw = (E - E_min_values[bw_index]) / 0.000284;
    if (bin_index_bw < 0 || bin_index_bw >= num_bins_bw) return 0.0;

    int bin_index_E0 = static_cast<int>((E0 - E_min_values[bw_index]) / bin_width_pen);
    bin_index_E0 = std::max(0, std::min(bin_index_E0, num_bins_pen - 1));

    int bin_index_pen = (E - E_min_values[bw_index]) / bin_width_pen;
    bin_index_pen = std::min(bin_index_pen, num_bins_pen - 1);

    double E_bin = E_min_values[bw_index] + (bin_index_bw + 0.5) * 0.000284;
    double Gamma_eff;

    if (E_bin >= 11.46) {
        double Gamma_pen_l0 = Gamma0 * T0_neutron_values[bin_index_pen] / T0_neutron_values[bin_index_E0];
        double Gamma_pen_l1 = Gamma0 * T1_neutron_values[bin_index_pen] / T1_neutron_values[bin_index_E0];
        double Gamma_pen_l2 = Gamma0 * T2_neutron_values[bin_index_pen] / T2_neutron_values[bin_index_E0];
        double Gamma_pen_l3 = Gamma0 * T3_neutron_values[bin_index_pen] / T3_neutron_values[bin_index_E0];
        double Gamma_pen_l4 = Gamma0 * T4_neutron_values[bin_index_pen] / T4_neutron_values[bin_index_E0];
        double Gamma_pen_l5 = Gamma0 * T5_neutron_values[bin_index_pen] / T5_neutron_values[bin_index_E0];

/*
    if (E_bin >= 11.23) {
        double Gamma_pen_l0 = Gamma0 * T0_values[bin_index_pen] / T0_values[bin_index_E0];
        double Gamma_pen_l1 = Gamma0 * T1_values[bin_index_pen] / T1_values[bin_index_E0];
        double Gamma_pen_l2 = Gamma0 * T2_values[bin_index_pen] / T2_values[bin_index_E0];            double Gamma_pen_l3 = Gamma0 * T3_values[bin_index_pen] / T3_values[bin_index_E0];
        double Gamma_pen_l4 = Gamma0 * T4_values[bin_index_pen] / T4_values[bin_index_E0];
        double Gamma_pen_l5 = Gamma0 * T5_values[bin_index_pen] / T5_values[bin_index_E0];
      */
        if (l == 0) Gamma_eff = Gamma_pen_l0;
        else if (l == 1) Gamma_eff = Gamma_pen_l1;
        else if (l == 2) Gamma_eff = Gamma_pen_l2;
        else if (l == 3) Gamma_eff = Gamma_pen_l3;
        else if (l == 4) Gamma_eff = Gamma_pen_l4;
        else if (l == 5) Gamma_eff = Gamma_pen_l5;
        else Gamma_eff = Gamma0;
    }

    double bw_modificada = Gamma_eff / ((E - E0)*(E - E0) + (Gamma_eff/2)*(Gamma_eff/2));

    return bw_modificada;
}

// Function to perform Gaussian smearing
double smearOut(double x, TF1 *bw, double sigma) {
    double Amp = 10 ;
    int nPoints = 1000; // Number of points for numerical integration
    double integral = 0.0;
    double step = 2.0 * sigma / nPoints; // Integration step size

    for (int i = 0; i < nPoints; ++i) {
        double e = x - sigma + i * step;
        integral += bw->Eval(e) * TMath::Gaus(x, e, sigma);
    }
    return Amp*integral * step;
}

TF1 *bw = nullptr;
double sigma_global = 0.0657; 

// Wrapper para la Breit-Wigner modificada con L=0
double BWModificada_L0(double *x, double *par) {
    return BWModificada(x, par, 2, 0);
}

// Wrapper para el smear out 
double smearOut_wrapper(double *x, double *p) {
    
    return smearOut(x[0], bw, sigma_global);
}

void Prueba_Convolucion() {

    const double E0 = 11.6;   
    const double Gamma = 0.2;  
    double Amp = 10;
    sigma_global = 0.0657;    

    bw = new TF1("bw", BWModificada_L0, 11, 12.5, 3);
    bw->SetParameters(Amp, E0, Gamma);
    bw->SetLineColor(kRed);
    bw->SetTitle("Neutron Channel BW Convolution;Energia (MeV);Intensity");
    bw->SetNpx(10000);

    TF1 *bw_smeared = new TF1("bw_smeared", smearOut_wrapper, 11, 12.5, 0);
    bw_smeared->SetLineColor(kBlue);
    bw_smeared->SetNpx(1000); // Puedes aumentar para más resolución, pero será más lento

    TCanvas *c = new TCanvas("c", "Breit-Wigner", 800, 600);
    bw->Draw();
    bw_smeared->Draw("SAME");

    auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->AddEntry(bw, "Modified BW ", "l");
    legend->AddEntry(bw_smeared, "BW + Gaussian Smear Out", "l");
    legend->Draw();
}