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
#include "binContents_proton_10000.C"
#include "penetrabilities_L_0.C"
#include "penetrabilities_L_1.C"
#include "penetrabilities_L_2.C"
#include "penetrabilities_L_3.C"
#include "penetrabilities_L_4.C"
#include "penetrabilities_L_5.C"

Float_t Ex = 0;
Int_t   detID = 0;
Float_t coinTime = 0;
Float_t e[24];
Float_t rdt[8];
Float_t x[24];
Float_t thetaCM = 0;
// NORMALIZATION PARAMETERS
const double solid_angle         = 1.744846665964482;
const double projectiles         = 7.19e12;
const double targets             = 7*10e22*0.00010752688172043011;
const double bin_width           = 0.1591;
const double normalization_factor = solid_angle*projectiles*targets*bin_width;
// L VALUES TAKEN FOR THE FIT
int l_values[6] = {0, 4, 0, 2, 4, 4};

double getPS(double E) {
    double frac      = (E - 11.1591) / (14.5 - 11.1591);
    double idx_exact = frac * 9999.0;
    int    i0        = (int)idx_exact;
    int    i1        = std::min(i0 + 1, 9999);
    double w         = idx_exact - i0;
    if (i0 < 0 || i0 >= 10000) return 0.0;
    return binContents_proton_10000[i0]*(1.0-w)
         + binContents_proton_10000[i1]*w;
}

// ── Breit-Wigner ─────────────────────────────────────────────────────
double BreitWigner(double E, double E0, double Gamma) {
    return Gamma / ((E-E0)*(E-E0) + pow(Gamma/2.0, 2));
}

// ── Gaussian ─────────────────────────────────────────────────────────
double Gaussian(double E, double mean, double sigma) {
    return (1.0/(sigma*sqrt(2.0*M_PI)))
           * exp(-0.5*pow((E-mean)/sigma, 2));
}

// ── Convoluted BW ────────────────────────────────────────────────────
double ConvolutedBW(double *x, double *par, double E_min, double E_max) {
    double E      = x[0];
    double Amp    = par[0];
    double E0     = par[1];
    double Gamma0 = par[2];
    double sigma  = par[3];

    int    n_points = 100;
    double lower    = std::max(E - 1.0*sigma, E_min);
    double upper    = std::min(E + 1.0*sigma, E_max);
    double step     = (upper - lower) / n_points;
    double integral = 0.0;

    for (int i = 0; i < n_points; ++i) {
        double t  = lower + (i + 0.5)*step;
        integral += BreitWigner(t, E0, Gamma0) * Gaussian(E, t, sigma) * step;
    }
    return Amp * integral;
}

// ── Modified BW with penetrability ───────────────────────────────────
double BWModificada(double *x, double *par, int l, int bw_index) {
    static const double E_min_values[6] = {11.1591,11.1591,11.1591,
                                           11.1591,11.1591,11.1591};
    static const double E_max_values[6] = {14.5,14.5,14.5,14.5,14.5,14.5};

    if (bw_index < 0 || bw_index > 5) return 0.0;

    double E      = x[0];
    double Amp    = par[0];
    double E0     = par[1];
    double Gamma0 = par[2];
    double sigma  = par[3];

    const int    num_bins_pen  = 10000;
    double bin_width_pen = (E_max_values[bw_index] - E_min_values[bw_index])
                           / num_bins_pen;

    int bin_index_bw = (int)((E  - E_min_values[bw_index]) / 0.000334);
    if (bin_index_bw < 0 || bin_index_bw >= num_bins_pen) return 0.0;

    int bin_index_E0  = (int)((E0 - E_min_values[bw_index]) / bin_width_pen);
    bin_index_E0 = std::max(0, std::min(bin_index_E0, num_bins_pen-1));

    int bin_index_pen = (int)((E  - E_min_values[bw_index]) / bin_width_pen);
    bin_index_pen = std::min(bin_index_pen, num_bins_pen-1);

    double *T[6] = {T0_values, T1_values, T2_values,
                    T3_values, T4_values, T5_values};

    if (T[l][bin_index_E0] < 1e-10) return 0.0;

    double E_bin = E_min_values[bw_index] + (bin_index_bw + 0.5)*0.000334;

    double Gamma_eff, E_eff;
    if (E_bin >= 11.1591) {
        Gamma_eff    = Gamma0 * T[l][bin_index_pen] / T[l][bin_index_E0];
        double E_val = -Gamma0 * (T[l][bin_index_pen] - T[l][bin_index_E0])
                       / (2.0 * T[l][bin_index_E0]);
        E_eff        = E0 - E_val;
    } else {
        Gamma_eff = Gamma0;
        E_eff     = E0;
    }

    double par_conv[4] = {Amp, E_eff, Gamma_eff, sigma};
    return ConvolutedBW(x, par_conv,
                        E_min_values[bw_index],
                        E_max_values[bw_index]);
}

// ── FWHM calculator ──────────────────────────────────────────────────
double calculateFWHM(TF1 *func, double E0, double E_min, double E_max) {
    double halfMax = func->Eval(E0) / 2.0;
    double E1 = E0, E2 = E0;
    const double step = 0.0001;
    while (E1 > E_min && func->Eval(E1) > halfMax) E1 -= step;
    while (E2 < E_max && func->Eval(E2) > halfMax) E2 += step;
    return E2 - E1;
}

double getPS_smooth(double E, int window = 2000) {
    const double E_min_shifted = 11.1591 - 0.041; 
    const double E_max_shifted = 14.5    - 0.041;  
    
    double frac      = (E - E_min_shifted) / (E_max_shifted - E_min_shifted);
    double idx_exact = frac * 9999.0;
    int    i_center  = (int)idx_exact;
    if (i_center < 0 || i_center >= 10000) return 0.0;

    int    i_min = std::max(0,    i_center - window/2);
    int    i_max = std::min(9999, i_center + window/2);
    double sum   = 0.0;
    for (int i = i_min; i <= i_max; i++)
        sum += binContents_proton_10000[i];
    return sum / (i_max - i_min + 1);
}

// ── Total fit function ───────────────────────────────────────────────
double TotalFunction(double *x, double *par) {
    double E     = x[0];
    double total = 0.0;

    for (int i = 0; i < 6; i++)
        total += BWModificada(&E, par + i*4, l_values[i], i);

    // Phase space with linear interpolation — no staircase noise
    total += par[24] * getPS_smooth(E) * 10e30 / normalization_factor;

    return total;
}



void compare_l_values_11p4() {

    TFile *f    = new TFile("h082_10BDP_trace_run013_015-019_025-041.root");
    TTree *tree = (TTree*)f->Get("tree");

    tree->SetBranchAddress("Ex",      &Ex);
    tree->SetBranchAddress("e",        e);
    tree->SetBranchAddress("rdt",      rdt);
    tree->SetBranchAddress("detID",   &detID);
    tree->SetBranchAddress("coinTime",&coinTime);
    tree->SetBranchAddress("x",       &x);
    tree->SetBranchAddress("thetaCM", &thetaCM);

    tree->SetBranchStatus("*",        0);
    tree->SetBranchStatus("Ex",       1);
    tree->SetBranchStatus("x",        1);
    tree->SetBranchStatus("rdt",      1);
    tree->SetBranchStatus("thetaCM",  1);
    tree->SetBranchStatus("coinTime", 1);
    tree->SetBranchStatus("detID",    1);
    tree->SetBranchStatus("e",        1);

    auto cutProtonRecoil1 = new TCutG("CUTPROTONRECOIL1", 6);
    cutProtonRecoil1->SetPoint(0, 3220.905, 1444.03);
    cutProtonRecoil1->SetPoint(1, 3311.691, 1170.896);
    cutProtonRecoil1->SetPoint(2, 3262.632, 1016.418);
    cutProtonRecoil1->SetPoint(3, 3164.516, 1311.94);
    cutProtonRecoil1->SetPoint(4, 3219.777, 1448.507);
    cutProtonRecoil1->SetPoint(5, 3220.905, 1444.03);

    auto cutProtonRecoil2 = new TCutG("CUTPROTONRECOIL2", 6);
    cutProtonRecoil2->SetPoint(0, 3119.599, 1380.306);
    cutProtonRecoil2->SetPoint(1, 3046.031, 1344.799);
    cutProtonRecoil2->SetPoint(2, 3112.911, 1077.659);
    cutProtonRecoil2->SetPoint(3, 3197.625, 1119.928);
    cutProtonRecoil2->SetPoint(4, 3120.713, 1381.996);
    cutProtonRecoil2->SetPoint(5, 3119.599, 1380.306);

    auto cutProtonRecoil3 = new TCutG("CUTPROTONRECOIL3", 6);
    cutProtonRecoil3->SetPoint(0, 3261.741, 1338.105);
    cutProtonRecoil3->SetPoint(1, 3168.548, 1306.482);
    cutProtonRecoil3->SetPoint(2, 3249.519, 1100.931);
    cutProtonRecoil3->SetPoint(3, 3348.823, 1144.413);
    cutProtonRecoil3->SetPoint(4, 3264.797, 1342.058);
    cutProtonRecoil3->SetPoint(5, 3261.741, 1338.105);

    auto cutProtonRecoil4 = new TCutG("CUTPROTONRECOIL4", 6);
    cutProtonRecoil4->SetPoint(0, 3196.594, 1393.424);
    cutProtonRecoil4->SetPoint(1, 3121.597, 1296.035);
    cutProtonRecoil4->SetPoint(2, 3229.732, 1090.438);
    cutProtonRecoil4->SetPoint(3, 3311.706, 1178.358);
    cutProtonRecoil4->SetPoint(4, 3194.85,  1389.366);
    cutProtonRecoil4->SetPoint(5, 3196.594, 1393.424);

    float pol4[24][5] = {
        {22.2632,  -0.443764,  -95.8043,   1.28447,  80.6196},
        {25.6633,   0.939019,  -86.6471,  -1.27647,  70.5687},
        {39.4859,  -1.42573,  -115.269,    2.50988,  89.4895},
        {26.2832,  -4.68976,   -88.0606,  12.0906,   71.8652},
        {30.0488,   4.98238,   -80.8835,  -8.06761,  59.2426},
        {36.7207,   3.9082,    -95.6975,  -5.3619,   71.3784},
        {20.5705,   3.91639,   -27.3175,  -4.67329,  15.0882},
        {25.4729,   2.61548,   -32.0358,  -0.697981, 17.4898},
        {20.7718,   3.45079,   -30.4761,  -5.86632,  20.0606},
        {47.5046,   2.57078,  -158.759,   -6.83235, 135.341 },
        {40.6197,  -0.42643,   -96.5613,   0.173386, 72.1999},
        {0,0,0,0,0},
        {21.3337,   3.09734,   -31.0887,  -4.00537,  17.6514},
        {23.8342,   0.712722,  -24.1388,   0.185379,  7.58745},
        {18.5734,   2.8507,    -23.9463,  -4.71869,  13.1247},
        {25.1945,   3.23653,   -27.8029,  -4.67166,  12.2644},
        {24.0988,  -7.92557,   -74.3601,  10.7257,   59.0229},
        {30.4273,  -1.36079,  -113.128,    1.84176, 100.789 },
        {19.3952,   4.53279,   -71.1824,  -6.23276,  58.6395},
        {32.0866,  -0.850581, -104.966,    3.29582,  82.7438},
        {26.6236,   6.76763,   -79.3186, -10.5532,   64.7973},
        {26.5402,   4.55678,   -51.5956,  -5.42489,  35.7765},
        {0,0,0,0,0},
        {28.9171,  -1.07817,   -82.8819,   2.867,    74.3833}
    };

    TH1F *exTotalH = new TH1F("exTotalH_shifted","",
                               22, 11.0-0.041, 14.5-0.041);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (i % 1000000 == 0) std::cout << "Entry: " << i << std::endl;

        if (x[detID] < -0.95 || x[detID] > 0.95 ||
            thetaCM < 10 || thetaCM > 45 || e[detID] < 1) continue;

        if (!cutProtonRecoil1->IsInside(rdt[0], rdt[1]) &&
            !cutProtonRecoil2->IsInside(rdt[2], rdt[3]) &&
            !cutProtonRecoil3->IsInside(rdt[4], rdt[5]) &&
            !cutProtonRecoil4->IsInside(rdt[6], rdt[7])) continue;

        double poly_correction =
            pol4[detID][0]
          + pol4[detID][1] * x[detID]
          + pol4[detID][2] * x[detID]*x[detID]
          + pol4[detID][3] * x[detID]*x[detID]*x[detID]
          + pol4[detID][4] * x[detID]*x[detID]*x[detID]*x[detID];

        if (coinTime - poly_correction < -20 ||
            coinTime - poly_correction >  15) continue;

        exTotalH->Fill(Ex - 0.041);
    }

    double binEfficiency10k[] = {
        1, 1, 0.1339, 0.1603, 0.16,   0.1543, 0.1485, 0.1435,
        0.1377, 0.1298, 0.1244, 0.1197, 0.1129, 0.1109, 0.1079,
        0.103,  0.0957, 0.0920, 0.0861, 0.0827, 0.0774, 0.0780, 0.0757
    };
    for (int i = 1; i <= 22; i++) {
        double eff    = binEfficiency10k[i-1];
        double c      = exTotalH->GetBinContent(i);
        double e_err  = exTotalH->GetBinError(i);
        double c_real = c     * 10e30 / (normalization_factor * eff);
        double e_stat = e_err * 10e30 / (normalization_factor * eff);
        double e_sist = c_real * 0.1;
        exTotalH->SetBinContent(i, c_real);
        exTotalH->SetBinError(i, sqrt(e_stat*e_stat + e_sist*e_sist));
    }

    // l VALUES FOR THE 1ST STATE
    int l_bw1_test[4] = {0, 1, 2, 3};

    TCanvas *c = new TCanvas("c_lcomp","l comparison E_{x}~11.4 MeV", 1400, 1100);
    c->Divide(2, 2, 0.005, 0.005);

    for (int ipad = 0; ipad < 4; ipad++) {

        int l1 = l_bw1_test[ipad];  

        c->cd(ipad + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.12);

        TH1F *hh = (TH1F*)exTotalH->Clone(Form("hh_pad%d", ipad));
        hh->SetTitle("");
        hh->SetXTitle("E_{x} (MeV)");
        hh->SetYTitle("d#sigma/d#Omega#upointdE [#mub/sr#upointMeV]");
        hh->GetXaxis()->CenterTitle(true);
        hh->GetYaxis()->CenterTitle(true);
        hh->GetXaxis()->SetTitleFont(62);
        hh->GetYaxis()->SetTitleFont(62);
        hh->SetStats(0);
        hh->GetXaxis()->SetRangeUser(11.2, 14.3);
        hh->GetYaxis()->SetRangeUser(0, 9.5);

        TF1 *bw1 = new TF1(Form("bw1_p%d",ipad),[l1](double *x,double *p){return BWModificada(x,p,l1,0);},11.1591,14.5,4);
        TF1 *bw2 = new TF1(Form("bw2_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[1], 1);},11.1591,14.5,4);
        TF1 *bw3 = new TF1(Form("bw3_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[2], 2);},11.1591,14.5,4);
        TF1 *bw4 = new TF1(Form("bw4_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[3], 3);},11.1591,14.5,4);
        TF1 *bw5 = new TF1(Form("bw5_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[4], 4);},11.1591,14.5,4);
        TF1 *bw6 = new TF1(Form("bw6_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[5], 5);},11.1591,14.5,4);

        bw1->SetParameters(50, 11.45, 0.15, 0.0657);
        bw2->SetParameters(50, 12.1,  0.4,  0.0657);
        bw3->SetParameters(50, 12.5,  0.5,  0.0657);
        bw4->SetParameters(50, 13.3,  0.2,  0.0657);
        bw5->SetParameters(50, 13.5,  0.2,  0.0657);
        bw6->SetParameters(50, 14.0,  0.2,  0.0657);

        hh->Fit(bw1, "RQ0"); hh->Fit(bw2, "RQ0+");
        hh->Fit(bw3, "RQ0+"); hh->Fit(bw4, "RQ0+");
        hh->Fit(bw5, "RQ0+"); hh->Fit(bw6, "RQ0+");

        TF1 *total = new TF1(Form("total_p%d", ipad),
            [l1](double *x, double *par) -> double {
                double val = 0.0;
                val += BWModificada(x, par +  0, l1, 0);
                val += BWModificada(x, par +  4, l_values[1],  1);
                val += BWModificada(x, par +  8, l_values[2],  2);
                val += BWModificada(x, par + 12, l_values[3],  3);
                val += BWModificada(x, par + 16, l_values[4],  4);
                val += BWModificada(x, par + 20, l_values[5],  5);
                val += par[24] * getPS_smooth(x[0]) * 10e30 / normalization_factor;
                return val;
            },
            11.1591, 14.3, 25);

        Double_t par[25];
        bw1->GetParameters(&par[0]);  bw2->GetParameters(&par[4]);
        bw3->GetParameters(&par[8]);  bw4->GetParameters(&par[12]);
        bw5->GetParameters(&par[16]); bw6->GetParameters(&par[20]);
        par[24] = 1.0;
        total->SetParameters(par);
/*
        // Límites idénticos a tu código original
        total->SetParLimits(0,  0.,    0.25);
        total->SetParLimits(1,  11.40, 11.405);
        total->SetParLimits(2,  0.015, 0.03);
        total->FixParameter(3,  0.0657);
        total->SetParLimits(4,  0.,    100.);
        total->SetParLimits(5,  12.1,  12.12);
        total->SetParLimits(6,  0.27,  0.29);
        total->FixParameter(7,  0.0657);
        total->SetParLimits(8,  0.07,  100.);
        total->SetParLimits(9,  12.52, 12.60);
        total->SetParLimits(10, 0.27,  0.29);
        total->FixParameter(11, 0.0657);
        total->SetParLimits(12, 0.,    100.);
        total->SetParLimits(13, 13.16, 13.22);
        total->SetParLimits(14, 0.200, 0.280);
        total->FixParameter(15, 0.0657);
        total->SetParLimits(16, 0.,    100.);
        total->SetParLimits(17, 13.53, 13.56);
        total->SetParLimits(18, 0.200, 0.280);
        total->FixParameter(19, 0.0657);
        total->SetParLimits(20, 0.,    100.);
        total->SetParLimits(21, 13.78, 13.82);
        total->SetParLimits(22, 0.170, 0.230);
        total->FixParameter(23, 0.0657);
        total->SetParLimits(24, 0.2,   100000);
*/

        // Límites idénticos a tu código original
        total->FixParameter(0, 0.25);
        total->FixParameter(1,  11.42);
        total->FixParameter(2,  0.030);
        total->FixParameter(3,  0.0657);
        total->SetParLimits(4,  0, 1.0);
        total->FixParameter(5,  12.07);
        total->FixParameter(6,  0.221);
        total->FixParameter(7,  0.0657);
        total->SetParLimits(8,  0, 1.0);
        total->FixParameter(9,  12.52);
        total->FixParameter(10, 0.288);
        total->FixParameter(11, 0.0657);
        total->SetParLimits(12, 0, 1.0);
        total->FixParameter(13, 13.23);
        total->FixParameter(14, 0.230);
        total->FixParameter(15, 0.0657);
        total->SetParLimits(16, 0, 1.0);
        total->FixParameter(17, 13.53);
        total->FixParameter(18, 0.240);
        total->FixParameter(19, 0.0657);
        total->SetParLimits(20, 0, 1.0);
        total->FixParameter(21, 13.8);
        total->FixParameter(22, 0.193);
        total->FixParameter(23, 0.0657);
        total->FixParameter(24, 570);

        total->SetLineColor(kRed);
        total->SetLineWidth(2);
        total->SetNpx(10000);

        hh->Fit(total, "R+");

        TF1 *psFunc = new TF1(Form("ps_p%d", ipad),
            [](double *x, double *par){
                return par[0] * getPS_smooth(x[0]) * 10e30 / normalization_factor;
            },
            11.1591-0.041, 14.5-0.041, 1);
        psFunc->SetParameter(0, total->GetParameter(24));
        psFunc->SetNpx(10000);
        psFunc->SetLineColor(kGreen+2);
        psFunc->SetLineWidth(2);
        psFunc->SetFillColorAlpha(kGreen, 0.4);
        psFunc->SetFillStyle(1001);

        Double_t par_total[25];
        total->GetParameters(par_total);

        TF1 *nb1 = new TF1(Form("nb1_p%d",ipad),[l1](double *x,double *p){return BWModificada(x,p,l1,0);},11.1591,14.5,4);
        TF1 *nb2 = new TF1(Form("nb2_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[1], 1);},11.1591,14.5,4);
        TF1 *nb3 = new TF1(Form("nb3_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[2], 2);},11.1591,14.5,4);
        TF1 *nb4 = new TF1(Form("nb4_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[3], 3);},11.1591,14.5,4);
        TF1 *nb5 = new TF1(Form("nb5_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[4], 4);},11.1591,14.5,4);
        TF1 *nb6 = new TF1(Form("nb6_p%d",ipad),[=](double *x,double *p){return BWModificada(x,p,l_values[5], 5);},11.1591,14.5,4);

        nb1->SetParameters(&par_total[0]);  nb1->SetNpx(10000); nb1->SetLineColor(kBlack); nb1->SetLineWidth(2);
        nb2->SetParameters(&par_total[4]);  nb2->SetNpx(10000); nb2->SetLineColor(kBlack); nb2->SetLineWidth(2);
        nb3->SetParameters(&par_total[8]);  nb3->SetNpx(10000); nb3->SetLineColor(kBlack); nb3->SetLineWidth(2);
        nb4->SetParameters(&par_total[12]); nb4->SetNpx(10000); nb4->SetLineColor(kBlack); nb4->SetLineWidth(2);
        nb5->SetParameters(&par_total[16]); nb5->SetNpx(10000); nb5->SetLineColor(kBlack); nb5->SetLineWidth(2);
        nb6->SetParameters(&par_total[20]); nb6->SetNpx(10000); nb6->SetLineColor(kBlack); nb6->SetLineWidth(2);

        hh->Draw("E1");
        psFunc->Draw("FC SAME");
        total->Draw("SAME");
        for (TF1 *bw : {nb1,nb2,nb3,nb4,nb5,nb6}) bw->Draw("SAME");

        TLegend *leg = new TLegend(0.55, 0.65, 0.89, 0.89);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.038);
        leg->AddEntry(hh,     "Data",        "lep");
        leg->AddEntry(nb1,    "Individual BW fits", "l");
        leg->AddEntry(total,  "Total fit",   "l");
        leg->AddEntry(psFunc, "Phase space", "f");
        leg->Draw();

        TLatex lat;
        lat.SetNDC();
        lat.SetTextFont(62);
        lat.SetTextSize(0.045);
        lat.DrawLatex(0.18, 0.84,
            Form("E_{x}(11.4):  l_{p} = %d", l1));
        lat.SetTextFont(42);
        lat.SetTextSize(0.038);
        lat.DrawLatex(0.18, 0.77,
        Form("#chi^{2}/ndf = %.2f / %d", total->GetChisquare(), total->GetNDF()));

        std::cout << "\n=== Pad " << ipad+1 << "  (l_bw1 = " << l1 << ") ===" << std::endl;
        std::cout << "  E0    = " << par_total[1] << " MeV" << std::endl;
        std::cout << "  Gamma = " << par_total[2] << " MeV" << std::endl;
        std::cout << "\n=== Pad " << ipad+1 << "  (l_bw1 = " << l1 << ") ===" << std::endl;
        std::cout << "  chi2/ndf = " << total->GetChisquare()
          << " / "           << total->GetNDF() << std::endl;

        std::cout << "\n  Contribucion por bin:" << std::endl;
        std::cout << Form("  %-10s  %-10s  %-10s  %-10s  %-10s",
                  "BinCenter","Obs","Pred","Error","chi2_bin") << std::endl;

        double chi2_manual = 0.0;
        for (int ib = 1; ib <= hh->GetNbinsX(); ib++) {
            double x_center = hh->GetBinCenter(ib);
            if (x_center < 11.1 || x_center > 14.3) continue;

            double obs  = hh->GetBinContent(ib);
            double err  = hh->GetBinError(ib);
            double pred = total->Eval(x_center);

            if (err <= 0) continue;

            double chi2_bin = pow((obs - pred) / err, 2);
            chi2_manual += chi2_bin;
            std::cout << Form("  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f", x_center, obs, pred, err, chi2_bin) << std::endl;
        }
        std::cout << "  chi2 manual total = " << chi2_manual << std::endl;
    }

    c->Update();

}
