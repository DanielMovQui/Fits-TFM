#include <iostream>
#include <vector>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>

using namespace std;

// Constantes
const double hbar = 197.327;    // MeV fm
const double Mn = 939.565;      // MeV
const double Mfrag = 16904.5;   // MeV (ejemplo para 17B)
const double pi = 3.14159265358979323846;

// Función para calcular la masa reducida
double reducedMass(double Mfrag) {
    return (Mfrag * Mn) / (Mfrag + Mn);
}

// Función para calcular el número de onda
double wavenumber(double E, double mu) {
    return sqrt(2 * mu * E) / hbar;
}

// Función para calcular la penetrabilidad
double penetrability(double k, double R, int l) {
    double J_l_half_2 = gsl_sf_bessel_jl(l + 0.5, k * R);
    double Y_l_half_2 = gsl_sf_bessel_yl(l + 0.5, k * R);
    return 1.0 / (J_l_half_2 * J_l_half_2 + Y_l_half_2 * Y_l_half_2);
}

// Función para calcular el factor de cambio de energía
double energyShiftFactor(double k, double R, int l) {
    double J_l_half_2 = gsl_sf_bessel_jl(l + 0.5, k * R);
    double Y_l_half_2 = gsl_sf_bessel_yl(l + 0.5, k * R);
    double J_lminus1_half_2 = gsl_sf_bessel_jl(l - 0.5, k * R);
    double Y_lminus1_half_2 = gsl_sf_bessel_yl(l - 0.5, k * R);
    return -l * J_l_half_2 * J_l_half_2 + Y_l_half_2 * Y_l_half_2 - k * R * (J_l_half_2 * J_lminus1_half_2 + Y_l_half_2 * Y_lminus1_half_2);
}

// Sección eficaz diferencial
double differentialCrossSection(double E, double R, int l, double Gamma0, double mu) {
    double k = wavenumber(E, mu);
    double s_l = R * penetrability(k, R, l) * energyShiftFactor(k, R, l);
    return s_l;
}

void Bessel() {
    // Parámetros de entrada
    double R = 5.0;          // Radio nuclear (fm)
    int l = 0;               // Momento angular
    double Gamma0 = 1.0;     // Ancho reducido (MeV)

    // Cálculo de la masa reducida
    double mu = reducedMass(Mfrag);

    // Rango de energías
    double Emin = 0.1;       // MeV
    double Emax = 5.0;       // MeV
    int numPoints = 100;     // Número de puntos en el gráfico
    double dE = (Emax - Emin) / numPoints;

    // Vectores para almacenar los datos
    vector<double> energies(numPoints);
    vector<double> crossSections(numPoints);

    // Calcular la sección eficaz para diferentes energías
    for (int i = 0; i < numPoints; ++i) {
        double E = Emin + i * dE;
        energies[i] = E;
        crossSections[i] = differentialCrossSection(E, R, l, Gamma0, mu);
    }

    // Crear un gráfico en ROOT
    TCanvas *c1 = new TCanvas("c1", "Sección Eficaz vs Energía", 800, 600);
    TGraph *graph = new TGraph(numPoints, &energies[0], &crossSections[0]);

    // Estilo del gráfico
    graph->SetTitle("Sección Eficaz vs Energía; Energía (MeV); Sección eficaz (fm^2)");
    graph->SetMarkerStyle(20);  // Estilo de los puntos
    graph->SetMarkerSize(0.8);
    graph->SetLineWidth(2);
    graph->SetLineColor(kBlue);

    // Dibujar el gráfico
    graph->Draw("APL");  // "APL" = puntos, líneas y ejes
    c1->Update();

    // Esperar interacción del usuario
    c1->WaitPrimitive();

}

