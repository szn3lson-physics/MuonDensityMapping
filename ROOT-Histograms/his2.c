#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraph.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>

void his2() {

    // --- Styl wykresu ---
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetNumberContours(600);

    const int RES = 1000;       // rozdzielczość siatki
    const double Rmax = 1;       // promień obszaru rysowania

    TCanvas *c1 = new TCanvas("c1","PSF Test",1200,1200);

    TH2F *hist = new TH2F("hist",
                          "Coincidance 4;X;Y",
                          RES, -Rmax, Rmax,
                          RES, -Rmax, Rmax);
    hist->SetStats(0);

    // --- Wczytujemy kąty (0..180°) z filtrowaniem błędnych wartości ---
    std::vector<double> angles;
    std::ifstream fin("Output/test.txt");
    double ang;
    while (fin >> ang) {
        if (ang < 0.0 || ang > 180.0) continue;
        angles.push_back(ang);
    }
    fin.close();

    if (angles.empty()) angles.push_back(45.0); // fallback

    // --- Parametr: szerokość przedziału równomiernego (stopnie) ---
    double spread_deg = 5;           // <-- ustaw tutaj np. 1.0 lub 9.0

    // Pomocniczna lambda: minimalna różnica kątowa w stopniach (0..360)
    auto minimal_angle_diff_deg = [](double a_deg, double b_deg) {
        double diff = fabs(a_deg - b_deg);
        if (diff > 360.0) diff = fmod(diff, 360.0);
        if (diff > 180.0) diff = 360.0 - diff;
        return diff;
    };

    // --- Dla każdego kąta generujemy dwa kierunki: theta i theta + 180 ---
    for (double theta_deg : angles) {

        double thetas_deg[2] = { theta_deg, fmod(theta_deg + 180.0, 360.0) };

        for (int idir = 0; idir < 2; ++idir) {

            double theta_check_deg = thetas_deg[idir];

            // iteracja po siatce punktów
            for (int ix = 0; ix < RES; ++ix) {
                double x = -Rmax + 2.0 * Rmax * ix / (RES - 1);

                for (int iy = 0; iy < RES; ++iy) {
                    double y = -Rmax + 2.0 * Rmax * iy / (RES - 1);

                    double r = sqrt(x*x + y*y);
                    if (r > Rmax) continue; // poza obszarem

                    // wektor kierunku theta (w radianach) do testu "przed detektorem"
                    double theta_rad = theta_check_deg * TMath::DegToRad();
                    double cosT = TMath::Cos(theta_rad);
                    double sinT = TMath::Sin(theta_rad);

                    // Punkt musi leżeć po "stronie źródła" dla tego kierunku
                    double dot = x * cosT + y * sinT;
                    if (dot <= 0) continue; // punkt za detektorem dla tego kierunku

                    // oblicz φ w stopniach w zakresie [0,360)
                    double phi_rad = TMath::ATan2(y, x); // -pi..pi
                    double phi_deg = phi_rad * TMath::RadToDeg();
                    if (phi_deg < 0) phi_deg += 360.0;

                    // minimalna różnica kątowa (stopnie)
                    double ddeg = minimal_angle_diff_deg(phi_deg, theta_check_deg);

                    // rozkład równomierny: jeśli w odległości <= spread_deg
                    if (ddeg <= spread_deg) {
                        hist->Fill(x, y, 1.0);
                    }
                }
            }
        } // koniec pętli idir
    } // koniec pętli po kątach

    // --- NORMALIZACJA OSY Z (SKALI KOLORÓW) do zakresu [0, 1] ---
    
    // 1. Obliczamy maksymalną wartość w histogramie
    double max_content = hist->GetMaximum();
    
    if (max_content > 0) {
        // 2. Skalujemy wszystkie komórki przez 1/max_content
        hist->Scale(1.0 / max_content);
    }
    
    // 3. Ustawiamy skalę kolorów (Z) na [0, 1]
    hist->SetMaximum(1.0);
    // hist->SetMinimum(0.0); // Opcjonalne

    // --- Rysowanie ---
    hist->SetContour(800);
    hist->Draw("COLZ");

    // Detektor: czarna kropka
    TGraph *det = new TGraph();
    det->SetPoint(0, 0.0, 0.0);
    det->SetMarkerStyle(20);
    det->SetMarkerSize(3.0);
    det->SetMarkerColor(kBlack);
    det->Draw("P SAME");
}