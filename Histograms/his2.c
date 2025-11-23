#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <fstream>
#include <vector>

void his2() {

    // --- Styl wykresu ---
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetNumberContours(600);

    const int RES = 500;        // rozdzielczość siatki (możesz zwiększyć)
    const double Rmax = 1.0;     // promień obszaru rysowania

    TCanvas *c1 = new TCanvas("c1","Angular Gaussian Map",1200,1200);

    TH2F *hist = new TH2F("hist",
                          "Angular Gaussian Map;X;Y",
                          RES, -Rmax, Rmax,
                          RES, -Rmax, Rmax);
    hist->SetStats(0);

    // --- Wczytujemy kąty (0..180°) z filtrowaniem błędnych wartości ---
    std::vector<double> angles;
    std::ifstream fin("Output/test.txt");
    double ang;
    while (fin >> ang) {
        if (ang < 0.0 || ang > 180.0) {
        continue; // ignorujemy kąt
    }
    angles.push_back(ang);
    }
    fin.close();

    if (angles.empty()) angles.push_back(45.0); // fallback, jeśli plik pusty


    // --- Parametry Gaussa ---
    double sigma_deg = 1.0;                      // sigma w stopniach
    double sigma = sigma_deg * TMath::DegToRad(); // sigma w radianach

    auto gaussian = [&](double dtheta_rad){
        return TMath::Exp(-0.5 * (dtheta_rad * dtheta_rad) / (sigma * sigma));
    };

    auto ang_diff = [&](double phi, double theta){
        double d = fabs(phi - theta);
        if (d > TMath::Pi()) d = 2.0 * TMath::Pi() - d;
        return d;
    };

    // --- Dla każdego kąta generujemy dwa wkłady: theta i theta + pi ---
    for (double theta_deg : angles) {

        // weź oba kierunki: theta i przeciwległy theta+180°
        double thetas_deg[2] = { theta_deg, theta_deg + 180.0 };

        for (int idir = 0; idir < 2; ++idir) {

            double theta = thetas_deg[idir] * TMath::DegToRad();
            double cosT = TMath::Cos(theta);
            double sinT = TMath::Sin(theta);

            // iteracja po siatce punktów
            for (int ix = 0; ix < RES; ++ix) {
                double x = -Rmax + 2.0 * Rmax * ix / (RES - 1);

                for (int iy = 0; iy < RES; ++iy) {
                    double y = -Rmax + 2.0 * Rmax * iy / (RES - 1);

                    double r = sqrt(x*x + y*y);
                    if (r > Rmax) continue; // poza obszarem

                    // tylko półprzestrzeń "przed" detektorem dla danego kierunku:
                    // iloczyn skalarny z wektorem kierunku > 0 oznacza, że punkt jest po stronie źródła
                    double dot = x * cosT + y * sinT;
                    if (dot <= 0) continue; // punkt za detektorem dla tego kierunku

                    double phi = TMath::ATan2(y, x);
                    double dphi = ang_diff(phi, theta);

                    // Gauss kątowy (mocniejszy dla mniejszych dphi)
                    double w_ang = gaussian(dphi);

                    // dodatkowe przyciemnienie/rozjaśnienie w zależności od r:
                    // chcemy, żeby intensywność rosła w stronę źródła (dalej od detektora),
                    // ale można to łatwo zmienić (np. r^1 lub r^2)
                    double radial_factor = r;            // liniowo rośnie z r (możesz zamienić na r*r itp.)
                    
                    // końcowa siła wpisu do histogramu
                    double strength = w_ang * radial_factor * 250.0; // skala do ładnego kontrastu

                    hist->Fill(x, y, strength);
                }
            }
        } // koniec pętli po idir (theta i theta+180)
    } // koniec pętli po kątach

    // --- Rysowanie ---
    hist->SetContour(800);
    hist->Draw("COLZ");

    // Detektor: czarna kropka w (0,0)
    TGraph *det = new TGraph();
    det->SetPoint(0, 0.0, 0.0);
    det->SetMarkerStyle(20);
    det->SetMarkerSize(3.0);
    det->SetMarkerColor(kBlack);
    det->Draw("P SAME");

    // Opcjonalnie: zapisz obraz
    // c1->SaveAs("angular_gauss_map.png");
}
